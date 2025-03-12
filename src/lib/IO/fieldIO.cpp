/*!
        @file    fieldIO.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate: 2015-03-24 18:19:19 #$

        @version $LastChangedRevision: 2492 $
*/

#include "fieldIO.h"

#include "Communicator/communicator.h"
#include "bridgeIO.h"
using Bridge::vout;

#include "Tools/file_utils.h"

//#include <errno.h>
#ifdef HAVE_ENDIAN
#include <endian.h>
#endif

const std::string FieldIO::class_name = "FieldIO";

//====================================================================
void FieldIO::deliver(Field *vlocal, Field *vglobal)
{
  const int nin  = vlocal->nin();
  const int nvol = vlocal->nvol();
  const int nex  = vlocal->nex();

  const int    Lx   = CommonParameters::Lx();
  const int    Ly   = CommonParameters::Ly();
  const int    Lz   = CommonParameters::Lz();
  const int    Lt   = CommonParameters::Lt();
  const long_t Lvol = CommonParameters::Lvol();

  const int Nx   = CommonParameters::Nx();
  const int Ny   = CommonParameters::Ny();
  const int Nz   = CommonParameters::Nz();
  const int Nt   = CommonParameters::Nt();
  const int Nvol = CommonParameters::Nvol();

  const int NPEx = CommonParameters::NPEx();
  const int NPEy = CommonParameters::NPEy();
  const int NPEz = CommonParameters::NPEz();
  const int NPEt = CommonParameters::NPEt();
  const int NPE  = CommonParameters::NPE();

  if (Communicator::is_primary()) {
    if ((nin != vglobal->nin()) ||
        (nex != vglobal->nex()) ||
        (Lvol != vglobal->nvol())) {
      vout.crucial(m_vl, "%s: %s: size mismatch.\n", class_name.c_str(), __func__);
      exit(EXIT_FAILURE);
    }
  }

  Index_lex gindex(Lx, Ly, Lz, Lt);

  Field vtmp(nin, nvol, nex);

  Communicator::sync();

  for (int iblock = 0; iblock < NPE; ++iblock) {
    int ipx = (iblock) % NPEx;
    int ipy = (iblock / NPEx) % NPEy;
    int ipz = (iblock / NPEx / NPEy) % NPEz;
    int ipt = (iblock / NPEx / NPEy / NPEz) % NPEt;

    if (Communicator::is_primary()) {
      for (int j = 0; j < nex; ++j) {
        for (int t = 0; t < Nt; ++t) {
          int t2 = t + Nt * ipt;

          for (int z = 0; z < Nz; ++z) {
            int z2 = z + Nz * ipz;

            for (int y = 0; y < Ny; ++y) {
              int y2 = y + Ny * ipy;

              for (int x = 0; x < Nx; ++x) {
                int x2 = x + Nx * ipx;

                int lsite = idx.site(x, y, z, t);
                int gsite = gindex.site(x2, y2, z2, t2);

                for (int i = 0; i < nin; ++i) {
                  vtmp.set(i, lsite, j, vglobal->cmp(i, gsite, j));
                }
              }
            }
          }
        }
      }
    }

    const int size = nin * nvol * nex;

    int ipe;
    int coord[4] = { ipx, ipy, ipz, ipt, };

    Communicator::grid_rank(&ipe, coord);

    //send_1to1(size, vlocal, &vtmp, ipe, 0, ipe);
    Communicator::send_1to1(size, vlocal->ptr(0), vtmp.ptr(0),
                            ipe, 0, ipe);

    Communicator::sync();
  }

  Communicator::sync();
}


//====================================================================
void FieldIO::gather(Field *vglobal, Field *vlocal)
{
  const int nin  = vlocal->nin();
  const int nvol = vlocal->nvol();
  const int nex  = vlocal->nex();

  const int    Lx   = CommonParameters::Lx();
  const int    Ly   = CommonParameters::Ly();
  const int    Lz   = CommonParameters::Lz();
  const int    Lt   = CommonParameters::Lt();
  const long_t Lvol = CommonParameters::Lvol();

  const int Nx   = CommonParameters::Nx();
  const int Ny   = CommonParameters::Ny();
  const int Nz   = CommonParameters::Nz();
  const int Nt   = CommonParameters::Nt();
  const int Nvol = CommonParameters::Nvol();

  const int NPEx = CommonParameters::NPEx();
  const int NPEy = CommonParameters::NPEy();
  const int NPEz = CommonParameters::NPEz();
  const int NPEt = CommonParameters::NPEt();
  const int NPE  = CommonParameters::NPE();

  if (Communicator::is_primary()) {
    if ((nin != vglobal->nin()) ||
        (nex != vglobal->nex()) ||
        (Lvol != vglobal->nvol())) {
      vout.crucial(m_vl, "%s: %s: size mismatch.\n", class_name.c_str(), __func__);
      exit(EXIT_FAILURE);
    }
  }

  Index_lex gindex(Lx, Ly, Lz, Lt);

  Field vtmp(nin, nvol, nex);

  Communicator::sync();

  for (int iblock = 0; iblock < NPE; ++iblock) {
    int ipx = (iblock) % NPEx;
    int ipy = (iblock / NPEx) % NPEy;
    int ipz = (iblock / NPEx / NPEy) % NPEz;
    int ipt = (iblock / NPEx / NPEy / NPEz) % NPEt;

    int ipe;
    int coord[4] = { ipx, ipy, ipz, ipt, };

    Communicator::grid_rank(&ipe, coord);

    const int size = nin * nvol * nex;

    // send_1to1(size, &vtmp, vlocal, 0, ipe, ipe);
    Communicator::send_1to1(size, vtmp.ptr(0), vlocal->ptr(0),
                            0, ipe, ipe);

    Communicator::sync();

    if (Communicator::is_primary()) {
      for (int j = 0; j < nex; ++j) {
        for (int t = 0; t < Nt; ++t) {
          int t2 = t + Nt * ipt;
          for (int z = 0; z < Nz; ++z) {
            int z2 = z + Nz * ipz;
            for (int y = 0; y < Ny; ++y) {
              int y2 = y + Ny * ipy;
              for (int x = 0; x < Nx; ++x) {
                int x2 = x + Nx * ipx;

                int site  = idx.site(x, y, z, t);
                int gsite = gindex.site(x2, y2, z2, t2);

                for (int i = 0; i < nin; ++i) {
                  vglobal->set(i, gsite, j, vtmp.cmp(i, site, j));
                }
              }
            }
          }
        }
      }
    }
  }

  Communicator::sync();
}


//====================================================================

bool FieldIO::is_bigendian()
{
#if defined(__BYTE_ORDER)
  return __BYTE_ORDER == __BIG_ENDIAN;
#else
  union
  {
    int  l;
    char c[sizeof(int)];
  }
  u;

  u.l = 1;

  return (u.c[sizeof(int) - 1] == 1) ? true : false;
#endif
}


//====================================================================

void FieldIO::convert_endian(void *ptr, size_t size, size_t nmemb)
{
  switch (size)
  {
  case 1: // bytes: do nothing.
    break;

  case 2:
  {         // uint16_t (short)
    uint16_t *p = (uint16_t *)ptr;

    for (unsigned int i = 0; i < nmemb; ++i) {
      uint16_t v = p[i];
      uint16_t w;

      w  = v >> 8 & 0x00ff;
      w |= v << 8 & 0xff00;

      p[i] = w;
    }

    break;
  }

  case 4:
  {         // uint32_t
    uint32_t *p = (uint32_t *)ptr;

    for (unsigned int i = 0; i < nmemb; ++i) {
      uint32_t v = p[i];
      uint32_t w;

      w  = v >> 24 & 0x000000ff;
      w |= v >> 8 & 0x0000ff00;
      w |= v << 8 & 0x00ff0000;
      w |= v << 24 & 0xff000000;

      p[i] = w;
    }

    break;
  }

  case 8:
  {         // uint64_t
    uint32_t *p = (uint32_t *)ptr;

    for (unsigned int i = 0; i < nmemb; ++i) {
      uint32_t v1 = *p;
      uint32_t v2 = *(p + 1);
      uint32_t w1, w2;

      w1  = v1 >> 24 & 0x000000ff;
      w1 |= v1 >> 8 & 0x0000ff00;
      w1 |= v1 << 8 & 0x00ff0000;
      w1 |= v1 << 24 & 0xff000000;

      w2  = v2 >> 24 & 0x000000ff;
      w2 |= v2 >> 8 & 0x0000ff00;
      w2 |= v2 << 8 & 0x00ff0000;
      w2 |= v2 << 24 & 0xff000000;

      *p       = w2;
      *(p + 1) = w1;

      p += 2;
    }

    break;
  }

  default:
//    return EINVAL;
    vout.crucial("%s: %s: unsupported word size.\n", class_name.c_str(), __func__);
    exit(EXIT_FAILURE);
  }

//  return 0;
}


//====================================================================
//============================================================END=====
