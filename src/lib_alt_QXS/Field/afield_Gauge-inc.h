/*!
        @file    afield-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$
        @version $LastChangedRevision: 2492 $
*/

#ifndef QXS_AFIELD_GAUGE_INC_INCLUDED
#define QXS_AFIELD_GAUGE_INC_INCLUDED

#include <cstdlib>

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/inline/afield_th-inc.h"

namespace QXS_Gauge {
//====================================================================
  template<typename REALTYPE>
  void set_boundary(AField<REALTYPE, QXS>& ulex,
                    const std::vector<int>& boundary)
  {
    typedef REALTYPE real_t;

#pragma omp barrier

    int Nx   = CommonParameters::Nx();
    int Ny   = CommonParameters::Ny();
    int Nz   = CommonParameters::Nz();
    int Nt   = CommonParameters::Nt();
    int Ndim = CommonParameters::Ndim();
    int Nvol = Nx * Ny * Nz * Nt;

    int Nin = ulex.nin();

    if (!ulex.check_size(Nin, Nvol, Ndim)) {
      vout.crucial("set_boundary: wrong size of input field\n");
      exit(EXIT_FAILURE);
    }

    AIndex_lex<real_t, QXS> index;

    int mu   = 0;
    int ipex = Communicator::ipe(mu);
    int npex = Communicator::npe(mu);

    if ((boundary[mu] != 1) && (ipex == npex - 1)) {
      real_t bc   = real_t(boundary[mu]);
      int    Nyzt = Ny * Nz * Nt;

      int ith, nth, is, ns;
      set_threadtask(ith, nth, is, ns, Nyzt);

      for (int iyzt = is; iyzt < ns; ++iyzt) {
        int site = Nx - 1 + Nx * iyzt;
        for (int in = 0; in < Nin; ++in) {
          int    idx = index.idx_G(in, site, mu);
          real_t uv  = ulex.cmp(idx);
          uv = uv * bc;
          ulex.set(idx, uv);
        }
      }
    }

    mu = 1;
    int ipey = Communicator::ipe(mu);
    int npey = Communicator::npe(mu);

    if ((boundary[mu] != 1) && (ipey == npey - 1)) {
      real_t bc   = real_t(boundary[mu]);
      int    Nxzt = Nx * Nz * Nt;

      int ith, nth, is, ns;
      set_threadtask(ith, nth, is, ns, Nxzt);

      for (int ixzt = is; ixzt < ns; ++ixzt) {
        int ix   = ixzt % Nx;
        int izt  = ixzt / Nx;
        int site = ix + Nx * (Ny - 1 + Ny * izt);
        for (int in = 0; in < Nin; ++in) {
          int    idx = index.idx_G(in, site, mu);
          real_t uv  = ulex.cmp(idx);
          uv = uv * bc;
          ulex.set(idx, uv);
        }
      }
    }

    mu = 2;
    int ipez = Communicator::ipe(mu);
    int npez = Communicator::npe(mu);

    if ((boundary[mu] != 1) && (ipez == npez - 1)) {
      real_t bc   = real_t(boundary[mu]);
      int    Nxy  = Nx * Ny;
      int    Nxyt = Nxy * Nt;

      int ith, nth, is, ns;
      set_threadtask(ith, nth, is, ns, Nxyt);

      for (int ixyt = is; ixyt < ns; ++ixyt) {
        int ixy  = ixyt % Nxy;
        int it   = ixyt / Nxy;
        int site = ixy + Nxy * (Nz - 1 + Nz * it);
        for (int in = 0; in < Nin; ++in) {
          int    idx = index.idx_G(in, site, mu);
          real_t uv  = ulex.cmp(idx);
          uv = uv * bc;
          ulex.set(idx, uv);
        }
      }
    }

    mu = 3;
    int ipet = Communicator::ipe(mu);
    int npet = Communicator::npe(mu);

    if ((boundary[mu] != 1) && (ipet == npet - 1)) {
      real_t bc   = real_t(boundary[mu]);
      int    Nxyz = Nx * Ny * Nz;

      int ith, nth, is, ns;
      set_threadtask(ith, nth, is, ns, Nxyz);

      for (int ixyz = is; ixyz < ns; ++ixyz) {
        int site = ixyz + Nxyz * (Nt - 1);
        for (int in = 0; in < Nin; ++in) {
          int    idx = index.idx_G(in, site, mu);
          real_t uv  = ulex.cmp(idx);
          uv = uv * bc;
          ulex.set(idx, uv);
        }
      }
    }

#pragma omp barrier
  }
} // namespace QXS_Gauge

//============================================================END=====
#endif
