/*!
        @file    fieldIO_Fortran.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2015-03-24 18:19:19 #$

        @version $LastChangedRevision: 2314 $
*/

#include "fieldIO_Fortran.h"

#include <fstream>

#include "bridgeIO.h"
using Bridge::vout;

const std::string FieldIO_Fortran::class_name = "FieldIO_Fortran";

//====================================================================
void FieldIO_Fortran::read_file(Field& v, const std::string& filename)
{
  const int nin_field = v.nin();
  const int nex_field = v.nex();

  int nin_file = m_format->nin();
  int nex_file = m_format->nex();

  if ((nin_file == 0) || (nex_file == 0)) {
    nin_file = nin_field;
    nex_file = nex_field;
  }

  const long_t Lvol = CommonParameters::Lvol();
  const long_t size = nin_file * Lvol * nex_file;

  Field vtmp;

  if (Communicator::is_primary()) {
    vout.detailed(m_vl, "reading field data from %s", filename.c_str());

    vtmp.reset(nin_field, Lvol, nex_field);

    std::fstream infile;
    infile.open(filename.c_str(), std::ios::in | std::ios::binary);
    if (!infile) {
      vout.crucial(m_vl, "Error at %s: file open failure. %s may not exist.\n", class_name.c_str(), filename.c_str());
      Communicator::abort();
    }

    bool do_byteswap = false;

    // Skipping 4B data size area (for size less than 2 GB)
    // infile.seekg(4);

    // read header: record length (bytes)
    // assume header is 4bytes.
    uint32_t flen;
    infile.read((char *)&flen, sizeof(uint32_t) * 1);
    if (!infile) {
      vout.crucial(m_vl, "Error at %s: file read error for header bytes.\n", class_name.c_str());
      Communicator::abort();
    }

    // find byteorder
    if (flen == sizeof(double) * size) {
      // file byteorder is equal to machine byteorder
      vout.detailed(m_vl, "endian matched.\n");
    } else {
      uint32_t flen_s = flen;
      byte_swap(&flen_s, sizeof(uint32_t), 1);

      if (flen_s == sizeof(double) * size) {
        // file byteorder is different from machine byteorder: need swap.
        vout.detailed(m_vl, "different endian. do swap.\n");
        do_byteswap = true;
      } else {
        vout.crucial(m_vl, "Error at %s: size mismatch or format unidentified.\n", class_name.c_str());
        Communicator::abort();
      }
    }

    // read content
    const int block_size = nin_file;
    char      buf[sizeof(double) * block_size];

    for (int j = 0; j < nex_file; ++j) {
      for (long_t isite = 0; isite < Lvol; ++isite) {
        // read 1 block
        infile.read(buf, sizeof(double) * block_size);

        if (!infile) {
          if (infile.eof()) {
            vout.crucial(m_vl, "Error at %s: file size too small.\n", __func__);
          } else {
            vout.crucial(m_vl, "Error at %s: io error.\n", __func__);
          }

          Communicator::abort();
        }

        if (do_byteswap) {
          byte_swap(buf, sizeof(double), block_size);
        }

        double *ptr = (double *)buf;

        for (int i = 0; i < nin_file; ++i) {
          int s, t;
          m_format->file_to_field(s, t, i, j);

          vtmp.set(s, isite, t, ptr[i]);
        }
      }
    }

    // read trailer
    uint32_t ftail;
    infile.read((char *)&ftail, sizeof(uint32_t) * 1);

    if (flen != ftail) {
      vout.crucial(m_vl, "Error at %s: record info mismatch.\n", __func__);
      Communicator::abort();
    }
    if (!infile) {
      vout.crucial(m_vl, "Error at %s: file read failed.\n", __func__);
      Communicator::abort();
    }

    infile.close();
  }

  FieldIO::deliver(&v, &vtmp);

  vout.detailed(m_vl, "read successful\n");
}


//====================================================================
void FieldIO_Fortran::write_file(Field& v, const std::string& filename)
{
  const int nin_field = v.nin();
  const int nex_field = v.nex();

  int nin_file = m_format->nin();
  int nex_file = m_format->nex();

  if ((nin_file == 0) || (nex_file == 0)) {
    nin_file = nin_field;
    nex_file = nex_field;
  }

  const long_t Lvol = CommonParameters::Lvol();

  Field vtmp;
  if (Communicator::is_primary()) {
    vtmp.reset(nin_field, Lvol, nex_field);
  }

  const size_t count = nin_file * Lvol * nex_file;

  FieldIO::gather(&vtmp, &v);

  if (Communicator::is_primary()) {
    vout.detailed(m_vl, "writing field data to %s\n", filename.c_str());

    std::ofstream outfile(filename.c_str(), std::ios::out | std::ios::binary);
    if (!outfile) {
      vout.crucial(m_vl, "Error at %s: file open of %s failed\n", class_name.c_str(), filename.c_str());
      Communicator::abort();
    }

    const uint32_t flen = sizeof(double) * count;

    // record header: record length (bytes)
    outfile.write((char *)&flen, sizeof(uint32_t) * 1);

    if (!outfile) {
      vout.crucial(m_vl, "Error at %s: io error. write record header failed.\n", __func__);
      Communicator::abort();
    }

    // record content
    const int block_size = nin_file;
    char      buf[sizeof(double) * block_size];


    for (int j = 0; j < nex_file; ++j) {
      for (long_t isite = 0; isite < Lvol; ++isite) {
        double *ptr = (double *)buf;

        for (int i = 0; i < nin_file; ++i) {
          int s, t;
          m_format->file_to_field(s, t, i, j);

          ptr[i] = vtmp.cmp(s, isite, t);
        }

        outfile.write(buf, sizeof(double) * block_size);

        if (!outfile) {
          vout.crucial(m_vl, "Error at %s: io error.\n", __func__);
          Communicator::abort();
        }
      }
    }

    // record footer: record length (bytes)
    outfile.write((char *)&flen, sizeof(uint32_t) * 1);

    if (!outfile) {
      vout.crucial(m_vl, "Error at %s: write record footer failed.\n", __func__);
      Communicator::abort();
    }

    outfile.close();
  }

  vout.detailed(m_vl, "write succeeded.\n");
}


//====================================================================
//============================================================END=====
