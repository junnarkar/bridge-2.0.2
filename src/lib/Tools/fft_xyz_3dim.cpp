/*!
        @file    fft_xyz_3dim.cpp

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2023-04-17 00:20:18 #$

        @version $LastChangedRevision: 2512 $
*/

#ifdef USE_FFTWLIB

#include "fft_xyz_3dim.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = FFT_xyz_3dim::register_factory();
}
#endif

const std::string FFT_xyz_3dim::class_name = "FFT_xyz_3dim";

//====================================================================
void FFT_xyz_3dim::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  string str_fft_direction;

  int err = 0;
  err += params.fetch_string("FFT_direction", str_fft_direction);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(str_fft_direction);
}


//====================================================================
void FFT_xyz_3dim::get_parameters(Parameters& params) const
{
  params.set_string("FFT_direction", m_is_forward ? "Forward" : "Backward");

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void FFT_xyz_3dim::set_parameters(const string& str_fft_direction)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  FFT_direction = %s\n", str_fft_direction.c_str());

  //- range check

  //- store values
  if (str_fft_direction == "Forward") {
    m_is_forward = true;
  } else if (str_fft_direction == "Backward") {
    m_is_forward = false;
  } else {
    vout.crucial(m_vl, "Error at %s: unsupported FFT direction \"%s\"\n", class_name.c_str(), str_fft_direction.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void FFT_xyz_3dim::init()
{
  //- global lattice size
  const int Lx = CommonParameters::Lx();
  const int Ly = CommonParameters::Ly();
  const int Lz = CommonParameters::Lz();

#ifdef USE_OPENMP
  int threads_ok = fftw_init_threads();
#endif

#ifdef USE_MPI
  const int NPE_x = CommonParameters::NPEx();
  const int NPE_y = CommonParameters::NPEy();
  // const int NPE_z = CommonParameters::NPEz();
  const int NPE_t = CommonParameters::NPEt();

  if ((NPE_x * NPE_y * NPE_t) != 1) {
    vout.crucial(m_vl, "Error at %s: FFTW supports parallelization only in z-direction.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }


  fftw_mpi_init();


  //- allocate m_in,out = m_in,out[Nz][Ly][Lx]
  const ptrdiff_t Lx_p = CommonParameters::Lx();
  const ptrdiff_t Ly_p = CommonParameters::Ly();
  const ptrdiff_t Lz_p = CommonParameters::Lz();

  ptrdiff_t fftw_size_p = fftw_mpi_local_size_3d(Lz_p, Ly_p, Lx_p,
                                                 Communicator_impl::world(),
                                                 &m_Nz_p, &m_z_start_p);

  m_in  = fftw_alloc_complex(fftw_size_p);
  m_out = fftw_alloc_complex(fftw_size_p);

  if (!m_in || !m_out) {
    vout.crucial(m_vl, "Error at %s: failed to allocate memory %d [Byte].\n",
                 class_name.c_str(), (int)fftw_size_p);
    exit(EXIT_FAILURE);
  }
#else
  //- allocate m_in,out = m_in,out[Nz][Ly][Lx]
  const size_t fftw_size = sizeof(fftw_complex) * Lx * Ly * Lz;
  m_in  = (fftw_complex *)fftw_malloc(fftw_size);
  m_out = (fftw_complex *)fftw_malloc(fftw_size);

  if (!m_in || !m_out) {
    vout.crucial(m_vl, "Error at %s: failed to allocate memory %d [Byte].\n",
                 class_name.c_str(), (int)fftw_size);
    exit(EXIT_FAILURE);
  }
#endif
}


//====================================================================
void FFT_xyz_3dim::tidy_up()
{
  if (m_in) fftw_free(m_in);
  if (m_out) fftw_free(m_out);
  if (m_plan) fftw_destroy_plan(m_plan);
}


//====================================================================
void FFT_xyz_3dim::fft(Field& field)
{
  //- global lattice size
  const int Lx   = CommonParameters::Lx();
  const int Ly   = CommonParameters::Ly();
  const int Lz   = CommonParameters::Lz();
  const int Lt   = CommonParameters::Lt();
  const int Lxyz = Lx * Ly * Lz;

  //- local size
  const int Nz = CommonParameters::Nz();

  const int Nin = field.nin();
  const int Nex = field.nex();


  //- setup FFTW plan
#ifdef USE_OPENMP
  const int Nthread = ThreadManager::get_num_threads();
  fftw_plan_with_nthreads(Nthread);
#endif
#ifdef USE_MPI
  const ptrdiff_t Lx_p = CommonParameters::Lx();
  const ptrdiff_t Ly_p = CommonParameters::Ly();
  const ptrdiff_t Lz_p = CommonParameters::Lz();

  if (m_is_forward) {
    m_plan = fftw_mpi_plan_dft_3d(Lz_p, Ly_p, Lx_p, m_in, m_out,
                                  Communicator_impl::world(),
                                  FFTW_FORWARD, FFTW_ESTIMATE);
  } else {
    m_plan = fftw_mpi_plan_dft_3d(Lz_p, Ly_p, Lx_p, m_in, m_out,
                                  Communicator_impl::world(),
                                  FFTW_BACKWARD, FFTW_ESTIMATE);
  }
#else
  if (m_is_forward) {
    m_plan = fftw_plan_dft_3d(Lz, Ly, Lx, m_in, m_out,
                              FFTW_FORWARD, FFTW_ESTIMATE);
  } else {
    m_plan = fftw_plan_dft_3d(Lz, Ly, Lx, m_in, m_out,
                              FFTW_BACKWARD, FFTW_ESTIMATE);
  }
#endif


  // ####  Execution main part  ####
  //- Nin is devided by 2, because of complex(i.e. real and imag)
  for (int in2 = 0; in2 < Nin / 2; ++in2) {
    for (int t_global = 0; t_global < Lt; t_global++) {
      for (int ex = 0; ex < Nex; ++ex) {
        //- input data
        for (int z = 0; z < Nz; z++) {
          for (int y_global = 0; y_global < Ly; y_global++) {
            for (int x_global = 0; x_global < Lx; x_global++) {
              int isite_xyz_local = x_global + Lx * (y_global + Ly * z);

              int isite  = m_index.site(x_global, y_global, z, t_global);
              int i_real = 2 * in2;
              int i_imag = 2 * in2 + 1;

              m_in[isite_xyz_local][0] = field.cmp(i_real, isite, ex);
              m_in[isite_xyz_local][1] = field.cmp(i_imag, isite, ex);
            }
          }
        }


        fftw_execute(m_plan);


        //- output data
        for (int z = 0; z < Nz; z++) {
          for (int y_global = 0; y_global < Ly; y_global++) {
            for (int x_global = 0; x_global < Lx; x_global++) {
              int isite_xyz_local = x_global + Lx * (y_global + Ly * z);

              int isite  = m_index.site(x_global, y_global, z, t_global);
              int i_real = 2 * in2;
              int i_imag = 2 * in2 + 1;

              field.set(i_real, isite, ex, m_out[isite_xyz_local][0]);
              field.set(i_imag, isite, ex, m_out[isite_xyz_local][1]);
            }
          }
        }
      }
    }
  }
  //- end of global loops

  //- normailzation for FFTW_BACKWARD
  if (!m_is_forward) {
    scal(field, 1.0 / Lxyz);
  }
}


//====================================================================
void FFT_xyz_3dim::fft(Field& field_out, const Field& field_in)
{
  //- global lattice size
  const int Lx   = CommonParameters::Lx();
  const int Ly   = CommonParameters::Ly();
  const int Lz   = CommonParameters::Lz();
  const int Lt   = CommonParameters::Lt();
  const int Lxyz = Lx * Ly * Lz;

  //- local size
  const int Nz = CommonParameters::Nz();

  const int Nin = field_in.nin();
  const int Nex = field_in.nex();


  //- setup FFTW plan
#ifdef USE_OPENMP
  const int Nthread = ThreadManager::get_num_threads();
  fftw_plan_with_nthreads(Nthread);
#endif
#ifdef USE_MPI
  const ptrdiff_t Lx_p = CommonParameters::Lx();
  const ptrdiff_t Ly_p = CommonParameters::Ly();
  const ptrdiff_t Lz_p = CommonParameters::Lz();

  if (m_is_forward) {
    m_plan = fftw_mpi_plan_dft_3d(Lz_p, Ly_p, Lx_p, m_in, m_out,
                                  Communicator_impl::world(),
                                  FFTW_FORWARD, FFTW_ESTIMATE);
  } else {
    m_plan = fftw_mpi_plan_dft_3d(Lz_p, Ly_p, Lx_p, m_in, m_out,
                                  Communicator_impl::world(),
                                  FFTW_BACKWARD, FFTW_ESTIMATE);
  }
#else
  if (m_is_forward) {
    m_plan = fftw_plan_dft_3d(Lz, Ly, Lx, m_in, m_out,
                              FFTW_FORWARD, FFTW_ESTIMATE);
  } else {
    m_plan = fftw_plan_dft_3d(Lz, Ly, Lx, m_in, m_out,
                              FFTW_BACKWARD, FFTW_ESTIMATE);
  }
#endif


  // ####  Execution main part  ####
  //- Nin is devided by 2, because of complex(i.e. real and imag)
  for (int in2 = 0; in2 < Nin / 2; ++in2) {
    for (int t_global = 0; t_global < Lt; t_global++) {
      for (int ex = 0; ex < Nex; ++ex) {
        //- input data
        for (int z = 0; z < Nz; z++) {
          for (int y_global = 0; y_global < Ly; y_global++) {
            for (int x_global = 0; x_global < Lx; x_global++) {
              int isite_xyz_local = x_global + Lx * (y_global + Ly * z);

              int isite  = m_index.site(x_global, y_global, z, t_global);
              int i_real = 2 * in2;
              int i_imag = 2 * in2 + 1;

              m_in[isite_xyz_local][0] = field_in.cmp(i_real, isite, ex);
              m_in[isite_xyz_local][1] = field_in.cmp(i_imag, isite, ex);
            }
          }
        }


        fftw_execute(m_plan);


        //- output data
        for (int z = 0; z < Nz; z++) {
          for (int y_global = 0; y_global < Ly; y_global++) {
            for (int x_global = 0; x_global < Lx; x_global++) {
              int isite_xyz_local = x_global + Lx * (y_global + Ly * z);

              int isite  = m_index.site(x_global, y_global, z, t_global);
              int i_real = 2 * in2;
              int i_imag = 2 * in2 + 1;

              field_out.set(i_real, isite, ex, m_out[isite_xyz_local][0]);
              field_out.set(i_imag, isite, ex, m_out[isite_xyz_local][1]);
            }
          }
        }
      }
    }
  }
  //- end of global loops

  //- normailzation for FFTW_BACKWARD
  if (!m_is_forward) {
    scal(field_out, 1.0 / Lxyz);
  }
}


//====================================================================
void FFT_xyz_3dim::fft(Field& field_out, const Field& field_in, const Direction dir)
{
  // save state
  bool backup_fwbw = m_is_forward;

  // find direction and set
  if (dir == FORWARD) {
    m_is_forward = true;
  } else if (dir == BACKWARD) {
    m_is_forward = false;
  } else {
    vout.crucial(m_vl, "%s: unknown direction %d. failed.\n", class_name.c_str(), dir);
    exit(EXIT_FAILURE);
  }

  // delegate to another method
  fft(field_out, field_in);

  // restore state
  m_is_forward = backup_fwbw;
}


//==========================================================
//==================================================END=====
#endif
