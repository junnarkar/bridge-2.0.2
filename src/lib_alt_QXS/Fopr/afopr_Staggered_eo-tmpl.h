/*!
      @file    afopr_Staggered_eo-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2492 $
*/

#include "lib/ResourceManager/threadManager.h"

template<typename AFIELD>
const std::string AFopr_Staggered_eo<AFIELD>::class_name
  = "AFopr_Staggered_eo<AFIELD>";
//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  // switches
  int req_comm = 1;  // set 1 if communication forced any time
  //int req_comm = 0;  // set 0 if communication only when necessary

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  } else {
    m_vl = CommonParameters::Vlevel();
  }

  vout.general(m_vl, "%s: construction\n", class_name.c_str());

  m_Nc = CommonParameters::Nc();
  if (m_Nc != 3) {
    vout.crucial("%s: only applicable to Nc = 3\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_Nvc  = m_Nc * 2;
  m_Ndf  = 2 * m_Nc * m_Nc;
  m_Nx   = CommonParameters::Nx();
  m_Ny   = CommonParameters::Ny();
  m_Nz   = CommonParameters::Nz();
  m_Nt   = CommonParameters::Nt();
  m_Nst  = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();

  m_Nx2  = m_Nx / 2;
  m_Nst2 = m_Nst / 2;

  m_Nx2v  = m_Nx2 / VLENX;
  m_Nyv   = m_Ny / VLENY;
  m_Nst2v = m_Nst2 / VLEN;

  // condition check
  if (m_Nx % 2 != 0) {
    vout.crucial(m_vl, "%s: Nx must be even.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
  if (m_Nx % (2 * VLENX) != 0) {
    vout.crucial(m_vl, "%s: Nx must be mulriple of 2*VLENX.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
  if (m_Ny % (VLENY) != 0) {
    vout.crucial(m_vl, "%s: Ny must be multiple of VLENY.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  vout.general(m_vl, "  VLENX = %2d  Nx2v  = %d\n", VLENX, m_Nx2v);
  vout.general(m_vl, "  VLENY = %2d  Nyv   = %d\n", VLENY, m_Nyv);
  vout.general(m_vl, "  VLEN  = %2d  Nst2v = %d\n", VLEN, m_Nst2v);

  m_Leo.resize(m_Ny * m_Nz * m_Nt);

  int ipe3 = Communicator::ipe(3);
  int ipe2 = Communicator::ipe(2);
  int ipe1 = Communicator::ipe(1);
  for (int t = 0; t < m_Nt; ++t) {
    for (int z = 0; z < m_Nz; ++z) {
      for (int y = 0; y < m_Ny; ++y) {
        int t2 = ipe3 * m_Nt + t;
        int z2 = ipe2 * m_Nz + z;
        int y2 = ipe1 * m_Ny + y;
        m_Leo[y + m_Ny * (z + m_Nz * t)] = (y2 + z2 + t2) % 2;
      }
    }
  }

  m_Nsize[0] = m_Nx2v;
  m_Nsize[1] = m_Nyv;
  m_Nsize[2] = m_Nz;
  m_Nsize[3] = m_Nt;

  do_comm_any = 0;
  for (int mu = 0; mu < m_Ndim; ++mu) {
    do_comm[mu] = 1;
    if ((req_comm == 0) && (Communicator::npe(mu) == 1)) do_comm[mu] = 0;
    do_comm_any += do_comm[mu];
    vout.general(m_vl, "  do_comm[%d] = %d\n", mu, do_comm[mu]);
  }

  m_Nbdsize.resize(m_Ndim);
  m_Nbdsize[0] = m_Nvc * ((m_Ny * m_Nz * m_Nt + 1) / 2);
  m_Nbdsize[1] = m_Nvc * m_Nx2 * m_Nz * m_Nt;
  m_Nbdsize[2] = m_Nvc * m_Nx2 * m_Ny * m_Nt;
  m_Nbdsize[3] = m_Nvc * m_Nx2 * m_Ny * m_Nz;

  setup_channels();

  // gauge configuration.
  m_Ueo.reset(m_Ndf, m_Nst, m_Ndim);
  m_Ulex.reset(m_Ndf, m_Nst, m_Ndim);

  // working vectors.
  m_w1.reset(m_Nvc, m_Nst2, 1);
  m_w2.reset(m_Nvc, m_Nst2, 1);
  m_v1.reset(m_Nvc, m_Nst2, 1);
  m_v2.reset(m_Nvc, m_Nst2, 1);

  //  m_z1.reset(m_Nvc, m_Nst, 1);  // used in convert/revese.

  // staggered phase and parity
  m_stg_phase.reset(1, m_Nst, m_Ndim, Element_type::REAL);
  m_parity.reset(1, m_Nst, 1, Element_type::REAL);

  set_staggered_phase();

  vout.increase_indent();

  set_parameters(params);

  m_shift = new ShiftAField_eo<AFIELD>(m_Nvc);

  vout.decrease_indent();

  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::tidyup()
{
  delete m_shift;
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::setup_channels()
{
  chsend_up.resize(m_Ndim);
  chrecv_up.resize(m_Ndim);
  chsend_dn.resize(m_Ndim);
  chrecv_dn.resize(m_Ndim);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    int Nvsize = m_Nbdsize[mu] * sizeof(real_t);

    chsend_dn[mu].send_init(Nvsize, mu, -1);
    chsend_up[mu].send_init(Nvsize, mu, 1);
#ifdef USE_MPI
    chrecv_up[mu].recv_init(Nvsize, mu, 1);
    chrecv_dn[mu].recv_init(Nvsize, mu, -1);
#else
    void *buf_up = (void *)chsend_dn[mu].ptr();
    chrecv_up[mu].recv_init(Nvsize, mu, 1, buf_up);
    void *buf_dn = (void *)chsend_up[mu].ptr();
    chrecv_dn[mu].recv_init(Nvsize, mu, -1, buf_dn);
#endif

    if (do_comm[mu] == 1) {
      chset_send.append(chsend_up[mu]);
      chset_send.append(chsend_dn[mu]);
      chset_recv.append(chrecv_up[mu]);
      chset_recv.append(chrecv_dn[mu]);
    }
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

  double           mq;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("quark_mass", mq);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "%s: fetch error, input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(real_t(mq), bc);
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::set_parameters(const real_t mq,
                                                const std::vector<int> bc)
{
  assert(bc.size() == m_Ndim);

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    m_mq = mq;
    m_boundary.resize(m_Ndim);
    for (int mu = 0; mu < m_Ndim; ++mu) {
      m_boundary[mu] = bc[mu];
    }
  }

  vout.general(m_vl, "%s: set parameters\n", class_name.c_str());
  vout.general(m_vl, "  mq   = %8.4f\n", m_mq);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_double("quark_mass", double(m_mq));
  params.set_int_vector("boundary_condition", m_boundary);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::set_staggered_phase()
{
  int ipex = Communicator::ipe(0);
  int ipey = Communicator::ipe(1);
  int ipez = Communicator::ipe(2);
  int ipet = Communicator::ipe(3);

  real_t *stgph = m_stg_phase.ptr(0);
  real_t *prty  = m_parity.ptr(0);

  AIndex_lex<real_t, AFIELD::IMPL> index;

  for (int t = 0; t < m_Nt; ++t) {
    int t2 = t + ipet * m_Nt;
    for (int z = 0; z < m_Nz; ++z) {
      int z2 = z + ipez * m_Nz;
      for (int y = 0; y < m_Ny; ++y) {
        int y2 = y + ipey * m_Ny;
        for (int x = 0; x < m_Nx; ++x) {
          int x2 = x + ipex * m_Nx;
          int is = index.site(x, y, z, t);

          stgph[index.idx(0, 1, is, 0)] = 1.0;
          stgph[index.idx(0, 1, is, 1)] = 1.0;
          stgph[index.idx(0, 1, is, 2)] = 1.0;
          stgph[index.idx(0, 1, is, 3)] = 1.0;

          prty[index.idx(0, 1, is, 0)] = 1.0;

          if ((x2 % 2) == 1) {
            stgph[index.idx(0, 1, is, 1)] = -1.0;
          }
          if (((x2 + y2) % 2) == 1) {
            stgph[index.idx(0, 1, is, 2)] = -1.0;
          }
          if (((x2 + y2 + z2) % 2) == 1) {
            stgph[index.idx(0, 1, is, 3)] = -1.0;
          }
          if (((x2 + y2 + z2 + t2) % 2) == 1) {
            prty[index.idx(0, 1, is, 0)] = -1.0;
          }
        }
      }
    }
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::set_config(Field *u)
{
  int nth = ThreadManager::get_num_threads();

  vout.detailed(m_vl, "%s: set_config is called: num_threads = %d\n",
                class_name.c_str(), nth);

  if (nth > 1) {
    set_config_impl(u);
  } else {
    set_config_omp(u);
  }

  vout.detailed(m_vl, "%s: set_config finished\n", class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::set_config_omp(Field *u)
{
  vout.detailed(m_vl, "  set_config_omp is called.\n");

#pragma omp parallel
  {
    set_config_impl(u);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::set_config_impl(Field *U)
{
  AIndex_lex<real_t, AFIELD::IMPL> index_lex;
  AIndex_eo<real_t, AFIELD::IMPL>  index_eo;

  real_t *Ueo = m_Ueo.ptr(0);

  int Nsize_lex[4];
  Nsize_lex[0] = m_Nsize[0] * 2;
  Nsize_lex[1] = m_Nsize[1];
  Nsize_lex[2] = m_Nsize[2];
  Nsize_lex[3] = m_Nsize[3];

  convert_gauge(index_lex, m_Ulex, *U);

  QXS_Gauge::set_boundary(m_Ulex, m_boundary);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    // staggered phase is multiplied to gauge field
    real_t *stgph = m_stg_phase.ptr(index_lex.idx(0, 1, 0, mu));
    real_t *up    = m_Ulex.ptr(index_lex.idx_G(0, 0, mu));
    BridgeQXS::mult_staggered_phase(up, stgph, Nsize_lex, m_Ndf);
  }
#pragma omp barrier

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nst);

  for (int ex = 0; ex < m_Ndim; ++ex) {
    for (int site = is; site < ns; ++site) {
      for (int in = 0; in < m_Ndf; ++in) {
        int iv1 = index_lex.idx(in, m_Ndf, site, ex);
        int iv2 = index_eo.idx(in, m_Ndf, site, ex);
        Ueo[iv2] = m_Ulex.cmp(iv1);
      }
    }
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::normalize_fopr(AFIELD& v)
{
#pragma omp barrier
  scal(v, real_t(1.0) / (m_mq * m_mq));
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::normalize_fprop(AFIELD& v)
{
#pragma omp barrier
  scal(v, m_mq * m_mq);
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::set_mode(std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    D(v, w);
  } else if (m_mode == "Ddag") {
    Ddag(v, w);
  } else if (m_mode == "DdagD") {
    DdagD(v, w);
    //  } else if (m_mode == "H") {
    //    H(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undeifined.\n", class_name.c_str());
    abort();
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult_dag(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    Ddag(v, w);
  } else if (m_mode == "Ddag") {
    D(v, w);
  } else if (m_mode == "DdagD") {
    DdagD(v, w);
    //} else if (m_mode == "H") {
    //   H(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undeifined.\n", class_name.c_str());
    abort();
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult(AFIELD& v, const AFIELD& w,
                                      const std::string mode)
{
#pragma omp barrier

  if (mode == "Dee") {
    copy(v, w);
    scal(v, m_mq);
  } else if (mode == "Doo") {
    copy(v, w);
    scal(v, m_mq);
  } else if (mode == "Dee_inv") {
    copy(v, w);
    scal(v, real_t(1.0) / m_mq);
  } else if (mode == "Doo_inv") {
    copy(v, w);
    scal(v, real_t(1.0) / m_mq);
  } else if (mode == "Deo") {
    Meo(v, w, 0);
  } else if (mode == "Doe") {
    Meo(v, w, 1);
  } else {
    vout.crucial(m_vl, "%s: illegal mode is given to mult with mode\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult_dag(AFIELD& v, const AFIELD& w,
                                          const std::string mode)
{
#pragma omp barrier

  if (mode == "Dee") {
    copy(v, w);
    scal(v, m_mq);
  } else if (mode == "Doo") {
    copy(v, w);
    scal(v, m_mq);
  } else if (mode == "Dee_inv") {
    copy(v, w);
    scal(v, real_t(1.0) / m_mq);
  } else if (mode == "Doo_inv") {
    copy(v, w);
    scal(v, real_t(1.0) / m_mq);
  } else if (mode == "Deo") {
    Meo(v, w, 0);
    scal(v, real_t(-1.0));
  } else if (mode == "Doe") {
    Meo(v, w, 1);
    scal(v, real_t(-1.0));
  } else {
    vout.crucial(m_vl, "%s: illegal mode is given to mult with mode\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::H(AFIELD& v, const AFIELD& w)
{
  D(v, w);
  mult_gm5(v);
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::DdagD(AFIELD& v, const AFIELD& w)
{
  D(m_v2, w);
  Ddag(v, m_v2);
  // H(m_v2, w);
  // H(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::D(AFIELD& v, const AFIELD& w)
{
  //  Meo(m_v1, w, w, 1, 0);
  //  Meo(v, m_v1, w, 0, 1);

  Meo(m_v1, w, w, 1, 0);
  Meo(v, m_v1, w, 0, 0);

  // aypx(real_t(-1.0/(m_mq * m_mq)), v, w);

  scal(v, real_t(-1.0));
  axpy(v, m_mq * m_mq, w);
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::Ddag(AFIELD& v, const AFIELD& w)
{
  //  Meo(m_v1, w, w, 1, 0);
  //  Meo(v, m_v1, w, 0, 1);

  Meo(m_v1, w, w, 1, 0);
  Meo(v, m_v1, w, 0, 0);

  scal(v, real_t(-1.0));
  axpy(v, m_mq * m_mq, w);
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::Meo(AFIELD& v, const AFIELD& w,
                                     int ieo)
{
  mult_Meo_qxs(v, w, w, ieo, 0);
  //mult_Meo_alt(v, w, w, ieo, 0);
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::Meo(AFIELD& v, const AFIELD& w,
                                     const AFIELD& x,
                                     int ieo, int iflag)
{
  mult_Meo_qxs(v, w, x, ieo, iflag);
  //mult_Meo_alt(v, w, x, ieo, iflag);
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult_Meo_qxs(AFIELD& v,
                                              const AFIELD& w,
                                              const AFIELD& x,
                                              int ieo, int iflag)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
  real_t *xp = const_cast<AFIELD *>(&x)->ptr(0);
  real_t *up = m_Ueo.ptr(0);

  int ith = ThreadManager::get_thread_id();

#pragma omp barrier

  if (do_comm_any > 0) {
    if (ith == 0) chset_recv.start();

    real_t *buf1xp = (real_t *)chsend_dn[0].ptr();
    real_t *buf1xm = (real_t *)chsend_up[0].ptr();

    real_t *buf1yp = (real_t *)chsend_dn[1].ptr();
    real_t *buf1ym = (real_t *)chsend_up[1].ptr();

    real_t *buf1zp = (real_t *)chsend_dn[2].ptr();
    real_t *buf1zm = (real_t *)chsend_up[2].ptr();

    real_t *buf1tp = (real_t *)chsend_dn[3].ptr();
    real_t *buf1tm = (real_t *)chsend_up[3].ptr();

    BridgeQXS::mult_staggered_eo_1(buf1xp, buf1xm, buf1yp, buf1ym,
                                   buf1zp, buf1zm, buf1tp, buf1tm,
                                   up, wp,
                                   m_Nsize, do_comm, &m_Leo[0], ieo);

#pragma omp barrier

    if (ith == 0) chset_send.start();
  }

  int jd = 1;
  BridgeQXS::mult_staggered_eo_bulk(vp, up, wp, xp, m_mq, jd,
                                    m_Nsize, do_comm,
                                    &m_Leo[0], ieo, iflag);

  if (do_comm_any > 0) {
    if (ith == 0) chset_recv.wait();

#pragma omp barrier

    real_t *buf2xp = (real_t *)chrecv_up[0].ptr();
    real_t *buf2xm = (real_t *)chrecv_dn[0].ptr();

    real_t *buf2yp = (real_t *)chrecv_up[1].ptr();
    real_t *buf2ym = (real_t *)chrecv_dn[1].ptr();

    real_t *buf2zp = (real_t *)chrecv_up[2].ptr();
    real_t *buf2zm = (real_t *)chrecv_dn[2].ptr();

    real_t *buf2tp = (real_t *)chrecv_up[3].ptr();
    real_t *buf2tm = (real_t *)chrecv_dn[3].ptr();

    BridgeQXS::mult_staggered_eo_2(vp, up, wp,
                                   buf2xp, buf2xm, buf2yp, buf2ym,
                                   buf2zp, buf2zm, buf2tp, buf2tm,
                                   m_mq, m_Nsize, do_comm,
                                   &m_Leo[0], ieo, iflag);

    if (ith == 0) chset_send.wait();
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult_Meo_alt(AFIELD& v,
                                              const AFIELD& w,
                                              const AFIELD& x,
                                              int ieo, int iflag)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
  real_t *xp = const_cast<AFIELD *>(&x)->ptr(0);

#pragma omp barrier

  clear(vp);
#pragma omp barrier

  for (int mu = 0; mu < m_Ndim; ++mu) {
    mult_up(mu, v, w, ieo);
    mult_dn(mu, v, w, ieo);
  }

#pragma omp barrier

  // note that mass normalization is adopted.
  if (iflag != 0) {
    real_t fac = -1.0 / (m_mq * m_mq);
    axpby(fac, vp, real_t(1.0), xp);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult_gm5(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier

  copy(v, w);
#pragma omp barrier

  /*
  real_t* ph = m_parity.ptr(0);
  real_t* vp = v.ptr(0);

  BridgeQXS::mult_staggered_phase(vp, ph, m_Nsize, m_Nvc);

#pragma omp barrier
  */
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult_gm5(AFIELD& v)
{
  real_t *ph = m_parity.ptr(0);
  real_t *vp = v.ptr(0);

#pragma omp barrier
  BridgeQXS::mult_staggered_phase(vp, ph, m_Nsize, m_Nvc);
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::clear(real_t *vp)
{
  BridgeQXS::mult_staggered_clear(vp, m_Nsize, m_Nvc);
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::axpby(real_t b, real_t *vp,
                                       real_t a, real_t *wp)
{
  BridgeQXS::mult_staggered_axpby(b, vp, a, wp, m_Nsize, m_Nvc);
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult_up(int idir,
                                         AFIELD& v, const AFIELD& w)
{
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult_dn(int idir,
                                         AFIELD& v, const AFIELD& w)
{
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult_up(int idir,
                                         AFIELD& v, const AFIELD& w,
                                         int ieo)
{
  m_shift->backward(m_w1, w, idir, ieo);

  real_t *wp = m_w1.ptr(0);
  real_t *up = m_Ueo.ptr(NDF * m_Nst2 * (ieo + 2 * idir));
  real_t *vp = m_w2.ptr(0);

  BridgeQXS::mult_staggered_mult_Gn(vp, up, wp, m_Nsize);
#pragma omp barrier

  //  copy(m_w2, m_w1);

  axpy(v, real_t(0.5), m_w2);
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult_dn(int idir,
                                         AFIELD& v, const AFIELD& w,
                                         int ieo)
{
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
  real_t *up = m_Ueo.ptr(NDF * m_Nst2 * (1 - ieo + 2 * idir));
  real_t *vp = m_w1.ptr(0);

  BridgeQXS::mult_staggered_mult_Gd(vp, up, wp, m_Nsize);
#pragma omp barrier

  //  copy(m_w1, w);
  //#pragma omp barrier

  m_shift->forward(m_w2, m_w1, idir, ieo);

  axpy(v, real_t(-0.5), m_w2);
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
double AFopr_Staggered_eo<AFIELD>::flop_count()
{
  return flop_count(m_mode);
}


//====================================================================
template<typename AFIELD>
double AFopr_Staggered_eo<AFIELD>::flop_count(const std::string mode)
{
  // The following counting explicitly depends on the implementation.
  // It will be recalculated when the code is modified.
  // The following is based on rev.1976.   [21 Jul 2019 H.Matsufuru]

  int    Lvol = CommonParameters::Lvol();
  double flop_site, flop;

  flop_site = static_cast<double>(m_Nvc * (2 + 8 * 2 * m_Nvc));
  //  #comp   aypx  dir FMA  #comp

  flop = flop_site * static_cast<double>(Lvol);

  if ((mode == "DdagD") || (mode == "DDdag")) flop *= 2.0;

  return flop;
}


//============================================================END=====
