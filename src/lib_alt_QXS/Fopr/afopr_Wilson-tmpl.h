/*!
      @file    afopr_Wilson-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2492 $
*/

#include "lib/ResourceManager/threadManager.h"

template<typename AFIELD>
const std::string AFopr_Wilson<AFIELD>::class_name = "AFopr_Wilson";

//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  // switch of coomunication
  int req_comm = 1;  // set 1 if communication forced any time
  //int req_comm = 0;  // set 0 if communication called in necessary

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  } else {
    m_vl = CommonParameters::Vlevel();
  }

  vout.general(m_vl, "%s: construction\n", class_name.c_str());

  m_repr = "Dirac";  // now only the Dirac repr is available.

  std::string repr;
  if (!params.fetch_string("gamma_matrix_type", repr)) {
    if (repr != "Dirac") {
      vout.crucial("  Error at %s: unsupported gamma-matrix type: %s\n",
                   class_name.c_str(), repr.c_str());
      exit(EXIT_FAILURE);
    }
  }

  m_Nc = CommonParameters::Nc();
  if (m_Nc != 3) {
    vout.crucial("%s: only applicable to Nc = 3\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_Nd  = CommonParameters::Nd();
  m_Nvc = 2 * m_Nc;
  m_Ndf = 2 * m_Nc * m_Nc;

  m_Nx   = CommonParameters::Nx();
  m_Ny   = CommonParameters::Ny();
  m_Nz   = CommonParameters::Nz();
  m_Nt   = CommonParameters::Nt();
  m_Nst  = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();

  m_Nxv  = m_Nx / VLENX;
  m_Nyv  = m_Ny / VLENY;
  m_Nstv = m_Nst / VLEN;

  if (VLENX * m_Nxv != m_Nx) {
    vout.crucial(m_vl, "%s: Nx must be multiple of VLENX.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
  if (VLENY * m_Nyv != m_Ny) {
    vout.crucial(m_vl, "%s: Ny must be multiple of VLENY.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  vout.general(m_vl, "  VLENX = %2d  Nxv  = %d\n", VLENX, m_Nxv);
  vout.general(m_vl, "  VLENY = %2d  Nyv  = %d\n", VLENY, m_Nyv);
  vout.general(m_vl, "  VLEN  = %2d  Nstv = %d\n", VLEN, m_Nstv);

  m_Nsize[0] = m_Nxv;
  m_Nsize[1] = m_Nyv;
  m_Nsize[2] = m_Nz;
  m_Nsize[3] = m_Nt;

  do_comm_any = 0;
  for (int mu = 0; mu < m_Ndim; ++mu) {
    do_comm[mu] = 1;
    if ((req_comm == 0) && (Communicator::npe(mu) == 1)) do_comm[mu] = 0;
    do_comm_any += do_comm[mu];
    vout.general("  do_comm[%d] = %d\n", mu, do_comm[mu]);
  }

  m_bdsize.resize(m_Ndim);
  int Nd2 = m_Nd / 2;
  m_bdsize[0] = m_Nvc * Nd2 * m_Ny * m_Nz * m_Nt;
  m_bdsize[1] = m_Nvc * Nd2 * m_Nx * m_Nz * m_Nt;
  m_bdsize[2] = m_Nvc * Nd2 * m_Nx * m_Ny * m_Nt;
  m_bdsize[3] = m_Nvc * Nd2 * m_Nx * m_Ny * m_Nz;

  setup_channels();

  // gauge configuration.
  m_U.reset(NDF, m_Nst, m_Ndim);

  // working vectors.
  int NinF = 2 * m_Nc * m_Nd;
  m_v2.reset(NinF, m_Nst, 1);

  vout.increase_indent();

  set_parameters(params);

  vout.decrease_indent();

  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::setup_channels()
{
  chsend_up.resize(m_Ndim);
  chrecv_up.resize(m_Ndim);
  chsend_dn.resize(m_Ndim);
  chrecv_dn.resize(m_Ndim);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    int    tag_up, tag_dn;
    size_t Nvsize = m_bdsize[mu] * sizeof(real_t);

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
void AFopr_Wilson<AFIELD>::tidyup()
{
  ThreadManager::assert_single_thread(class_name);

  // nothing to do at present.
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::set_parameters(const Parameters& params)
{
  double           kappa;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_int_vector("boundary_condition", bc);
  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(real_t(kappa), bc);
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::set_parameters(const real_t CKs,
                                          const std::vector<int> bc)
{
  assert(bc.size() == m_Ndim);

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    m_CKs = CKs;
    m_boundary.resize(m_Ndim);
    for (int mu = 0; mu < m_Ndim; ++mu) {
      m_boundary[mu] = bc[mu];
    }
  }

  vout.general(m_vl, "%s: set parameters\n", class_name.c_str());
  vout.general(m_vl, "  gamma-matrix type = %s\n", m_repr.c_str());
  vout.general(m_vl, "  kappa = %8.4f\n", m_CKs);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_double("hopping_parameter", double(m_CKs));
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_string("gamma_matrix_type", m_repr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::set_config(Field *u)
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
void AFopr_Wilson<AFIELD>::set_config_omp(Field *u)
{
  vout.detailed(m_vl, "  set_config_omp is called.\n");

#pragma omp parallel
  {
    set_config_impl(u);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::set_config_impl(Field *u)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) m_conf = u;

  AIndex_lex<real_t, AFIELD::IMPL> index;

  convert_gauge(index, m_U, *u);

  QXS_Gauge::set_boundary(m_U, m_boundary);
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::convert(AFIELD& v, const Field& w)
{
  AIndex_lex<real_t, AFIELD::IMPL> index_lex;
  convert_spinor(index_lex, m_v2, w);

  mult_gm4(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::reverse(Field& v, const AFIELD& w)
{
  mult_gm4(m_v2, w);

  AIndex_lex<real_t, AFIELD::IMPL> index_lex;
  reverse_spinor(index_lex, v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::mult_up(int mu, AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  if (mu == 0) {
    mult_xp(vp, wp);
  } else if (mu == 1) {
    mult_yp(vp, wp);
  } else if (mu == 2) {
    mult_zp(vp, wp);
  } else if (mu == 3) {
    mult_tp(vp, wp);
  } else {
    vout.crucial(m_vl, "%s: mult_up for %d direction is undefined.",
                 class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::mult_dn(int mu, AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  if (mu == 0) {
    mult_xm(vp, wp);
  } else if (mu == 1) {
    mult_ym(vp, wp);
  } else if (mu == 2) {
    mult_zm(vp, wp);
  } else if (mu == 3) {
    mult_tm(vp, wp);
  } else {
    vout.crucial(m_vl, "%s: mult_dn for %d direction is undefined.",
                 class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::set_mode(std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
std::string AFopr_Wilson<AFIELD>::get_mode() const
{
  return m_mode;
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::mult(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    return D(v, w);
  } else if (m_mode == "DdagD") {
    return DdagD(v, w);
  } else if (m_mode == "Ddag") {
    return Ddag(v, w);
  } else if (m_mode == "H") {
    return H(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::mult_dag(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    return Ddag(v, w);
  } else if (m_mode == "DdagD") {
    return DdagD(v, w);
  } else if (m_mode == "Ddag") {
    return D(v, w);
  } else if (m_mode == "H") {
    return H(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::D(AFIELD& v, const AFIELD& w)
{
  mult_D(v, w);
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::DdagD(AFIELD& v, const AFIELD& w)
{
  D(m_v2, w);
  mult_gm5(v, m_v2);
  D(m_v2, v);
  mult_gm5(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::Ddag(AFIELD& v, const AFIELD& w)
{
  mult_gm5(v, w);
  D(m_v2, v);
  mult_gm5(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::mult_gm5(AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

#pragma omp barrier

  //  mult_gm5(vp, wp);

  BridgeQXS::mult_wilson_gm5_dirac(vp, wp, m_Nsize);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::mult_gm4(AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  Vsimd_t wt[2];

  for (int site = is; site < ns; ++site) {
    for (int ic = 0; ic < NC; ++ic) {
      for (int id = 0; id < ND2; ++id) {
        int idx1 = 2 * (id + ND * ic) + NVCD * site;
        load_vec(wt, &wp[VLEN * idx1], 2);
        save_vec(&vp[VLEN * idx1], wt, 2);
      }

      for (int id = ND2; id < ND; ++id) {
        int idx1 = 2 * (id + ND * ic) + NVCD * site;
        load_vec(wt, &wp[VLEN * idx1], 2);
        scal_vec(wt, real_t(-1.0), 2);
        save_vec(&vp[VLEN * idx1], wt, 2);
      }
    }
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::mult_D(AFIELD& v, const AFIELD& w)
{
  real_t *v2 = v.ptr(0);
  real_t *v1 = const_cast<AFIELD *>(&w)->ptr(0);
  real_t *up = m_U.ptr(0);

  int ith = ThreadManager::get_thread_id();

#pragma omp barrier

  if (do_comm_any > 0) {
    if (ith == 0) chset_recv.start();

    real_t *buf1_xp = (real_t *)chsend_dn[0].ptr();
    real_t *buf1_xm = (real_t *)chsend_up[0].ptr();
    real_t *buf1_yp = (real_t *)chsend_dn[1].ptr();
    real_t *buf1_ym = (real_t *)chsend_up[1].ptr();
    real_t *buf1_zp = (real_t *)chsend_dn[2].ptr();
    real_t *buf1_zm = (real_t *)chsend_up[2].ptr();
    real_t *buf1_tp = (real_t *)chsend_dn[3].ptr();
    real_t *buf1_tm = (real_t *)chsend_up[3].ptr();

    BridgeQXS::mult_wilson_1_dirac(buf1_xp, buf1_xm, buf1_yp, buf1_ym,
                                   buf1_zp, buf1_zm, buf1_tp, buf1_tm,
                                   up, v1, &m_boundary[0], m_Nsize, do_comm);

#pragma omp barrier

    if (ith == 0) chset_send.start();
  }

  BridgeQXS::mult_wilson_bulk_dirac(v2, up, v1,
                                    m_CKs, &m_boundary[0], m_Nsize, do_comm);

  if (do_comm_any > 0) {
    if (ith == 0) chset_recv.wait();

#pragma omp barrier

    real_t *buf2_xp = (real_t *)chrecv_up[0].ptr();
    real_t *buf2_xm = (real_t *)chrecv_dn[0].ptr();
    real_t *buf2_yp = (real_t *)chrecv_up[1].ptr();
    real_t *buf2_ym = (real_t *)chrecv_dn[1].ptr();
    real_t *buf2_zp = (real_t *)chrecv_up[2].ptr();
    real_t *buf2_zm = (real_t *)chrecv_dn[2].ptr();
    real_t *buf2_tp = (real_t *)chrecv_up[3].ptr();
    real_t *buf2_tm = (real_t *)chrecv_dn[3].ptr();

    BridgeQXS::mult_wilson_2_dirac(v2, up, v1,
                                   buf2_xp, buf2_xm, buf2_yp, buf2_ym,
                                   buf2_zp, buf2_zm, buf2_tp, buf2_tm,
                                   m_CKs, &m_boundary[0], m_Nsize, do_comm);

    if (ith == 0) chset_send.wait();
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::mult_D_alt(AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  clear(vp);

  mult_xp(vp, wp);
  mult_xm(vp, wp);
  mult_yp(vp, wp);
  mult_ym(vp, wp);
  mult_zp(vp, wp);
  mult_zm(vp, wp);
  mult_tp(vp, wp);
  mult_tm(vp, wp);

  aypx(-m_CKs, vp, wp);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::H(AFIELD& v, const AFIELD& w)
{
  D(m_v2, w);
  mult_gm5(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::aypx(real_t a, real_t *v, real_t *w)
{
  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  Vsimd_t vt[NVCD], wt[NVCD];

  for (int site = is; site < ns; ++site) {
    load_vec(vt, &v[VLEN * NVCD * site], NVCD);
    load_vec(wt, &w[VLEN * NVCD * site], NVCD);
    aypx_vec(a, vt, wt, NVCD);
    save_vec(&v[VLEN * NVCD * site], vt, NVCD);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::clear(real_t *v)
{
  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  Vsimd_t vt[NVCD];
  clear_vec(vt, NVCD);

  for (int site = is; site < ns; ++site) {
    save_vec(&v[VLEN * NVCD * site], vt, NVCD);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::mult_xp(real_t *v2, real_t *v1)
{
  int idir = 0;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  Vsimd_t v2v[NVCD];

  real_t *buf1 = (real_t *)chsend_dn[0].ptr();
  real_t *buf2 = (real_t *)chrecv_up[0].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

#pragma omp barrier

  if (do_comm[0] > 0) {
    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nxv;
      int iyzt = site / m_Nxv;
      if (ix == 0) {
        int ibf = VLENY * NVC * ND2 * iyzt;
        mult_wilson_xp1(&buf1[ibf], &v1[VLEN * NVCD * site]);
      }
    }

#pragma omp barrier

#pragma omp master
    {
      chrecv_up[0].start();
      chsend_dn[0].start();
      chrecv_up[0].wait();
      chsend_dn[0].wait();
    }
#pragma omp barrier
  } // if(do_comm[0] == 1)

  for (int site = is; site < ns; ++site) {
    int ix   = site % m_Nxv;
    int iyzt = site / m_Nxv;

    Vsimd_t v2v[NVCD];
    clear_vec(v2v, NVCD);

    real_t zL[VLEN * NVCD];

    if ((ix < m_Nxv - 1) || (do_comm[0] == 0)) {
      int nei = ix + 1 + m_Nxv * iyzt;
      if (ix == m_Nxv - 1) nei = 0 + m_Nxv * iyzt;
      shift_vec2_xbw(zL, &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei], NVCD);
      mult_wilson_xpb(v2v, &u[VLEN * NDF * site], zL);
    } else {
      int ibf = VLENY * NVC * ND2 * iyzt;
      shift_vec0_xbw(zL, &v1[VLEN * NVCD * site], NVCD);
      mult_wilson_xpb(v2v, &u[VLEN * NDF * site], zL);
      mult_wilson_xp2(v2v, &u[VLEN * NDF * site], &buf2[ibf]);
    }

    add_vec(&v2[VLEN * NVCD * site], v2v, NVCD);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::mult_xm(real_t *v2, real_t *v1)
{
  int idir = 0;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  Vsimd_t v2v[NVCD];

  real_t *buf1 = (real_t *)chsend_up[0].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[0].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

#pragma omp barrier

  if (do_comm[0] > 0) {
    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nxv;
      int iyzt = site / m_Nxv;
      if (ix == m_Nxv - 1) {
        int ibf = VLENY * NVC * ND2 * iyzt;
        mult_wilson_xm1(&buf1[ibf], &u[VLEN * NDF * site],
                        &v1[VLEN * NVCD * site]);
      }
    }

#pragma omp barrier
#pragma omp master
    {
      chrecv_dn[0].start();
      chsend_up[0].start();
      chrecv_dn[0].wait();
      chsend_up[0].wait();
    }
#pragma omp barrier
  } // end of if(do_comm[0] > 0)

  for (int site = is; site < ns; ++site) {
    int ix   = site % m_Nxv;
    int iyzt = site / m_Nxv;

    real_t zL[VLEN * NVCD];
    real_t uL[VLEN * NDF];

    clear_vec(v2v, NVCD);

    if ((ix > 0) || (do_comm[0] == 0)) {
      int nei = ix - 1 + m_Nxv * iyzt;
      if (ix == 0) nei = m_Nxv - 1 + m_Nxv * iyzt;
      shift_vec2_xfw(zL, &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei], NVCD);
      shift_vec2_xfw(uL, &u[VLEN * NDF * site], &u[VLEN * NDF * nei], NDF);
      mult_wilson_xmb(v2v, uL, zL);
    } else {
      int ibf = VLENY * NVC * ND2 * iyzt;
      shift_vec0_xfw(zL, &v1[VLEN * NVCD * site], NVCD);
      shift_vec0_xfw(uL, &u[VLEN * NDF * site], NDF);
      mult_wilson_xmb(v2v, uL, zL);
      mult_wilson_xm2(v2v, &buf2[ibf]);
    }

    add_vec(&v2[VLEN * NVCD * site], v2v, NVCD);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::mult_yp(real_t *v2, real_t *v1)
{
  int idir = 1;
  int Nxy  = m_Nxv * m_Nyv;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_dn[1].ptr();
  real_t *buf2 = (real_t *)chrecv_up[1].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

#pragma omp barrier

  if (do_comm[1] > 0) {
    for (int site = is; site < ns; ++site) {
      int ix  = site % m_Nxv;
      int iy  = (site / m_Nxv) % m_Nyv;
      int izt = site / Nxy;
      if (iy == 0) {
        int ibf = VLENX * NVC * ND2 * (ix + m_Nxv * izt);
        mult_wilson_yp1(&buf1[ibf], &v1[VLEN * NVCD * site]);
      }
    }

#pragma omp barrier

#pragma omp master
    {
      chrecv_up[1].start();
      chsend_dn[1].start();
      chrecv_up[1].wait();
      chsend_dn[1].wait();
    }

#pragma omp barrier
  }  // end of if(do_comm[1] > 0)

  for (int site = is; site < ns; ++site) {
    int ix  = site % m_Nxv;
    int iy  = (site / m_Nxv) % m_Nyv;
    int izt = site / Nxy;

    Vsimd_t v2v[NVCD];
    clear_vec(v2v, NVCD);

    real_t zL[VLEN * NVCD];

    if ((iy < m_Nyv - 1) || (do_comm[1] == 0)) {
      int iy2 = (iy + 1) % m_Nyv;
      int nei = ix + m_Nxv * (iy2 + m_Nyv * izt);
      shift_vec2_ybw(zL, &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei], NVCD);
      mult_wilson_ypb(v2v, &u[VLEN * NDF * site], zL);
    } else {
      int ibf = VLENX * NVC * ND2 * (ix + m_Nxv * izt);
      shift_vec0_ybw(zL, &v1[VLEN * NVCD * site], NVCD);
      mult_wilson_ypb(v2v, &u[VLEN * NDF * site], zL);
      mult_wilson_yp2(v2v, &u[VLEN * NDF * site], &buf2[ibf]);
    }

    add_vec(&v2[VLEN * NVCD * site], v2v, NVCD);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::mult_ym(real_t *v2, real_t *v1)
{
  int idir = 1;
  int Nxy  = m_Nxv * m_Nyv;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_up[1].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[1].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

#pragma omp barrier

  if (do_comm[1] > 0) {
    for (int site = is; site < ns; ++site) {
      int ix  = site % m_Nxv;
      int iy  = (site / m_Nxv) % m_Nyv;
      int izt = site / Nxy;
      if (iy == m_Nyv - 1) {
        int ibf = VLENX * NVC * ND2 * (ix + m_Nxv * izt);
        mult_wilson_ym1(&buf1[ibf], &u[VLEN * NDF * site],
                        &v1[VLEN * NVCD * site]);
      }
    }

#pragma omp barrier

#pragma omp master
    {
      chrecv_dn[1].start();
      chsend_up[1].start();
      chrecv_dn[1].wait();
      chsend_up[1].wait();
    }

#pragma omp barrier
  }

  for (int site = is; site < ns; ++site) {
    int ix  = site % m_Nxv;
    int iy  = (site / m_Nxv) % m_Nyv;
    int izt = site / Nxy;

    Vsimd_t v2v[NVCD];
    clear_vec(v2v, NVCD);

    real_t zL[VLEN * NVCD];
    real_t uL[VLEN * NDF];

    if ((iy != 0) || (do_comm[idir] == 0)) {
      int iy2 = (iy - 1 + m_Nyv) % m_Nyv;
      int nei = ix + m_Nxv * (iy2 + m_Nyv * izt);
      shift_vec2_yfw(zL, &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei], NVCD);
      shift_vec2_yfw(uL, &u[VLEN * NDF * site], &u[VLEN * NDF * nei], NDF);
      mult_wilson_ymb(v2v, uL, zL);
    } else {
      int ibf = VLENX * NVC * ND2 * (ix + m_Nxv * izt);
      shift_vec0_yfw(zL, &v1[VLEN * NVCD * site], NVCD);
      shift_vec0_yfw(uL, &u[VLEN * NDF * site], NDF);
      mult_wilson_ymb(v2v, uL, zL);
      mult_wilson_ym2(v2v, &buf2[ibf]);
    }

    add_vec(&v2[VLEN * NVCD * site], v2v, NVCD);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::mult_zp(real_t *v2, real_t *v1)
{
  int idir = 2;
  int Nxy  = m_Nxv * m_Nyv;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_dn[2].ptr();
  real_t *buf2 = (real_t *)chrecv_up[2].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

#pragma omp barrier

  if (do_comm[2] > 0) {
    for (int site = is; site < ns; ++site) {
      int ixy = site % Nxy;
      int iz  = (site / Nxy) % m_Nz;
      int it  = site / (Nxy * m_Nz);
      if (iz == 0) {
        int ibf = VLEN * NVC * ND2 * (ixy + Nxy * it);
        mult_wilson_zp1(&buf1[ibf], &v1[VLEN * NVCD * site]);
      }
    }

#pragma omp barrier

#pragma omp master
    {
      chrecv_up[2].start();
      chsend_dn[2].start();
      chrecv_up[2].wait();
      chsend_dn[2].wait();
    }

#pragma omp barrier
  }

  for (int site = is; site < ns; ++site) {
    int ixy = site % Nxy;
    int iz  = (site / Nxy) % m_Nz;
    int it  = site / (Nxy * m_Nz);

    Vsimd_t v2v[NVCD];
    clear_vec(v2v, NVCD);

    if ((iz != m_Nz - 1) || (do_comm[2] == 0)) {
      int iz2 = (iz + 1) % m_Nz;
      int nei = ixy + Nxy * (iz2 + m_Nz * it);
      mult_wilson_zpb(v2v, &u[VLEN * NDF * site], &v1[VLEN * NVCD * nei]);
    } else {
      int ibf = VLEN * NVC * ND2 * (ixy + Nxy * it);
      mult_wilson_zp2(v2v, &u[VLEN * NDF * site], &buf2[ibf]);
    }

    add_vec(&v2[VLEN * NVCD * site], v2v, NVCD);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::mult_zm(real_t *v2, real_t *v1)
{
  int idir = 2;
  int Nxy  = m_Nxv * m_Nyv;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_up[2].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[2].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

#pragma omp barrier

  if (do_comm[2] > 0) {
    for (int site = is; site < ns; ++site) {
      int ixy = site % Nxy;
      int iz  = (site / Nxy) % m_Nz;
      int it  = site / (Nxy * m_Nz);
      if (iz == m_Nz - 1) {
        int ibf = VLEN * NVC * ND2 * (ixy + Nxy * it);
        mult_wilson_zm1(&buf1[ibf], &u[VLEN * NDF * site],
                        &v1[VLEN * NVCD * site]);
      }
    }

#pragma omp barrier

#pragma omp master
    {
      chrecv_dn[2].start();
      chsend_up[2].start();
      chrecv_dn[2].wait();
      chsend_up[2].wait();
    }

#pragma omp barrier
  }

  for (int site = is; site < ns; ++site) {
    int ixy = site % Nxy;
    int iz  = (site / Nxy) % m_Nz;
    int it  = site / (Nxy * m_Nz);

    Vsimd_t v2v[NVCD];
    clear_vec(v2v, NVCD);

    if ((iz > 0) || (do_comm[2] == 0)) {
      int iz2 = (iz - 1 + m_Nz) % m_Nz;
      int nei = ixy + Nxy * (iz2 + m_Nz * it);
      mult_wilson_zmb(v2v, &u[VLEN * NDF * nei], &v1[VLEN * NVCD * nei]);
    } else {
      int ibf = VLEN * NVC * ND2 * (ixy + Nxy * it);
      mult_wilson_zm2(v2v, &buf2[ibf]);
    }

    add_vec(&v2[VLEN * NVCD * site], v2v, NVCD);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::mult_tp(real_t *v2, real_t *v1)
{
  int idir = 3;
  int Nxyz = m_Nxv * m_Nyv * m_Nz;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_dn[3].ptr();
  real_t *buf2 = (real_t *)chrecv_up[3].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

#pragma omp barrier

  if (do_comm[3] > 0) {
    for (int site = is; site < ns; ++site) {
      int ixyz = site % Nxyz;
      int it   = site / Nxyz;
      if (it == 0) {
        mult_wilson_tp1_dirac(&buf1[VLEN * NVC * ND2 * ixyz],
                              &v1[VLEN * NVCD * site]);
      }
    }

#pragma omp barrier

#pragma omp master
    {
      chrecv_up[3].start();
      chsend_dn[3].start();
      chrecv_up[3].wait();
      chsend_dn[3].wait();
    }

#pragma omp barrier
  }

  for (int site = is; site < ns; ++site) {
    int ixyz = site % Nxyz;
    int it   = site / Nxyz;

    Vsimd_t v2v[NVCD];
    clear_vec(v2v, NVCD);

    if ((it < m_Nt - 1) || (do_comm[3] == 0)) {
      int it2 = (it + 1) % m_Nt;
      int nei = ixyz + Nxyz * it2;
      mult_wilson_tpb_dirac(v2v, &u[VLEN * NDF * site],
                            &v1[VLEN * NVCD * nei]);
    } else {
      mult_wilson_tp2_dirac(v2v, &u[VLEN * NDF * site],
                            &buf2[VLEN * NVC * ND2 * ixyz]);
    }

    add_vec(&v2[VLEN * NVCD * site], v2v, NVCD);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson<AFIELD>::mult_tm(real_t *v2, real_t *v1)
{
  int idir = 3;
  int Nxyz = m_Nxv * m_Nyv * m_Nz;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_up[3].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[3].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

#pragma omp barrier

  if (do_comm[3] > 0) {
    for (int site = is; site < ns; ++site) {
      int ixyz = site % Nxyz;
      int it   = site / Nxyz;
      if (it == m_Nt - 1) {
        mult_wilson_tm1_dirac(&buf1[VLEN * NVC * ND2 * ixyz],
                              &u[VLEN * NDF * site], &v1[VLEN * NVCD * site]);
      }
    }

#pragma omp barrier

#pragma omp master
    {
      chrecv_dn[3].start();
      chsend_up[3].start();
      chrecv_dn[3].wait();
      chsend_up[3].wait();
    }
#pragma omp barrier
  }

  for (int site = is; site < ns; ++site) {
    int ixyz = site % Nxyz;
    int it   = site / Nxyz;

    Vsimd_t v2v[NVCD];
    clear_vec(v2v, NVCD);

    if ((it > 0) || (do_comm[3] == 0)) {
      int it2 = (it - 1 + m_Nt) % m_Nt;
      int nei = ixyz + Nxyz * it2;
      mult_wilson_tmb_dirac(v2v, &u[VLEN * NDF * nei],
                            &v1[VLEN * NVCD * nei]);
    } else {
      mult_wilson_tm2_dirac(v2v, &buf2[VLEN * NVC * ND2 * ixyz]);
    }

    add_vec(&v2[VLEN * NVCD * site], v2v, NVCD);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
double AFopr_Wilson<AFIELD>::flop_count(const std::string mode)
{
  // The following counting explicitly depends on the implementation.
  // It will be recalculated when the code is modified.

  int    Lvol = CommonParameters::Lvol();
  double flop_site, flop;

  if (m_repr == "Dirac") {
    flop_site = static_cast<double>(
      m_Nc * m_Nd * (4 + 6 * (4 * m_Nc + 2) + 2 * (4 * m_Nc + 1)));
  } else if (m_repr == "Chiral") {
    flop_site = static_cast<double>(
      m_Nc * m_Nd * (4 + 8 * (4 * m_Nc + 2)));
  } else {
    vout.crucial(m_vl, "%s: input repr is undefined.\n");
    abort();
  }

  flop = flop_site * static_cast<double>(Lvol);
  if ((mode == "DdagD") || (mode == "DDdag")) flop *= 2.0;

  return flop;
}


//============================================================END=====
