/*
        @file    asolver_BiCGStab-tmpl.h
        @brief
        @author  Hideo Matsufuru (matufuru)
        $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2#$
        @version $LastChangedRevision: 2492 $
*/

template<typename AFIELD>
const std::string ASolver_BiCGStab<AFIELD>::class_name
  = "ASolver_BiCGStab";
//#define DEBUG_BICGSTAB
//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab<AFIELD>::init(void)
{
  ThreadManager::assert_single_thread(class_name);

  vout.detailed(m_vl, "%s: being setup.\n", class_name.c_str());

  int nin  = m_fopr->field_nin();
  int nvol = m_fopr->field_nvol();
  int nex  = m_fopr->field_nex();
  vout.detailed(m_vl, "  Field size: Nin = %d  Nvol = %d  Nex = %d\n",
                nin, nvol, nex);

  m_r.reset(nin, nvol, nex);
  m_x.reset(nin, nvol, nex);
  m_p.reset(nin, nvol, nex);
  m_s.reset(nin, nvol, nex);
  m_v.reset(nin, nvol, nex);
  m_t.reset(nin, nvol, nex);
  m_rh.reset(nin, nvol, nex);

  m_ecrit           = 1.e-32;
  m_nconv           = -1;
  m_Omega_tolerance = 0.7;

  m_initial_mode = InitialGuess::RHS;

  vout.detailed(m_vl, "%s: setup finished.\n", class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab<AFIELD>::tidyup(void)
{
  // ThreadManager::assert_single_thread(class_name);
  // nothing is to be deleted.
}


//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab<AFIELD>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  int    Niter, Nrestart;
  double Stop_cond;
  double Omega_tolerance;

  int err = 0;
  err += params.fetch_int("maximum_number_of_iteration", Niter);
  err += params.fetch_int("maximum_number_of_restart", Nrestart);
  err += params.fetch_double("convergence_criterion_squared", Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
  err = params.fetch_double("Omega_tolerance", Omega_tolerance);
  if (err) {
    Omega_tolerance = 0.7;
  }

  // - default initial mode is RHS
  // - negative Stop_cond is for fixed nubmer of iterations
  int Niter2 = Niter * Nrestart;
  set_parameters(Niter2, Stop_cond, InitialGuess::RHS);
  set_parameters_BiCGStab_series(Omega_tolerance);
}


//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab<AFIELD>::set_parameters(const int Niter,
                                              const real_t Stop_cond)
{
  set_parameters(Niter, Stop_cond, InitialGuess::RHS);
}


//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab<AFIELD>::set_parameters(const int Niter,
                                              const real_t Stop_cond,
                                              const InitialGuess init_guess_mode)

{
  ThreadManager::assert_single_thread(class_name);

  m_Niter     = Niter;
  m_Stop_cond = Stop_cond;
  std::string prec = "double";
  if (sizeof(real_t) == 4) prec = "float";
  m_initial_mode = init_guess_mode;

  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Precision: %s\n", prec.c_str());
  vout.general(m_vl, "  Niter     = %d\n", m_Niter);
  vout.general(m_vl, "  Stop_cond = %16.8e\n", m_Stop_cond);
  vout.general(m_vl, "  init_guess_mode: %d\n", m_initial_mode);
}


//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab<AFIELD>::set_parameters_BiCGStab_series(const real_t Omega_tolerance)
{
  ThreadManager::assert_single_thread(class_name);

  //- print input parameters
  vout.general(m_vl, "  Omega_tolerance = %8.2e\n", Omega_tolerance);

  //- range check
  // NB. Omega_tolerance == 0.0 is allowed.

  //- store values
  m_Omega_tolerance = Omega_tolerance;
}


//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab<AFIELD>::solve(AFIELD& xq, const AFIELD& b,
                                     int& Nconv, real_t& diff)
{
  vout.paranoiac(m_vl, "%s: solver start:\n", class_name.c_str());
#pragma omp master
  {
    m_init_mult = 0;
  }
  copy(m_s, b);
  real_t sr     = norm2(m_s);
  real_t snorm  = real_t(1.0) / sr;
  real_t snorm2 = sqrt(snorm);
  vout.detailed(m_vl, "  sr =    %22.15e\n", sr);
  vout.detailed(m_vl, "  snorm = %22.15e\n", snorm);

  real_t alpha_prev, rho_prev, omega_prev;
  scal(m_s, snorm2);

  if (m_initial_mode == InitialGuess::RHS) {
    {
      vout.detailed(m_vl, "%s: initial_mode = RHS\n", class_name.c_str());
#ifdef DEBUG_BICGSTAB
      double s2 = m_s.norm2();

      vout.detailed(m_vl, "%s: m_x: %d %d %d\n", class_name.c_str(), m_x.nin(), m_x.nvol(), m_x.nex());
      vout.detailed(m_vl, "%s: m_s: %d %d %d:  nrom2=%23.16e\n", class_name.c_str(), m_s.nin(), m_s.nvol(), m_s.nex(), s2);
#endif
    }
    copy(m_x, m_s);
  } else if (m_initial_mode == InitialGuess::GIVEN) {
    vout.detailed(m_vl, "%s: initial_mode = GIVEN\n", class_name.c_str());
    copy(m_x, xq);
    scal(m_x, snorm2);
  } else if (m_initial_mode == InitialGuess::ZERO) {
    vout.detailed(m_vl, "%s: initial_mode = ZERO\n", class_name.c_str());
  } else {
    vout.crucial("%s: unkown init guess mode\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  int nconv = -1;

  real_t rr;
  solve_init(b, rr, alpha_prev, rho_prev, omega_prev, m_initial_mode);
  vout.detailed(m_vl, "  init: %22.15e\n", rr);

  for (int iter = 0; iter < m_Niter; ++iter) {
    int iflg = 0;
    solve_step(rr, iflg, alpha_prev, rho_prev, omega_prev);
#ifdef DEBUG_BICGSTAB
    {
      m_fopr->mult(xq, m_x);
      copy(m_t, b);
      m_t.scal(snorm2);
      axpy(xq, -1.0, m_t);
      double r2 = xq.norm2();
      vout.detailed(m_vl, "   iter=%d, alg_rr=%22.15e,  actual r2=%22.15e\n", iter, rr, r2);
    }
#endif
    if (iflg != 0) {
      copy(m_s, b);
      scal(m_s, snorm2);
      solve_init(b, rr, alpha_prev, rho_prev, omega_prev, InitialGuess::GIVEN);
    }
    vout.detailed(m_vl, "%6d  %22.15e\n", iter, rr);

    // if(rr*snorm < m_Stop_cond)
    if (rr < m_Stop_cond) {
      nconv = 2 * (iter + 1);  // counting fermion mult
      break;
    }
  }

  if (nconv == -1) {
    if ((m_Stop_cond < 0) || (Nconv < 0)) {
      vout.detailed(m_vl, "  truncating with iter= %d\n", 2 * m_Niter);
      nconv = 2 * m_Niter;
    } else {
      vout.crucial(m_vl, "Error at %s: not converged\n",
                   class_name.c_str());
      vout.crucial(m_vl, "  iter(final): %8d  %22.15e\n",
                   m_Niter, rr);
      exit(EXIT_FAILURE);
    }
  } else {
    vout.detailed(m_vl, "converged:\n");
  }
  vout.detailed(m_vl, "  nconv = %d\n", nconv);

  copy(xq, m_x);
  real_t sr2 = sqrt(sr);
  scal(xq, sr2);

  m_fopr->mult(m_s, xq);

  axpy(m_s, real_t(-1.0), b);
  real_t diff2 = norm2(m_s);

#pragma omp master
  {
    m_nconv = nconv;
    Nconv   = nconv + m_init_mult; // include solve_init();
    diff    = sqrt(diff2 / sr);
  }
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab<AFIELD>::solve_init(const AFIELD& b,
                                          real_t& rr,
                                          real_t& alpha_prev, real_t& rho_prev, real_t& omega_prev,
                                          const InitialGuess init_mode)
{
  if (init_mode == InitialGuess::ZERO) {
    m_x.set(0.0);
    m_v.set(0.0);
    copy(m_r, m_s);
  } else {
    m_fopr->mult(m_v, m_x);
    copy(m_r, m_s);
    axpy(m_r, real_t(-1.0), m_v);
#pragma omp master
    {
      m_init_mult++;
    }
  }
#ifdef DEBUG_BICGSTAB
  {
    double x2 = m_x.norm2();
    double v2 = m_v.norm2();
    vout.general(m_vl, "%s: x2=%23.15e\n", class_name.c_str(), x2);
    vout.general(m_vl, "%s: v2=%23.15e\n", class_name.c_str(), v2);
  }
#endif


  copy(m_rh, m_r);
  rr = norm2(m_r);
#pragma omp barrier
  m_p.set(real_t(0.0));
  m_v.set(real_t(0.0));

  rho_prev   = real_t(1.0);
  alpha_prev = real_t(1.0);
  omega_prev = real_t(1.0);
}


//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab<AFIELD>::solve_step(real_t& rr, int& iflg,
                                          real_t& alpha_prev, real_t& rho_prev, real_t& omega_prev)
{
  real_t rho = dot(m_rh, m_r);
  if (fabs(rho) < m_ecrit) {
    iflg = 1;
    vout.detailed(m_vl, "too small value of rho = %e\n", rho);
    return;
  }
  real_t bet = rho * alpha_prev / (rho_prev * omega_prev);
#ifdef DEBUG_BICGSTAB
  vout.general(m_vl, "%s: rho        = %23.15e\n", class_name.c_str(), rho);
  vout.general(m_vl, "%s: alpha_prev = %23.15e\n", class_name.c_str(), alpha_prev);
  vout.general(m_vl, "%s: rho_prev   = %23.15e\n", class_name.c_str(), rho_prev);
  vout.general(m_vl, "%s: omega_prev = %23.15e\n", class_name.c_str(), omega_prev);
  vout.general(m_vl, "%s: bet        = %23.15e\n", class_name.c_str(), bet);
#endif

  // p = r + bet * (p - omega_prev * v);
  axpy(m_p, -omega_prev, m_v);
  aypx(bet, m_p, m_r);

  m_fopr->mult(m_v, m_p);
  real_t aden = dot(m_rh, m_v);

#ifdef DEBUG_BICGSTAB
  {
    double p2 = m_p.norm2();
    double v2 = m_v.norm2();
    vout.general(m_vl, "%s: p2=%23.15e\n", class_name.c_str(), p2);
    vout.general(m_vl, "%s: v2=%23.15e\n", class_name.c_str(), v2);
    vout.general(m_vl, "%s: aden=%23.15e\n", class_name.c_str(), aden);
  }
#endif

  if (fabs(aden) < m_ecrit) {
    iflg = 1;
    vout.detailed(m_vl, "too small denominator aden = %e\n", aden);
    return;
  }
  real_t alpha = rho / aden;
  //    if(fabs(alpha) >100.0) {
  //      iflg = 1;
  //      vout.detailed(m_vl, "too large alpha = %e\n", alpha);
  //      return;
  //    }
  copy(m_s, m_r);
  axpy(m_s, -alpha, m_v);

  m_fopr->mult(m_t, m_s);

  real_t omega_n = dot(m_t, m_s);
  real_t omega_d = m_t.norm2();
  real_t s_norm2 = m_s.norm2();
#ifdef DEBUG_BICGSTAB
  vout.general(m_vl, "%s: alpha   = %23.15e\n", class_name.c_str(), alpha);
  vout.general(m_vl, "%s: omega_n = %23.15e\n", class_name.c_str(), omega_n);
  vout.general(m_vl, "%s: omega_d = %23.15e\n", class_name.c_str(), omega_d);
  vout.general(m_vl, "%s: s_norm2 = %23.15e\n", class_name.c_str(), s_norm2);
#endif
  if (fabs(omega_d) < m_ecrit) {
    iflg = 1;
    vout.detailed(m_vl, "too small denominator omega_d = %e\n", omega_d);
    return;
  }
  real_t omega   = omega_n / omega_d;
  real_t abs_rho = abs(omega_n) / sqrt(omega_d * s_norm2);
#ifdef DEBUG_BICGSTAB
  vout.general(m_vl, "%s: omega  =%23.15e\n", class_name.c_str(), omega);
  vout.general(m_vl, "%s: abs_rho=%23.15e\n", class_name.c_str(), abs_rho);
#endif
  if (abs_rho < m_Omega_tolerance) {
    vout.detailed(m_vl, "Omega_prescription: abs_rho = %23.16e, omega= %23.16e\n", abs_rho, omega);
    omega *= m_Omega_tolerance / abs_rho;
  }
  axpy(m_x, omega, m_s);
  axpy(m_x, alpha, m_p);

  copy(m_r, m_s);
  axpy(m_r, -omega, m_t);

  rr = norm2(m_r);
#ifdef DEBUG_BICGSTAB
  vout.general(m_vl, "%s: rr = %23.15e\n", class_name.c_str(), rr);
#endif
  rho_prev   = rho;
  alpha_prev = alpha;
  omega_prev = omega;
}


//====================================================================
template<typename AFIELD>
double ASolver_BiCGStab<AFIELD>::flop_count()
{
  int Nin  = m_fopr->field_nin();
  int Nvol = m_fopr->field_nvol();
  int Nex  = m_fopr->field_nex();
  int NPE  = CommonParameters::NPE();

  int ninit = m_init_mult;

  double flop_field  = static_cast<double>(Nin * Nvol * Nex) * NPE;
  double flop_vector = (6 + ninit * 4 + m_nconv * 11) * flop_field;
  double flop_fopr   = (1 + ninit + m_nconv) * m_fopr->flop_count();

  if (m_initial_mode == InitialGuess::ZERO) {
    flop_vector -= m_init_mult * 2 * flop_field;
  }
  double flop = flop_vector + flop_fopr;

  return flop;
}


//====================================================================
//============================================================END=====
