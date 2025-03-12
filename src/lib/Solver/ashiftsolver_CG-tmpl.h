/*!
        @file    ashiftsolver_CG-tmpl.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#include "ashiftsolver_CG.h"

template<typename FIELD, typename FOPR>
const std::string AShiftsolver_CG<FIELD, FOPR>::class_name
  = "AShiftsolver_CG";

//====================================================================
template<typename FIELD, typename FOPR>
void AShiftsolver_CG<FIELD, FOPR>::set_parameters(
  const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  int    Niter;
  double Stop_cond;

  int err = 0;
  err += params.fetch_int("maximum_number_of_iteration", Niter);
  err += params.fetch_double("convergence_criterion_squared", Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(Niter, Stop_cond);
}


//====================================================================
template<typename FIELD, typename FOPR>
void AShiftsolver_CG<FIELD, FOPR>::get_parameters(Parameters& params) const
{
  params.set_int("maximum_number_of_iteration", m_Niter);
  params.set_double("convergence_criterion_squared", m_Stop_cond);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename FIELD, typename FOPR>
void AShiftsolver_CG<FIELD, FOPR>::set_parameters(
  const int Niter,
  const double Stop_cond)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Niter     = %d\n", Niter);
  vout.general(m_vl, "  Stop_cond = %8.2e\n", Stop_cond);

  //- range check
  int err = 0;
  err += ParameterCheck::non_negative(Niter);
  err += ParameterCheck::square_non_zero(Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_Niter     = Niter;
  m_Stop_cond = Stop_cond;
}


//====================================================================
template<typename FIELD, typename FOPR>
void AShiftsolver_CG<FIELD, FOPR>::solve(
  std::vector<FIELD>& xq,
  const std::vector<double>& sigma,
  const FIELD& b,
  int& Nconv,
  double& diff)
{
  int Nshift = sigma.size();

  vout.paranoiac(m_vl, "  Shift CG solver start.\n");
  vout.paranoiac(m_vl, "    number of shift = %d\n", Nshift);
  vout.paranoiac(m_vl, "    values of shift:\n");
  for (int i = 0; i < Nshift; ++i) {
    vout.paranoiac(m_vl, "    %d  %12.8f\n", i, sigma[i]);
  }

  m_snorm = 1.0 / b.norm2();

  int Nconv2 = -1;

  reset_field(b, sigma, Nshift);

  copy(m_s, b);
  copy(m_r, b);

  double rr = 0.0;

  solve_init(rr);

  vout.detailed(m_vl, "    iter: %8d  %22.15e\n", 0, rr * m_snorm);

  bool is_converged = false;

  for (int iter = 0; iter < m_Niter; iter++) {
    solve_step(rr);

    Nconv2 += 1;

    vout.detailed(m_vl, "    iter: %8d  %22.15e  %4d\n",
                  (iter + 1), rr * m_snorm, m_Nshift2);

    if (rr * m_snorm < m_Stop_cond) {
      is_converged = true;
      break;
    }
  }

  if (!is_converged) {
    vout.crucial(m_vl, "Error at %s: not converged.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }


  std::vector<double> diffs(Nshift);
  for (int i = 0; i < Nshift; ++i) {
    diffs[i] = 0.0;
  }

  for (int i = 0; i < Nshift; ++i) {
    m_fopr->mult(m_s, m_x[i]);
    axpy(m_s, sigma[i], m_x[i]);
    axpy(m_s, -1.0, b);

    double diff1 = sqrt(m_s.norm2() * m_snorm);

    vout.paranoiac(m_vl, "    %4d  %22.15e\n", i, diff1);

    // if (diff1 > diff2) diff2 = diff1;
    diffs[i] = diff1;
  }

#pragma omp barrier
#pragma omp master
  {
    double diff2 = -1.0;

    for (int i = 0; i < Nshift; ++i) {
      if (diffs[i] > diff2) diff2 = diffs[i];
    }

    diff = diff2;

    Nconv = Nconv2;
  }
#pragma omp barrier

  for (int i = 0; i < Nshift; ++i) {
    copy(xq[i], m_x[i]);
  }

  vout.paranoiac(m_vl, "   diff(max) = %22.15e  \n", diff);
}


//====================================================================
template<typename FIELD, typename FOPR>
void AShiftsolver_CG<FIELD, FOPR>::solve_init(double& rr)
{
  int Nshift = m_p.size();

  vout.paranoiac(m_vl, "number of shift = %d\n", Nshift);

  for (int i = 0; i < Nshift; ++i) {
    copy(m_p[i], m_s);
    scal(m_x[i], 0.0);
  }

  copy(m_r, m_s);
  rr = m_r.norm2();

#pragma omp barrier
#pragma omp master
  {
    m_alpha_p = 0.0;
    m_beta_p  = 1.0;
  }
#pragma omp barrier
}


//====================================================================
template<typename FIELD, typename FOPR>
void AShiftsolver_CG<FIELD, FOPR>::solve_step(double& rr)
{
  m_fopr->mult(m_s, m_p[0]);
  axpy(m_s, m_sigma0, m_p[0]);

  double rr_p = rr;
  double pa_p = dot(m_s, m_p[0]);
  double beta = -rr_p / pa_p;

  axpy(m_x[0], -beta, m_p[0]);
  axpy(m_r, beta, m_s);
  rr = m_r.norm2();

  double alpha = rr / rr_p;

  aypx(alpha, m_p[0], m_r);

#pragma omp barrier
#pragma omp master
  {
    m_pp[0] = rr;
  }
#pragma omp barrier

  double alpha_h = 1.0 + m_alpha_p * beta / m_beta_p;

  for (int ish = 1; ish < m_Nshift2; ++ish) {
    double zeta = (alpha_h - m_csh2[ish] * beta) / m_zeta1[ish]
                  + (1.0 - alpha_h) / m_zeta2[ish];
    zeta = 1.0 / zeta;
    double zr      = zeta / m_zeta1[ish];
    double beta_s  = beta * zr;
    double alpha_s = alpha * zr * zr;

    axpy(m_x[ish], -beta_s, m_p[ish]);
    scal(m_p[ish], alpha_s);
    axpy(m_p[ish], zeta, m_r);

    double ppr = m_p[ish].norm2();

#pragma omp barrier
#pragma omp master
    {
      m_pp[ish] = ppr * m_snorm;

      m_zeta2[ish] = m_zeta1[ish];
      m_zeta1[ish] = zeta;
    }
#pragma omp barrier
  }

  int ish1 = m_Nshift2;

  for (int ish = m_Nshift2 - 1; ish >= 0; --ish) {
    vout.paranoiac(m_vl, "%4d %16.8e\n", ish, m_pp[ish]);
    if (m_pp[ish] > m_Stop_cond) {
      ish1 = ish + 1;
      break;
    }
  }

#pragma omp barrier
#pragma omp master
  {
    m_Nshift2 = ish1;

    m_alpha_p = alpha;
    m_beta_p  = beta;
  }
#pragma omp barrier
}


//====================================================================
template<typename FIELD, typename FOPR>
void AShiftsolver_CG<FIELD, FOPR>::reset_field(const FIELD& b,
                                               const std::vector<double>& sigma,
                                               const int Nshift)
{
#pragma omp barrier
#pragma omp master
  {
    int Nin  = b.nin();
    int Nvol = b.nvol();
    int Nex  = b.nex();

    m_p.resize(Nshift);
    m_x.resize(Nshift);
    m_zeta1.resize(Nshift);
    m_zeta2.resize(Nshift);
    m_csh2.resize(Nshift);
    m_pp.resize(Nshift);

    for (int i = 0; i < Nshift; ++i) {
      m_p[i].reset(Nin, Nvol, Nex);
      m_x[i].reset(Nin, Nvol, Nex);
      m_zeta1[i] = 1.0;
      m_zeta2[i] = 1.0;
      m_csh2[i]  = sigma[i] - sigma[0];

      m_pp[i] = 0.0;
    }

    m_s.reset(Nin, Nvol, Nex);
    m_r.reset(Nin, Nvol, Nex);

    m_sigma0 = sigma[0];

    m_Nshift2 = Nshift;
  }
#pragma omp barrier
}


//====================================================================
template<typename FIELD, typename FOPR>
double AShiftsolver_CG<FIELD, FOPR>::flop_count()
{
  vout.general(m_vl, "Warning at %s: flop_count() not yet implemented.\n",
               class_name.c_str());
  return 0.0;
}


//============================================================END=====
