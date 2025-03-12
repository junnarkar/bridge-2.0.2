/*!
        @file    force_F_Rational.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#include "force_F_Rational.h"
#include "lib/Solver/shiftsolver_CG.h"

const std::string Force_F_Rational::class_name = "Force_F_Rational";

//====================================================================
void Force_F_Rational::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  int    Np, n_exp, d_exp;
  double x_min, x_max;
  int    Niter;
  double Stop_cond;

  int err = 0;
  err += params.fetch_int("number_of_poles", Np);
  err += params.fetch_int("exponent_numerator", n_exp);
  err += params.fetch_int("exponent_denominator", d_exp);
  err += params.fetch_double("lower_bound", x_min);
  err += params.fetch_double("upper_bound", x_max);
  err += params.fetch_int("maximum_number_of_iteration", Niter);
  err += params.fetch_double("convergence_criterion_squared", Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(Np, n_exp, d_exp, x_min, x_max, Niter, Stop_cond);
}


//====================================================================
void Force_F_Rational::get_parameters(Parameters& params) const
{
  params.set_int("number_of_poles", m_Np);
  params.set_int("exponent_numerator", m_n_exp);
  params.set_int("exponent_denominator", m_d_exp);
  params.set_double("lower_bound", m_x_min);
  params.set_double("upper_bound", m_x_max);
  params.set_int("maximum_number_of_iteration", m_Niter);
  params.set_double("convergence_criterion_squared", m_Stop_cond);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Force_F_Rational::set_parameters(const int Np, const int n_exp, const int d_exp,
                                      const double x_min, const double x_max,
                                      const int Niter, const double Stop_cond)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Np        = %d\n", Np);
  vout.general(m_vl, "  n_exp     = %d\n", n_exp);
  vout.general(m_vl, "  d_exp     = %d\n", d_exp);
  vout.general(m_vl, "  x_min     = %12.8f\n", x_min);
  vout.general(m_vl, "  x_max     = %12.8f\n", x_max);
  vout.general(m_vl, "  Niter     = %d\n", Niter);
  vout.general(m_vl, "  Stop_cond = %8.2e\n", Stop_cond);

  //- range check
  int err = 0;
  err += ParameterCheck::non_zero(Np);
  err += ParameterCheck::non_zero(n_exp);
  err += ParameterCheck::non_zero(d_exp);
  // NB. x_min,x_max=0 is allowed.
  err += ParameterCheck::non_zero(Niter);
  err += ParameterCheck::square_non_zero(Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_Np        = Np;
  m_n_exp     = n_exp;
  m_d_exp     = d_exp;
  m_x_min     = x_min;
  m_x_max     = x_max;
  m_Niter     = Niter;
  m_Stop_cond = Stop_cond;

  //- post-process
  m_cl.resize(m_Np);
  m_bl.resize(m_Np);

  //- Rational approximation
  const double x_min2 = m_x_min * m_x_min;
  const double x_max2 = m_x_max * m_x_max;

  Math_Rational rational;
  rational.set_parameters(m_Np, m_n_exp, m_d_exp, x_min2, x_max2);
  rational.get_parameters(m_a0, m_bl, m_cl);

  vout.general(m_vl, " a0 = %18.14f\n", m_a0);
  for (int i = 0; i < m_Np; i++) {
    vout.general(m_vl, " bl[%d] = %18.14f  cl[%d] = %18.14f\n",
                 i, m_bl[i], i, m_cl[i]);
  }
}


//====================================================================
void Force_F_Rational::force_udiv(Field& force_, const Field& eta_)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  Field_G force(Nvol, Ndim);
  Field_F eta(eta_);

  force_udiv_impl(force, eta);

  copy(force_, force); // force_ = force;
}


//====================================================================
void Force_F_Rational::force_udiv_impl(Field_G& force, const Field_F& eta)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  const int NinF  = eta.nin();
  const int NvolF = eta.nvol();
  const int NexF  = eta.nex();

  //- Shiftsolver
  const int Nshift = m_Np;

  std::vector<Field> psi(Nshift);

  for (int i = 0; i < Nshift; ++i) {
    psi[i].reset(NinF, NvolF, NexF);
  }

  vout.general(m_vl, "    Shift solver in force calculation\n");
  vout.general(m_vl, "      Number of shift values = %d\n", m_cl.size());
  m_fopr->set_mode("DdagD");

  Shiftsolver_CG solver(m_fopr, m_Niter, m_Stop_cond);
  int            Nconv;
  double         diff;
  solver.solve(psi, m_cl, eta, Nconv, diff);
  vout.general(m_vl, "      diff(max) = %22.15e  \n", diff);

  force.set(0.0);

  for (int i = 0; i < Nshift; ++i) {
    Field_G force1(Nvol, Ndim);
    m_force->force_udiv(force1, psi[i]);
    scal(force1, m_bl[i]);    // force1 *= m_bl[i];
    axpy(force, 1.0, force1); // force  += force1;
  }
}


//====================================================================
void Force_F_Rational::force_core1(Field&, const Field&, const Field&)
{
  vout.crucial(m_vl, "Error at %s: not implemented.\n", __func__);
  exit(EXIT_FAILURE);
}


//====================================================================
void Force_F_Rational::force_udiv1(Field&, const Field&, const Field&)
{
  vout.crucial(m_vl, "Error at %s: not implemented.\n", __func__);
  exit(EXIT_FAILURE);
}


//====================================================================
//===========================================================END======
