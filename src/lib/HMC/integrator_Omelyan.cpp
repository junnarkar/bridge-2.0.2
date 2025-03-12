/*!
        @file    integrator_Omelyan.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#include "integrator_Omelyan.h"

const std::string Integrator_Omelyan::class_name = "Integrator_Omelyan";

//====================================================================
void Integrator_Omelyan::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  int    level;
  int    Nstep;
  double lambda_Omelyan;

  int err = 0;
  err += params.fetch_int("level", level);
  err += params.fetch_int("number_of_steps", Nstep);
  err += params.fetch_double("lambda_Omelyan", lambda_Omelyan);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(level, Nstep, lambda_Omelyan);
}


//====================================================================
void Integrator_Omelyan::get_parameters(Parameters& params) const
{
  params.set_int("level", m_level);
  params.set_int("number_of_steps", m_Nstep);
  params.set_double("lambda_Omelyan", m_lambda);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Integrator_Omelyan::set_parameters(const int level, const int Nstep, const double lambda_Omelyan)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Level: %4d\n", level);
  vout.general(m_vl, "  Nstep    = %4d\n", Nstep);
  vout.general(m_vl, "  lambda_Omelyan = %10.6f\n", lambda_Omelyan);

  //- range check
  // NB. level,Estep,Nstep,lambda_Omelyan == 0 is allowed.

  //- store values
  m_Nstep  = Nstep;
  m_level  = level;
  m_lambda = lambda_Omelyan;
}


//====================================================================
void Integrator_Omelyan::set_parameter_level(const int level)
{
  m_level = level;
}


//====================================================================
void Integrator_Omelyan::set_parameter_Nstep(const int Nstep)
{
  m_Nstep = Nstep;
}


//====================================================================
void Integrator_Omelyan::set_parameter_Nsteps(const std::vector<int>& Nsteps)
{
  if (Nsteps.size() > 0) {
    set_parameter_Nstep(Nsteps[0]);

    // transfer to lower levels
    if (m_update_U && (Nsteps.size() > 1)) {
      std::vector<int> next_steps(Nsteps.begin() + 1, Nsteps.end());
      m_update_U->set_parameter_Nsteps(next_steps);
    }
  }
}


//====================================================================
void Integrator_Omelyan::set_parameter_lambda(const double lambda_omelyan)
{
  m_lambda = lambda_omelyan;
}


//====================================================================
void Integrator_Omelyan::evolve(const double step_size, Field_G& iP, Field_G& U)
{
  const int Nin  = U.nin();
  const int Nvol = U.nvol();
  const int Nex  = U.nex();
  const int Nc   = CommonParameters::Nc();

  vout.detailed(m_vl, "%s: level %d: Nstep = %d, step_size = %8.6f\n", class_name.c_str(), m_level, m_Nstep, step_size);

  vout.general(m_vl, "Integration level-%d start.\n", m_level);

  if (m_Nstep > 0) {
    double estep = step_size / m_Nstep;

    // initial lambda step
    m_update_p->evolve(estep * m_lambda, iP, U);

    // molecular dynamics steps
    for (int istep = 1; istep <= m_Nstep; ++istep) {
      m_update_U->evolve(estep * 0.5, iP, U);

      m_update_p->evolve(estep * (1.0 - 2.0 * m_lambda), iP, U);

      m_update_U->evolve(estep * 0.5, iP, U);

      if (istep < m_Nstep) {
        m_update_p->evolve(estep * m_lambda * 2.0, iP, U);
      }
    }

    // last lambda step
    m_update_p->evolve(estep * m_lambda, iP, U);
  } else {
    vout.general(m_vl, "Nstep is zero. skip.\n");
  }

  vout.general(m_vl, "Integration level-%d finished.\n", m_level);
}


//====================================================================
//============================================================END=====
