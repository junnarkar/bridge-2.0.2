/*!
        @file    action_F_Standard_lex_alt.cpp
        @brief
        @author  Matsufuru, Hideo (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$
        @version $LastChangedRevision: 2492 $
*/

#include "Action/Fermion/action_F_Standard_lex_alt.h"

#include "lib/Tools/timer.h"

#include "Field/index_lex_alt.h"
#include "Field/afield-inc.h"

template<typename AFIELD>
const std::string Action_F_Standard_lex_alt<AFIELD>::class_name
  = "Action_F_Standard_lex_alt";

//====================================================================
template<typename AFIELD>
void Action_F_Standard_lex_alt<AFIELD>::init()
{
  m_vl = CommonParameters::Vlevel();
}


//====================================================================
template<typename AFIELD>
void Action_F_Standard_lex_alt<AFIELD>::tidyup()
{
  // do nothing.
}


//====================================================================
template<typename AFIELD>
void Action_F_Standard_lex_alt<AFIELD>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);
}


//====================================================================
template<typename AFIELD>
void Action_F_Standard_lex_alt<AFIELD>::set_config(Field *U)
{
  m_U = U;

  m_fopr->set_config(U);
  m_fopr_force->set_config(U);
}


//====================================================================
template<typename AFIELD>
double Action_F_Standard_lex_alt<AFIELD>::langevin(RandomNumbers *rand)
{
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  int NinF     = m_fopr->field_nin();
  int NvolF    = m_fopr->field_nvol();
  int NexF     = m_fopr->field_nex();
  int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  assert(NvolF == Nvol);
  m_psf.reset(NinF, NvolF, NexF);

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  Field xi(NinF, NvolF, NexF);
  rand->gauss_lex_global(xi);

  AFIELD xiA(NinF, NvolF, NexF);

  Index_lex_alt<real_t> index_alt;

#pragma omp parallel
  {
    if (m_fopr->needs_convert()) {
      m_fopr->convert(xiA, xi);
    } else {
      convert_spinor(index_alt, xiA, xi);
    }
  }

  m_fopr->set_config(m_U);
  m_fopr->set_mode("Ddag");

  m_fopr->mult(m_psf, xiA);

  double xi2   = xi.norm();
  double H_psf = xi2 * xi2;

  vout.general(m_vl, "    H_Fstandard  = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);

  return H_psf;
}


//====================================================================
template<typename AFIELD>
double Action_F_Standard_lex_alt<AFIELD>::calcH()
{
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  int NinF     = m_fopr->field_nin();
  int NvolF    = m_fopr->field_nvol();
  int NexF     = m_fopr->field_nex();
  int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  AFIELD v1(NinF, NvolF, NexF);
  int    Nconv;
  double diff;

  m_fprop_H->set_config(m_U);
  m_fprop_H->set_mode("DdagD");
  m_fprop_H->invert(v1, m_psf, Nconv, diff);

  vout.general(m_vl, "    Fprop_H: %6d %18.15e\n", Nconv, diff);

  const double H_psf = dot(v1, m_psf);

  vout.general(m_vl, "    H_Fstandard  = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);

  return H_psf;
}


//====================================================================
template<typename AFIELD>
void Action_F_Standard_lex_alt<AFIELD>::force(Field& force)
{
  const int Nin  = m_U->nin();
  const int Nvol = m_U->nvol();
  const int Nex  = m_U->nex();
  const int Nc   = CommonParameters::Nc();
  const int Ndim = CommonParameters::Ndim();

  assert(force.nin() == Nin);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Nex);

  AFIELD force1(Nin, Nvol, Nex);

  const int NinF  = m_fopr->field_nin();
  const int NvolF = m_fopr->field_nvol();
  const int NexF  = m_fopr->field_nex();

  Timer  timer;
  double elapsed_time;
  timer.start();

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  AFIELD eta(NinF, NvolF, NexF);
  int    Nconv;
  double diff;

  m_fprop_MD->set_config(m_U);
  m_fprop_MD->set_mode("DdagD");
  m_fprop_MD->invert(eta, m_psf, Nconv, diff);

  vout.general(m_vl, "    Fprop_MD: %6d %18.15e\n", Nconv, diff);

  m_fopr->set_config(m_U);
  m_fopr_force->set_config(m_U);

  m_fopr_force->force_core(force1, eta);

  Index_lex_alt<real_t> index_lex;
#pragma omp parallel
  {
    reverse_gauge(index_lex, force, force1);
  }

  double Fave, Fmax, Fdev;
  force.stat(Fave, Fmax, Fdev);
  vout.general(m_vl, "    Fave = %12.6f  Fmax = %12.6f  Fdev = %12.6f\n",
               Fave, Fmax, Fdev);

  timer.stop();
  elapsed_time = timer.elapsed_sec();
  vout.general(m_vl, "    elapsed time: %12.6f [sec]\n", elapsed_time);
}


//============================================================END=====
