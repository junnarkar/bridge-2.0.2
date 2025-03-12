/*!
      @file    action_F_Rational_alt.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-02-28 16:09:41 #$
      @version $LastChangedRevision: 2492 $
*/

//#include "Action/Fermion/action_F_Rational_alt.h"
//#include "Field/index_lex_alt.h"
//#include "Field/afield-inc.h"

template<typename AFIELD>
const std::string Action_F_Rational_alt<AFIELD>::class_name
  = "Action_F_Rational_alt";

//====================================================================
template<typename AFIELD>
void Action_F_Rational_alt<AFIELD>::init()
{
  m_vl = CommonParameters::Vlevel();
}


//====================================================================
template<typename AFIELD>
void Action_F_Rational_alt<AFIELD>::tidyup()
{
  // do nothing.
}


//====================================================================
template<typename AFIELD>
void Action_F_Rational_alt<AFIELD>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);
}


//====================================================================
template<typename AFIELD>
void Action_F_Rational_alt<AFIELD>::set_config(Field *U)
{
  m_U = U;
  m_fopr_langev->set_config(U);
  m_fopr_H->set_config(U);
  m_fopr_force_MD->set_config(U);
}


//====================================================================
template<typename AFIELD>
double Action_F_Rational_alt<AFIELD>::langevin(RandomNumbers *rand)
{
  const int NinF     = m_fopr_langev->field_nin();
  const int NvolF    = m_fopr_langev->field_nvol();
  const int NexF     = m_fopr_langev->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  m_psf.reset(NinF, NvolF, NexF);

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  Field xi(NinF, NvolF, NexF);
  rand->gauss_lex_global(xi);

  AFIELD xiA(NinF, NvolF, NexF);

  Index_lex_alt<real_t, AFIELD::IMPL> index_alt;

#pragma omp parallel
  {
    if (m_fopr_langev->needs_convert()) {
      m_fopr_langev->convert(xiA, xi);
    } else {
      convert(index_alt, xiA, xi);
    }
  }

  m_fopr_langev->set_config(m_U);
  m_fopr_langev->mult(m_psf, xiA);

  double xi2   = xi.norm();
  double H_psf = xi2 * xi2;

  vout.general(m_vl, "    H_Frational  = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);

  return H_psf;
}


//====================================================================
template<typename AFIELD>
double Action_F_Rational_alt<AFIELD>::calcH()
{
  const int NinF     = m_fopr_H->field_nin();
  const int NvolF    = m_fopr_H->field_nvol();
  const int NexF     = m_fopr_H->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  AFIELD v1(NinF, NvolF, NexF);
  m_fopr_H->set_config(m_U);
  m_fopr_H->mult(v1, m_psf);

  double H_psf = dot(v1, m_psf);

  vout.general(m_vl, "    H_Frational  = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);

  return H_psf;
}


//====================================================================
template<typename AFIELD>
void Action_F_Rational_alt<AFIELD>::force(Field& force)
{
  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  Timer  timer;
  double elapsed_time;
  timer.start();

  const int Nin  = m_U->nin();
  const int Nvol = m_U->nvol();
  const int Nex  = m_U->nex();

  assert(force.nin() == Nin);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Nex);

  AFIELD force1(Nin, Nvol, Nex);

  m_fopr_force_MD->set_config(m_U);
  m_fopr_force_MD->force_core(force1, m_psf);

  Index_lex_alt<real_t, AFIELD::IMPL> index_lex;
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


//====================================================================
//============================================================END=====
