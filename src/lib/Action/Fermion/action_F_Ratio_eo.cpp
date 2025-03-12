/*!
        @file    action_F_Ratio_eo.cpp

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#include "action_F_Ratio_eo.h"

const std::string Action_F_Ratio_eo::class_name = "Action_F_Ratio_eo";

//====================================================================
void Action_F_Ratio_eo::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }
}


//====================================================================
void Action_F_Ratio_eo::get_parameters(Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Action_F_Ratio_eo::set_parameters()
{
  vout.general(m_vl, "%s:\n", class_name.c_str());
}


//====================================================================
void Action_F_Ratio_eo::set_config(Field *U)
{
  m_U = U;

  //- NB. only solver part is even-odd preconditioned.
  m_fopr_prec->set_config(U);
  m_fopr_prec_force->set_config(U);
  m_fopr->set_config(U);
  m_fopr_force->set_config(U);
}


//====================================================================
double Action_F_Ratio_eo::langevin(RandomNumbers *rand)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  const int NinF     = m_fopr_prec->field_nin();
  const int NvolF    = m_fopr_prec->field_nvol();
  const int NexF     = m_fopr_prec->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  assert(NvolF == Nvol);
  m_psf.reset(NinF, NvolF, NexF);

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  Field xi(NinF, NvolF, NexF);
  rand->gauss_lex_global(xi);

  m_fopr_prec->set_config(m_U);
  m_fopr->set_config(m_U);

  Field v2(NinF, NvolF, NexF);
  m_fopr->set_mode("H");
  m_fopr->mult_dag(v2, xi);

  Field  v1(NinF, NvolF, NexF);
  int    Nconv;
  double diff;
  m_fprop_H_prec->set_config(m_U);
  m_fprop_H_prec->invert_DdagD(v1, v2, Nconv, diff);
  vout.general(m_vl, "    Nconv = %d  diff  = %.8e\n", Nconv, diff);

  m_fopr_prec->set_mode("H");
  m_fopr_prec->mult(m_psf, v1);

  const double H_psf = xi.norm2();

  vout.general(m_vl, "    H_Fratio     = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);

  return H_psf;
}


//====================================================================
double Action_F_Ratio_eo::calcH()
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  const int NinF     = m_fopr_prec->field_nin();
  const int NvolF    = m_fopr_prec->field_nvol();
  const int NexF     = m_fopr_prec->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  m_fopr_prec->set_config(m_U);
  m_fopr->set_config(m_U);

  Field v1(NinF, NvolF, NexF);
  m_fopr_prec->set_mode("H");
  m_fopr_prec->mult_dag(v1, m_psf);

  Field  v2(NinF, NvolF, NexF);
  int    Nconv;
  double diff;
  m_fprop_H->set_config(m_U);
  m_fprop_H->invert_DdagD(v2, v1, Nconv, diff);
  vout.general(m_vl, "    Nconv = %d  diff  = %.8e\n", Nconv, diff);

  const double H_psf = dot(v1, v2);

  vout.general(m_vl, "    H_Fratio     = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);

  return H_psf;
}


//====================================================================
void Action_F_Ratio_eo::force(Field& force)
{
  const int Nin  = m_U->nin();
  const int Nvol = m_U->nvol();
  const int Nex  = m_U->nex();
  const int Nc   = CommonParameters::Nc();
  const int Ndim = CommonParameters::Ndim();

  assert(force.nin() == Nin);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Nex);

  const int NinF  = m_fopr_prec->field_nin();
  const int NvolF = m_fopr_prec->field_nvol();
  const int NexF  = m_fopr_prec->field_nex();

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  m_fopr_prec->set_config(m_U);
  m_fopr->set_config(m_U);
  m_fopr_prec_force->set_config(m_U);
  m_fopr_force->set_config(m_U);

  Field v1(NinF, NvolF, NexF);
  m_fopr_prec->set_mode("H");
  m_fopr_prec->mult_dag(v1, m_psf);

  Field  v2(NinF, NvolF, NexF);
  int    Nconv;
  double diff;
  m_fprop_MD->set_config(m_U);
  m_fprop_MD->invert_DdagD(v2, v1, Nconv, diff);
  vout.general(m_vl, "    Solver: Nconv = %6d  diff  = %12.6e\n", Nconv, diff);

  m_fopr_force->force_core(force, v2);

  Field force_tmp(Nin, Nvol, Nex);
  m_fopr_prec_force->set_mode("Hdag");
  m_fopr_prec_force->force_core1(force_tmp, v2, m_psf);
  axpy(force, -1.0, force_tmp);

  m_fopr_prec_force->set_mode("H");
  m_fopr_prec_force->force_core1(force_tmp, m_psf, v2);
  axpy(force, -1.0, force_tmp);

  double Fave, Fmax, Fdev;
  force.stat(Fave, Fmax, Fdev);
  vout.general(m_vl, "    Fratio_ave = %12.6f  Fratio_max = %12.6f  Fratio_dev = %12.6f\n",
               Fave, Fmax, Fdev);
}


//====================================================================
//============================================================END=====
