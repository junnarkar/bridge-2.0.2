/*!
      @file    action_F_Ratio_lex_alt-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-02-28 16:09:41 #$
      @version $LastChangedRevision: 2492 $
*/

template<typename AFIELD>
const std::string Action_F_Ratio_lex_alt<AFIELD>::class_name
  = "Action_F_Ratio_lex_alt";

//====================================================================
template<typename AFIELD>
void Action_F_Ratio_lex_alt<AFIELD>::init()
{
  m_vl = CommonParameters::Vlevel();
}


//====================================================================
template<typename AFIELD>
void Action_F_Ratio_lex_alt<AFIELD>::tidyup()
{
  // do nothing.
}


//====================================================================
template<typename AFIELD>
void Action_F_Ratio_lex_alt<AFIELD>::set_parameters(const Parameters& params)
{
  Bridge::VerboseLevel vl = params.get_VerboseLevel();
  set_parameter_verboselevel(vl);
}


//====================================================================
template<typename AFIELD>
void Action_F_Ratio_lex_alt<AFIELD>::set_parameters()
{
  vout.general(m_vl, "%s:\n", class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void Action_F_Ratio_lex_alt<AFIELD>::set_config(Field *U)
{
  m_U = U;

  m_fopr_prec->set_config(U);
  m_fopr_prec_force->set_config(U);

  m_fopr->set_config(U);
  m_fopr_force->set_config(U);
}


//====================================================================
template<typename AFIELD>
double Action_F_Ratio_lex_alt<AFIELD>::langevin(RandomNumbers *rand)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  const int NinF     = m_fopr_prec->field_nin();
  const int NvolF    = m_fopr_prec->field_nvol();
  const int NexF     = m_fopr_prec->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  assert(NvolF == Nvol);
  m_psf.reset(NinF, NvolF, NexF);

  Timer  timer;
  double elapsed_time;
  timer.start();

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  Field xi(NinF, NvolF, NexF);
  rand->gauss_lex_global(xi);

  AFIELD xiA(NinF, NvolF, NexF);

  Index_lex_alt<real_t, AFIELD::IMPL> index_alt;

#pragma omp parallel
  {
    if (m_fopr->needs_convert()) {
      m_fopr->convert(xiA, xi);
    } else {
      convert(index_alt, xiA, xi);
    }
  }

  m_fopr_prec->set_config(m_U);
  m_fopr->set_config(m_U);

  AFIELD v2(NinF, NvolF, NexF);
  m_fopr->set_mode("H");
  m_fopr->mult_dag(v2, xiA);

  AFIELD v1(NinF, NvolF, NexF);
  m_fopr_prec->set_mode("H");
  m_fopr_prec->mult_dag(v1, v2);

  int    Nconv;
  double diff;
  m_fprop_H_prec->set_config(m_U);
  m_fprop_H_prec->set_mode("DdagD");
  m_fprop_H_prec->invert(m_psf, v1, Nconv, diff);
  vout.general(m_vl, "    Fprop_H_prec: %6d %18.15e\n", Nconv, diff);

  const double xi2   = xi.norm();
  const double H_psf = xi2 * xi2;

  timer.stop();
  elapsed_time = timer.elapsed_sec();

  vout.general(m_vl, "    H_Fratio     = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);
  vout.general(m_vl, "    elapsed time: %12.6f [sec]\n", elapsed_time);

  return H_psf;
}


//====================================================================
template<typename AFIELD>
double Action_F_Ratio_lex_alt<AFIELD>::calcH()
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  const int NinF     = m_fopr_prec->field_nin();
  const int NvolF    = m_fopr_prec->field_nvol();
  const int NexF     = m_fopr_prec->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  Timer  timer;
  double elapsed_time;
  timer.start();

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  m_fopr_prec->set_config(m_U);
  m_fopr->set_config(m_U);

  AFIELD v1(NinF, NvolF, NexF);
  m_fopr_prec->set_mode("H");
  m_fopr_prec->mult(v1, m_psf);

  AFIELD v2(NinF, NvolF, NexF);
  int    Nconv;
  double diff;

  m_fprop_H->set_config(m_U);
  m_fprop_H->reset_performance();
  m_fprop_H->set_mode("DdagD");
  m_fprop_H->invert(v2, v1, Nconv, diff);

  double flop_count, elapsed_time2, gflops;
  m_fprop_H->get_performance(flop_count, elapsed_time2);
  if (elapsed_time2 > 0.0) {
    gflops = (flop_count / elapsed_time2) * 1.0e-9;
  } else {
    gflops = 0.0;
  }

  vout.general(m_vl, "    Fprop_H: %6d  %12.4e  %12.4e [GFlops]\n",
               Nconv, diff, gflops);

  double H_psf = dot(v1, v2);

  timer.stop();
  elapsed_time = timer.elapsed_sec();

  vout.general(m_vl, "    H_Fratio     = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);
  vout.general(m_vl, "    elapsed time: %12.6f [sec]\n", elapsed_time);

  return H_psf;
}


//====================================================================
template<typename AFIELD>
void Action_F_Ratio_lex_alt<AFIELD>::force(Field& force)
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

  int NinF  = m_fopr_prec->field_nin();
  int NvolF = m_fopr_prec->field_nvol();
  int NexF  = m_fopr_prec->field_nex();

  Timer  timer;
  double elapsed_time;
  timer.start();

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  m_fopr_prec->set_config(m_U);
  m_fopr->set_config(m_U);
  m_fopr_prec_force->set_config(m_U);
  m_fopr_force->set_config(m_U);

  AFIELD v1(NinF, NvolF, NexF);
  m_fopr_prec->set_mode("H");
#pragma omp parallel
  {
    m_fopr_prec->mult(v1, m_psf);
  }

  AFIELD v2(NinF, NvolF, NexF);
  int    Nconv;
  double diff;
  m_fprop_MD->set_config(m_U);
  m_fprop_MD->set_mode("DdagD");
  m_fprop_MD->invert(v2, v1, Nconv, diff);
  vout.general(m_vl, "    Fprop_MD: %6d %18.15e\n", Nconv, diff);

  m_fopr_force->force_core(force1, v2);

  AFIELD force_tmp(Nin, Nvol, Nex);
  m_fopr_prec_force->set_mode("H");
  m_fopr_prec_force->force_core1(force_tmp, v2, m_psf);

#pragma omp parallel
  {
    axpy(force1, -1.0, force_tmp);
  }

  m_fopr_prec_force->set_mode("Hdag");
  m_fopr_prec_force->force_core1(force_tmp, m_psf, v2);

#pragma omp parallel
  {
    axpy(force1, -1.0, force_tmp);
  }

  Index_lex_alt<real_t, AFIELD::IMPL> index_lex;
#pragma omp parallel
  {
    reverse_gauge(index_lex, force, force1);
  }

  double Fave, Fmax, Fdev;
  force.stat(Fave, Fmax, Fdev);

  timer.stop();
  elapsed_time = timer.elapsed_sec();

  vout.general(m_vl, "    Fave = %12.6f  Fmax = %12.6f  Fdev = %12.6f\n",
               Fave, Fmax, Fdev);
  vout.general(m_vl, "    elapsed time: %12.6f [sec]\n", elapsed_time);

  vout.paranoiac(m_vl, "    Fratio::force end\n");
}


//============================================================END=====
