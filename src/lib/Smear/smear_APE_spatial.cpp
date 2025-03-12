/*!
        @file    smear_APE_spatial.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#include "smear_APE_spatial.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Smear_APE_spatial::register_factory();
}
#endif

const std::string Smear_APE_spatial::class_name = "Smear_APE_spatial";

//====================================================================
void Smear_APE_spatial::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double rho;

  int err = 0;
  err += params.fetch_double("rho", rho);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(rho);
}


//====================================================================
void Smear_APE_spatial::get_parameters(Parameters& params) const
{
  params.set_double("rho", m_rho);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Smear_APE_spatial::set_parameters(const double rho)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  rho = %10.6F\n", rho);

  //- range check
  // NB. rho == 0 is allowed.

  //- store values
  m_rho = rho;
}


//====================================================================
void Smear_APE_spatial::smear(Field_G& Usmear, const Field_G& U)
{
  const int Nvol = CommonParameters::Nvol();

  assert(U.nvol() == Nvol);
  assert(U.nex() == m_Ndim);

  assert(Usmear.nvol() == Nvol);
  assert(Usmear.nex() == m_Ndim);

  const int Ndim_spc = m_Ndim - 1;

  Staple_lex staple;

  double plaq = staple.plaq_s(U);
  vout.general(m_vl, "  plaq_s(org  ) = %12.8f\n", plaq);

  plaq = staple.plaq_t(U);
  vout.general(m_vl, "  plaq_t(org  ) = %12.8f\n", plaq);

  Usmear.set(0.0);

  for (int mu = 0; mu < Ndim_spc; ++mu) {
    Field_G c_tmp;
    c_tmp.set(0.0);

    Field_G u_tmp;
    u_tmp.setpart_ex(0, U, mu);

    Field_G u_tmp2;

    for (int nu = 0; nu < Ndim_spc; ++nu) {
      if (nu != mu) {
        staple.upper(u_tmp2, U, mu, nu);
        c_tmp.addpart_ex(0, u_tmp2, 0, m_rho);

        staple.lower(u_tmp2, U, mu, nu);
        c_tmp.addpart_ex(0, u_tmp2, 0, m_rho);
      }
    }

    m_proj->project(u_tmp2, m_rho, c_tmp, u_tmp);
    Usmear.setpart_ex(mu, u_tmp2, 0);
  }

  const int mu = m_Ndim - 1; // temporal link: unsmeared.
  Usmear.setpart_ex(mu, U, mu);

  plaq = staple.plaq_s(Usmear);
  vout.general(m_vl, "  plaq_s(smear) = %12.8f\n", plaq);

  plaq = staple.plaq_t(Usmear);
  vout.general(m_vl, "  plaq_t(smear) = %12.8f\n", plaq);
}


//====================================================================
//============================================================END=====
