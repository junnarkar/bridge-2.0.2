/*!
        @file    forceSmear_APE.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#include "forceSmear_APE.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = ForceSmear_APE::register_factory();
}
#endif

const std::string ForceSmear_APE::class_name = "ForceSmear_APE";

//====================================================================
void ForceSmear_APE::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double rho1;

  int err = 0;
  err += params.fetch_double("rho_uniform", rho1);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(rho1);
}


//====================================================================
void ForceSmear_APE::get_parameters(Parameters& params) const
{
  params.set_double_vector("rho", m_rho);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void ForceSmear_APE::set_parameters(const double rho1)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  rho = %8.4f\n", rho1);

  //- range check
  // NB. rho == 0 is allowed.

  //- store values
  // m_rho.resize(m_Ndim * m_Ndim);  // already resized in init.
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      m_rho[mu + nu * m_Ndim] = rho1;
    }
  }
}


//====================================================================
void ForceSmear_APE::set_parameters(const std::vector<double>& rho)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  rho[%d] = %8.4f\n", mu, rho[mu]);
  }

  //- range check
  // NB. rho == 0 is allowed.
  assert(rho.size() == m_Ndim * m_Ndim);

  //- store values
  // m_rho.resize(m_Ndim * m_Ndim);  // already resized in init.
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      m_rho[mu + nu * m_Ndim] = rho[mu + nu * m_Ndim];
    }
  }
}


//====================================================================
void ForceSmear_APE::init()
{
  m_Ndim = CommonParameters::Ndim();
  m_Nvol = CommonParameters::Nvol();

  m_rho.resize(m_Ndim * m_Ndim);
  m_U.resize(m_Ndim);
  m_iTheta.resize(m_Ndim);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      m_rho[mu + nu * m_Ndim] = 0.0;
    }
  }
}


//====================================================================
void ForceSmear_APE::force_udiv(Field_G& Sigma, const Field_G& Sigmap, const Field_G& U)
{
  const int Nc = CommonParameters::Nc();

  assert(Sigmap.nin() == (2 * Nc * Nc));
  assert(Sigmap.nvol() == m_Nvol);
  assert(Sigmap.nex() == m_Ndim);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_U[mu].setpart_ex(0, U, mu);
  }

  for (int mu = 0; mu < m_Ndim; ++mu) {
    Field_G C;
    C.set(0.0);

    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;
      double rho = m_rho[mu + m_Ndim * nu];

      Field_G c_tmp;
      staple(c_tmp, m_U[mu], m_U[nu], mu, nu);
      C.addpart_ex(0, c_tmp, 0, rho);
    }

    Field_G sigmap_tmp;
    sigmap_tmp.setpart_ex(0, Sigmap, mu);

    double alpha = m_rho[mu + m_Ndim * mu];

    Field_G Xi;
    m_proj->force_recursive(Xi, m_iTheta[mu],
                            alpha, sigmap_tmp, C, m_U[mu]);
    Sigma.setpart_ex(mu, Xi, 0);
  }

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;
      double rho = m_rho[mu + m_Ndim * nu];

      Field_G sigma_tmp;
      force_each(sigma_tmp, m_U[mu], m_U[nu],
                 m_iTheta[mu], m_iTheta[nu], mu, nu);
      Sigma.addpart_ex(mu, sigma_tmp, 0, rho);
    }
  }
}


//====================================================================
void ForceSmear_APE::force_each(Field_G& Sigma_mu,
                                const Field_G& V_mu, const Field_G& V_nu,
                                const Field_G& iTheta_mu,
                                const Field_G& iTheta_nu,
                                const int mu, const int nu)
{
  Sigma_mu.set(0.0);

  Field_G vt1;
  m_shift.backward(vt1, V_nu, mu);

  Field_G vt2;
  m_shift.backward(vt2, V_mu, nu);

  Field_G vt3;
  mult_Field_Gnd(vt3, 0, vt1, 0, vt2, 0);
  multadd_Field_Gnd(Sigma_mu, 0, vt3, 0, iTheta_nu, 0, 1.0);

  mult_Field_Gdn(vt3, 0, iTheta_mu, 0, V_nu, 0);
  mult_Field_Gdn(vt2, 0, vt1, 0, vt3, 0);
  m_shift.forward(vt3, vt2, nu);
  axpy(Sigma_mu, 1.0, vt3); // Sigma_mu += vt3;

  mult_Field_Gdn(vt3, 0, V_mu, 0, iTheta_nu, 0);
  mult_Field_Gdn(vt2, 0, vt1, 0, vt3, 0);
  m_shift.forward(vt3, vt2, nu);
  axpy(Sigma_mu, 1.0, vt3); // Sigma_mu += vt3;

  m_shift.backward(vt1, iTheta_nu, mu);
  m_shift.backward(vt2, V_mu, nu);
  mult_Field_Gnd(vt3, 0, vt1, 0, vt2, 0);
  multadd_Field_Gnd(Sigma_mu, 0, vt3, 0, V_nu, 0, 1.0);

  mult_Field_Gdd(vt2, 0, vt1, 0, V_mu, 0);
  mult_Field_Gnn(vt3, 0, vt2, 0, V_nu, 0);
  m_shift.forward(vt2, vt3, nu);
  axpy(Sigma_mu, 1.0, vt2); // Sigma_mu += vt2;

  m_shift.backward(vt1, V_nu, mu);
  m_shift.backward(vt2, iTheta_mu, nu);
  mult_Field_Gnd(vt3, 0, vt1, 0, vt2, 0);
  multadd_Field_Gnd(Sigma_mu, 0, vt3, 0, V_nu, 0, 1.0);
}


//====================================================================
void ForceSmear_APE::staple(Field_G& c,
                            const Field_G& u_mu, const Field_G& u_nu,
                            const int mu, const int nu)
{
  //- upper direction
  Field_G v1;

  m_shift.backward(v1, u_mu, nu);

  Field_G v2;
  mult_Field_Gnn(v2, 0, u_nu, 0, v1, 0);

  m_shift.backward(v1, u_nu, mu);
  mult_Field_Gnd(c, 0, v2, 0, v1, 0);

  //- lower direction
  m_shift.backward(v2, u_nu, mu);
  mult_Field_Gnn(v1, 0, u_mu, 0, v2, 0);
  mult_Field_Gdn(v2, 0, u_nu, 0, v1, 0);
  m_shift.forward(v1, v2, nu);
  c.addpart_ex(0, v1, 0);
}


//====================================================================
//============================================================END=====
