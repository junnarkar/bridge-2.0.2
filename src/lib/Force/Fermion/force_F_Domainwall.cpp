/*!
        @file    force_F_Domainwall.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#include "force_F_Domainwall.h"

const std::string Force_F_Domainwall::class_name = "Force_F_Domainwall";

//====================================================================
void Force_F_Domainwall::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double           mq, M0, b, c;
  int              Ns;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("quark_mass", mq);
  err += params.fetch_double("domain_wall_height", M0);
  err += params.fetch_int("extent_of_5th_dimension", Ns);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  int err2 = 0;
  err2 += params.fetch_double("coefficient_b", b);
  err2 += params.fetch_double("coefficient_c", c);

  if (err2) {
    vout.general(m_vl, "  coefficients b, c are not provided:"
                       " set to Shamir's form.\n");
    b = 1.0;
    c = 0.0;
  }

  set_parameters(mq, M0, Ns, bc, b, c);
}


//====================================================================
void Force_F_Domainwall::get_parameters(Parameters& params) const
{
  params.set_double("quark_mass", m_mq);
  params.set_double("domain_wall_height", m_M0);
  params.set_int("extent_of_5th_dimension", m_Ns);
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_double_vector("coefficient_b", m_b);
  params.set_double_vector("coefficient_c", m_c);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Force_F_Domainwall::set_parameters(const double mq,
                                        const double M0,
                                        const int Ns,
                                        const std::vector<int> bc,
                                        const double b,
                                        const double c)
{
  const int Ndim = CommonParameters::Ndim();

  //- range check
  int err = 0;
  err += ParameterCheck::non_zero(mq);
  err += ParameterCheck::non_zero(M0);
  err += ParameterCheck::non_zero(Ns);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  assert(bc.size() == Ndim);

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    m_mq = mq;
    m_M0 = M0;
    m_Ns = Ns;

    m_boundary.resize(Ndim);
    m_boundary = bc;

    if (m_b.size() != m_Ns) {
      m_b.resize(m_Ns);
      m_c.resize(m_Ns);
    }
    for (int is = 0; is < m_Ns; ++is) {
      m_b[is] = b;
      m_c[is] = c;
    }
  }

  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  mq   = %12.8f\n", mq);
  vout.general(m_vl, "  M0   = %12.8f\n", M0);
  vout.general(m_vl, "  Ns   = %4d\n", Ns);
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
  }
  vout.general(m_vl, "  coefficients:\n");
  for (int is = 0; is < m_Ns; ++is) {
    vout.general(m_vl, "  b[%2d] = %16.10f  c[%2d] = %16.10f\n",
                 is, m_b[is], is, m_c[is]);
  }

  //-- post-process
  //- Domain-wall operator
  m_fopr_dw->set_parameters(m_mq, m_M0, m_Ns, m_boundary, b, c);

  //- Wilson opr/force
  const double kappa = 1.0 / (8.0 - 2.0 * m_M0);
  m_fopr_w->set_parameters(kappa, m_boundary);
  m_force_w->set_parameters(kappa, m_boundary);
}


//====================================================================
void Force_F_Domainwall::force_core1(Field& force, const Field& zeta,
                                     const Field& eta)
{
  vout.crucial("Error at %s: unimplemented.\n", __func__);
  exit(EXIT_FAILURE);
}


//====================================================================
void Force_F_Domainwall::force_udiv(Field& force_, const Field& eta)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  Field_F zeta(Nvol, m_Ns);

  m_fopr_dw->set_mode("H");
  m_fopr_dw->mult(zeta, eta);

  Field_G force(Nvol, Ndim);
  force_udiv1(force, zeta, eta);
  copy(force_, force); // force_ = force;

  force_udiv1(force, eta, zeta);
  axpy(force_, 1.0, force); // force_ += force;
}


//====================================================================
void Force_F_Domainwall::force_udiv1(Field& force_, const Field& zeta, const Field& eta)
{
  const int Nc   = CommonParameters::Nc();
  const int Nd   = CommonParameters::Nd();
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();
  const int NinF = 2 * Nc * Nd;

  assert(eta.nin() == NinF);
  assert(eta.nvol() == Nvol);
  assert(eta.nex() == m_Ns);

  assert(zeta.nin() == NinF);
  assert(zeta.nvol() == Nvol);
  assert(zeta.nex() == m_Ns);

  force_.set(0.0);

  for (int is = 0; is < m_Ns; ++is) {
    Field_F zeta4;
    zeta4.setpart_ex(0, zeta, m_Ns - 1 - is);

    Field_F eta4;
    eta4.setpart_ex(0, eta, is);

    scal(eta4, 4.0 - m_M0);

    Field_G force(Nvol, Ndim);
    m_force_w->force_udiv1(force, zeta4, eta4);
    axpy(force_, 1.0, force); // force_ += force;
  }
}


//====================================================================
//============================================================END=====
