/*!
        @file    force_F_CloverTerm.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#include "force_F_CloverTerm.h"

const std::string Force_F_CloverTerm::class_name = "Force_F_CloverTerm";

//====================================================================
void Force_F_CloverTerm::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double           kappa, cSW;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_double("clover_coefficient", cSW);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(kappa, cSW, bc);
}


//====================================================================
void Force_F_CloverTerm::get_parameters(Parameters& params) const
{
  params.set_double("hopping_parameter", m_kappa);
  params.set_double("clover_coefficient", m_cSW);
  params.set_int_vector("boundary_condition", m_boundary);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Force_F_CloverTerm::set_parameters(const double kappa, const double cSW, const std::vector<int> bc)
{
  const int Ndim = CommonParameters::Ndim();

  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", kappa);
  vout.general(m_vl, "  cSW   = %12.8f\n", cSW);
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
  }

  //- range check
  // NB. kappa,cSW == 0 is allowed.
  assert(bc.size() == Ndim);

  //- store values
  m_kappa = kappa;
  m_cSW   = cSW;

  m_boundary.resize(Ndim);
  m_boundary = bc;
}


//====================================================================
void Force_F_CloverTerm::init(const std::string repr)
{
  m_repr = repr;

  m_Ndim = CommonParameters::Ndim();
  const int Nvol = CommonParameters::Nvol();

  m_fopr_csw = new Fopr_CloverTerm(m_repr);
  m_Cud      = new Field_G(Nvol, m_Ndim * m_Ndim);

  m_boundary.resize(CommonParameters::Ndim());
}


//====================================================================
void Force_F_CloverTerm::tidyup()
{
  delete m_Cud;
  delete m_fopr_csw;
}


// override default force_core()
//====================================================================
void Force_F_CloverTerm::force_core(Field& force, const Field& eta)
{
  vout.crucial(m_vl, "Error at %s: force_core() is not available.\n", class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
void Force_F_CloverTerm::force_udiv(Field& force, const Field& eta)
{
  vout.crucial(m_vl, "Error at %s: force_udiv() is not available.\n", class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
void Force_F_CloverTerm::force_udiv1(Field& force_,
                                     const Field& zeta_, const Field& eta_)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  Field_G force(Nvol, Ndim);
  Field_F zeta(zeta_);
  Field_F eta(eta_);

  force_udiv1_impl(force, zeta, eta);

  copy(force_, force); // force_ = force;
}


//====================================================================
void Force_F_CloverTerm::force_udiv1_impl(Field_G& force,
                                          const Field_F& zeta, const Field_F& eta)
{
  const int Nc   = CommonParameters::Nc();
  const int Nd   = CommonParameters::Nd();
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  force.set(0.0);

  Field_F eta2;
  m_fopr_csw->mult_gm5(eta2, eta);

  for (int mu = 0; mu < Ndim; ++mu) {
    for (int nu = 0; nu < Ndim; ++nu) {
      if (nu == mu) continue;

      Field_F eta3;
      m_fopr_csw->mult_isigma(eta3, eta2, mu, nu);

      Field_G Umu;
      Umu.setpart_ex(0, *m_U, mu);

      Field_G Unu;
      Unu.setpart_ex(0, *m_U, nu);

      int ex = 0;

      // R(1) and R(5)
      Field_F vt1;
      mult_Field_Gd(vt1, 0, *m_Cud, index_dir(mu, nu), eta3, ex);

      Field_G force1;
      tensorProd_Field_F(force1, zeta, vt1);

      Field_G force2;
      copy(force2, force1); // force2 = force1;

      // R(2)
      Field_F vt3;
      mult_Field_Gd(vt3, 0, Umu, 0, eta3, ex);

      ShiftField_lex shift;
      shift.backward(vt1, vt3, nu);

      Field_F vt2;
      shift.backward(vt2, zeta, nu);

      Field_G Utmp;
      shift.backward(Utmp, Unu, mu);

      mult_Field_Gn(vt3, 0, Utmp, 0, vt1, ex);

      Field_F vt4;
      mult_Field_Gn(vt4, 0, Unu, 0, vt2, ex);
      tensorProd_Field_F(force1, vt4, vt3);
      axpy(force2, 1.0, force1); // force2 += force1;

      // R(4) and R(8)
      shift.backward(vt1, eta3, mu);

      Field_F zeta_mu;
      shift.backward(zeta_mu, zeta, mu);
      mult_Field_Gn(vt4, 0, *m_Cud, index_dir(mu, nu), zeta_mu, ex);
      tensorProd_Field_F(force1, vt4, vt1);
      axpy(force2, 1.0, force1); // force2 += force1;

      // R(3)
      shift.backward(vt1, eta3, nu);
      mult_Field_Gn(vt3, 0, Unu, 0, vt1, ex);
      mult_Field_Gn(vt4, 0, Umu, 0, zeta_mu, ex);
      shift.backward(vt1, vt3, mu);
      shift.backward(vt2, vt4, nu);
      mult_Field_Gn(vt4, 0, Unu, 0, vt2, ex);
      tensorProd_Field_F(force1, vt4, vt1);
      axpy(force2, 1.0, force1); // force2 += force1;

      // R(6)
      shift.backward(Utmp, Unu, mu);

      Field_G Utmp2;
      mult_Field_Gdd(Utmp2, 0, Utmp, 0, Umu, 0);
      mult_Field_Gn(vt1, 0, Utmp2, 0, eta3, ex);
      mult_Field_Gd(vt2, 0, Unu, 0, zeta, ex);
      shift.forward(vt3, vt1, nu);
      shift.forward(vt4, vt2, nu);
      tensorProd_Field_F(force1, vt4, vt3);
      axpy(force2, -1.0, force1); // force2 -= force1;

      // R(7)
      mult_Field_Gd(vt1, 0, Unu, 0, eta3, ex);
      mult_Field_Gn(vt2, 0, Umu, 0, zeta_mu, ex);
      shift.backward(vt3, vt1, mu);
      shift.forward(vt1, vt3, nu);
      mult_Field_Gd(vt4, 0, Unu, 0, vt2, ex);
      shift.forward(vt2, vt4, nu);
      tensorProd_Field_F(force1, vt2, vt1);
      axpy(force2, -1.0, force1);           // force2 -= force1;

      scal(force2, -m_kappa * m_cSW / 8.0); // force2 *= -m_kappa * m_cSW / 8.0;
      force.addpart_ex(mu, force2, 0);
    }
  }
}


//====================================================================
void Force_F_CloverTerm::set_component()
{
  const int Nc   = CommonParameters::Nc();
  const int Nd   = CommonParameters::Nd();
  const int Nvol = CommonParameters::Nvol();

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;

      Staple_lex staple;

      Field_G Cmu_ud1;
      staple.upper(Cmu_ud1, *m_U, mu, nu);

      Field_G Cmu_ud2;
      staple.lower(Cmu_ud2, *m_U, mu, nu);

      axpy(Cmu_ud1, -1.0, Cmu_ud2);
      m_Cud->setpart_ex(index_dir(mu, nu), Cmu_ud1, 0);
    }
  }
}


//====================================================================
//============================================================END=====
