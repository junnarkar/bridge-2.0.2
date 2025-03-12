/*!
        @file    projection_Stout_SU3.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#include "projection_Stout_SU3.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Projection_Stout_SU3::register_factory();
}
#endif

const std::string Projection_Stout_SU3::class_name = "Projection_Stout_SU3";

//====================================================================
void Projection_Stout_SU3::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }
}


//====================================================================
void Projection_Stout_SU3::get_parameters(Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Projection_Stout_SU3::init()
{
  m_flop = 0;
  m_time = 0.0;

  assert(CommonParameters::Nc() == NC);

  // strict check
  if (CommonParameters::Nc() != NC) {
    vout.crucial(m_vl, "Error at %s: Nc = 3 is needed, but Nc = %d\n", class_name.c_str(), CommonParameters::Nc());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Projection_Stout_SU3::print_stat()
{
  const double gflops = 1.0e-9 * double(m_flop) / m_time;

  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  total time: %f\n", m_time);
  vout.general(m_vl, "  total flop: %d\n", m_flop);
  vout.general(m_vl, "  GFlops    : %f\n", gflops);
}


//====================================================================
void Projection_Stout_SU3::project(Field_G& U,
                                   const double alpha,
                                   const Field_G& Cst, const Field_G& Uorg)
{
  //  int id = 31;
  //  KEK_FopCountStart(id);
  const double time0 = Communicator::get_time();

  // in stout projection, parameter alpha is dummy.

  const int Nex  = Uorg.nex();
  const int Nvol = Uorg.nvol();
  const int NinG = Uorg.nin();

  assert(Cst.nex() == Nex);
  assert(Cst.nvol() == Nvol);
  assert(U.nex() == Nex);
  assert(U.nvol() == Nvol);

  Mat_SU_N iQ0(NC);
  iQ0.unit();

  for (int mu = 0; mu < Nex; ++mu) {
    for (int site = 0; site < Nvol; ++site) {
      Mat_SU_N ut(NC);
      Uorg.mat(ut, site, mu);

      Mat_SU_N ct(NC);
      Cst.mat(ct, site, mu);

      Mat_SU_N iQ1(NC);
      iQ1.mult_nd(ct, ut);
      iQ1.at();

      Mat_SU_N iQ2(NC);
      iQ2.mult_nn(iQ1, iQ1);

      Mat_SU_N iQ3(NC);
      iQ3.mult_nn(iQ1, iQ2);

      Mat_SU_N e_iQ(NC);

      double norm = iQ1.norm2();
      if (norm > 1.0e-10) {
        double u, w;
        set_uw(u, w, iQ2, iQ3);

        dcomplex f0, f1, f2;
        set_fj(f0, f1, f2, u, w);

        for (int cc = 0; cc < NC * NC; ++cc) {
          dcomplex qt = f0 * cmplx(iQ0.r(cc), iQ0.i(cc))
                        + f1 * cmplx(iQ1.i(cc), -iQ1.r(cc))
                        - f2 * cmplx(iQ2.r(cc), iQ2.i(cc));
          e_iQ.set_r(cc, real(qt));
          e_iQ.set_i(cc, imag(qt));
        }
      } else {
        //  vout.general(m_vl,"project: |iQ1|^2 too small: %lf. Set e_iQ=1.\n",norm);
        e_iQ.unit();
      }

      Mat_SU_N ut2(NC);
      ut2.mult_nn(e_iQ, ut);
      U.set_mat(site, mu, ut2);
    }
  }

  /*
  unsigned long count;
  double time;
  KEK_FopCountFinish(id,&count,&time);
  m_time += time;
  m_flop += count;
  */
  const double time1 = Communicator::get_time();
  m_time += time1 - time0;
}


//====================================================================
void Projection_Stout_SU3::exp_iQ(Field_G& e_iQ, const Field_G& iQ)
{
  const int Nvol = iQ.nvol();
  const int Nex  = iQ.nex();

  Mat_SU_N iQ0(NC);

  iQ0.unit();

  for (int mu = 0; mu < Nex; ++mu) {
    for (int site = 0; site < Nvol; ++site) {
      Mat_SU_N iQ1 = iQ.mat(site, mu);
      Mat_SU_N iQ2 = iQ1 * iQ1;
      Mat_SU_N iQ3 = iQ1 * iQ2;

      double norm = iQ1.norm2();
      if (norm > 1.0e-10) {
        double u, w;
        set_uw(u, w, iQ2, iQ3);

        dcomplex f0, f1, f2;
        set_fj(f0, f1, f2, u, w);

        for (int cc = 0; cc < NC * NC; ++cc) {
          dcomplex qt = f0 * cmplx(iQ0.r(cc), iQ0.i(cc))
                        + f1 * cmplx(iQ1.i(cc), -iQ1.r(cc))
                        - f2 * cmplx(iQ2.r(cc), iQ2.i(cc));
          e_iQ.set_ri(cc, site, mu, real(qt), imag(qt));
        }
      } else {
        //      vout.general(m_vl,"exp_iQ: |iQ1|^2 too small: %lf. Set e_iQ=1.\n",norm);
        e_iQ.set_mat(site, mu, iQ0);
      }
    }
  }
}


//====================================================================

/*!
<ul>
<li>See the implementation note "note_cloverHMC.pdf" (21 Mar 2012) by H.Matsufuru.
<li>Evaluate \f$\Sigma_\mu'(x)\exp(iQ_\mu(x))+C_\mu^\dagger i\Lambda_\mu(x)\f$ and \f$i\Lambda_\mu(x)U_\mu(x)\f$ in eq.(93)
<li>argument Sigmap \f$=\Sigma_\mu'(x)\f$ in eq.(93) in terms of (k)-th smearing.
<li>argument Cst \f$=C_\mu(x)\f$ of eq.(43) in terms of (k-1)-th smearing.
<li>argument Uorg \f$=U_\mu(x)\f$ in (k-1)-th smearing.
<li>argument Xi is a resultant \f$\Sigma_\mu'(x)\exp(iQ_\mu(x))+C_\mu^\dagger i\Lambda_\mu(x)\f$.
<li>argument iTheta is a resultant \f$i\Lambda_\mu(x)U_\mu(x)\f$.
<li>Comment by [Y.Taniguchi 2012.04.16]
</ul>
*/
//====================================================================
void Projection_Stout_SU3::force_recursive(Field_G& Xi, Field_G& iTheta,
                                           const double alpha, const Field_G& Sigmap,
                                           const Field_G& Cst, const Field_G& Uorg)
{
  // in stout projection, parameter alpha is dummy.

  //  int id = 31;
  //  KEK_FopCountStart(id);
  const double time0 = Communicator::get_time();

  const int Nvol = CommonParameters::Nvol();
  const int Nex  = Xi.nex();

  assert(Xi.nvol() == Nvol);
  assert(iTheta.nvol() == Nvol);
  assert(Sigmap.nvol() == Nvol);
  assert(Cst.nvol() == Nvol);
  assert(Uorg.nvol() == Nvol);
  assert(iTheta.nex() == Nex);
  assert(Sigmap.nex() == Nex);
  assert(Cst.nex() == Nex);
  assert(Uorg.nex() == Nex);

  Mat_SU_N iQ0(NC);
  iQ0.unit();

  for (int mu = 0; mu < Nex; ++mu) {
    for (int site = 0; site < Nvol; ++site) {
      //! C_tmp \f$=C_\mu(x)\f$
      Mat_SU_N C_tmp(NC);
      Cst.mat(C_tmp, site, mu);

      //! U_tmp \f$=U_\mu(x)\f$
      Mat_SU_N U_tmp(NC);
      Uorg.mat(U_tmp, site, mu);

      // Sigmap_tmp \f$=\Sigma_\mu'(x)\f$
      Mat_SU_N Sigmap_tmp(NC);
      Sigmap.mat(Sigmap_tmp, site, mu);

      //! iQ1 \f$=iQ_\mu\f$
      Mat_SU_N iQ1(NC);
      iQ1.mult_nd(C_tmp, U_tmp);
      iQ1.at();

      Mat_SU_N iQ2(NC);
      iQ2.mult_nn(iQ1, iQ1);

      Mat_SU_N iQ3(NC);
      iQ3.mult_nn(iQ1, iQ2);

      // In order to aviod 1Q1=0
      Mat_SU_N e_iQ(NC), iGamma(NC);

      double norm = iQ1.norm2();
      if (norm > 1.0e-10) {
        double u, w;
        set_uw(u, w, iQ2, iQ3);

        dcomplex f0, f1, f2;
        set_fj(f0, f1, f2, u, w);

        for (int cc = 0; cc < NC * NC; ++cc) {
          dcomplex qt = f0 * cmplx(iQ0.r(cc), iQ0.i(cc))
                        + f1 * cmplx(iQ1.i(cc), -iQ1.r(cc))
                        - f2 * cmplx(iQ2.r(cc), iQ2.i(cc));
          e_iQ.set(cc, real(qt), imag(qt));
        }

        double xi0   = func_xi0(w);
        double xi1   = func_xi1(w);
        double u2    = u * u;
        double w2    = w * w;
        double cos_w = cos(w);

        dcomplex emiu = cmplx(cos(u), -sin(u));
        dcomplex e2iu = cmplx(cos(2.0 * u), sin(2.0 * u));

        dcomplex r01 = cmplx(2.0 * u, 2.0 * (u2 - w2)) * e2iu
                       + emiu * cmplx(16.0 * u * cos_w + 2.0 * u * (3.0 * u2 + w2) * xi0,
                                      -8.0 * u2 * cos_w + 2.0 * (9.0 * u2 + w2) * xi0);

        dcomplex r11 = cmplx(2.0, 4.0 * u) * e2iu
                       + emiu * cmplx(-2.0 * cos_w + (3.0 * u2 - w2) * xi0,
                                      2.0 * u * cos_w + 6.0 * u * xi0);

        dcomplex r21 = cmplx(0.0, 2.0) * e2iu
                       + emiu * cmplx(-3.0 * u * xi0, cos_w - 3.0 * xi0);

        dcomplex r02 = cmplx(-2.0, 0.0) * e2iu
                       + emiu * cmplx(-8.0 * u2 * xi0,
                                      2.0 * u * (cos_w + xi0 + 3.0 * u2 * xi1));

        dcomplex r12 = emiu * cmplx(2.0 * u * xi0,
                                    -cos_w - xi0 + 3.0 * u2 * xi1);

        dcomplex r22 = emiu * cmplx(xi0, -3.0 * u * xi1);

        double fden = 1.0 / (2 * (9.0 * u2 - w2) * (9.0 * u2 - w2));

        dcomplex b10 = cmplx(2.0 * u, 0.0) * r01 + cmplx(3.0 * u2 - w2, 0.0) * r02
                       - cmplx(30.0 * u2 + 2.0 * w2, 0.0) * f0;
        dcomplex b11 = cmplx(2.0 * u, 0.0) * r11 + cmplx(3.0 * u2 - w2, 0.0) * r12
                       - cmplx(30.0 * u2 + 2.0 * w2, 0.0) * f1;
        dcomplex b12 = cmplx(2.0 * u, 0.0) * r21 + cmplx(3.0 * u2 - w2, 0.0) * r22
                       - cmplx(30.0 * u2 + 2.0 * w2, 0.0) * f2;

        dcomplex b20 = r01 - cmplx(3.0 * u, 0.0) * r02 - cmplx(24.0 * u, 0.0) * f0;
        dcomplex b21 = r11 - cmplx(3.0 * u, 0.0) * r12 - cmplx(24.0 * u, 0.0) * f1;
        dcomplex b22 = r21 - cmplx(3.0 * u, 0.0) * r22 - cmplx(24.0 * u, 0.0) * f2;

        b10 *= cmplx(fden, 0.0);
        b11 *= cmplx(fden, 0.0);
        b12 *= cmplx(fden, 0.0);
        b20 *= cmplx(fden, 0.0);
        b21 *= cmplx(fden, 0.0);
        b22 *= cmplx(fden, 0.0);

        Mat_SU_N B1(NC), B2(NC);
        for (int cc = 0; cc < NC * NC; ++cc) {
          dcomplex qt1 = b10 * cmplx(iQ0.r(cc), iQ0.i(cc))
                         + b11 * cmplx(iQ1.i(cc), -iQ1.r(cc))
                         - b12 * cmplx(iQ2.r(cc), iQ2.i(cc));
          B1.set(cc, real(qt1), imag(qt1));

          dcomplex qt2 = b20 * cmplx(iQ0.r(cc), iQ0.i(cc))
                         + b21 * cmplx(iQ1.i(cc), -iQ1.r(cc))
                         - b22 * cmplx(iQ2.r(cc), iQ2.i(cc));
          B2.set(cc, real(qt2), imag(qt2));
        }

        Mat_SU_N USigmap(NC);
        USigmap.mult_nn(U_tmp, Sigmap_tmp);

        Mat_SU_N tmp1(NC);
        tmp1.mult_nn(USigmap, B1);

        Mat_SU_N tmp2(NC);
        tmp2.mult_nn(USigmap, B2);

        dcomplex tr1 = cmplx(tmp1.r(0) + tmp1.r(4) + tmp1.r(8),
                             tmp1.i(0) + tmp1.i(4) + tmp1.i(8));
        dcomplex tr2 = cmplx(tmp2.r(0) + tmp2.r(4) + tmp2.r(8),
                             tmp2.i(0) + tmp2.i(4) + tmp2.i(8));

        Mat_SU_N iQUS(NC);
        iQUS.mult_nn(iQ1, USigmap);

        Mat_SU_N iUSQ(NC);
        iUSQ.mult_nn(USigmap, iQ1);

        for (int cc = 0; cc < NC * NC; ++cc) {
          dcomplex qt = tr1 * cmplx(iQ1.i(cc), -iQ1.r(cc))
                        - tr2 * cmplx(iQ2.r(cc), iQ2.i(cc))
                        + f1 * cmplx(USigmap.r(cc), USigmap.i(cc))
                        + f2 * cmplx(iQUS.i(cc), -iQUS.r(cc))
                        + f2 * cmplx(iUSQ.i(cc), -iUSQ.r(cc));
          iGamma.set(cc, -imag(qt), real(qt));
        }
      } else {
        // vout.general(m_vl,"force_recursive: |iQ1|^2 too small: %lf. Set e_iQ=1.\n",norm);
        iGamma.zero();
        e_iQ.unit();
      }

      //! iGamma \f$=i\Lambda\f$
      iGamma.at();

      Mat_SU_N iTheta_tmp(NC);
      iTheta_tmp.mult_nn(iGamma, U_tmp);

      //! iTheta \f$=i\Lambda U_\mu(x)\f$
      iTheta.set_mat(site, mu, iTheta_tmp);

      Mat_SU_N Xi_tmp(NC);
      Xi_tmp.mult_nn(Sigmap_tmp, e_iQ);
      Xi_tmp.multadd_dn(C_tmp, iGamma);

      //! Xi \f$=\Sigma_\mu'(x)\exp(iQ_\mu(x))+C_\mu^\dagger i\Lambda_\mu(x)\f$.
      Xi.set_mat(site, mu, Xi_tmp);
    }
  }

  /*
  unsigned long count;
  double time;
  KEK_FopCountFinish(id,&count,&time);
  m_time += time;
  m_flop += count;
  */
  const double time1 = Communicator::get_time();
  m_time += time1 - time0;
}


//====================================================================
void Projection_Stout_SU3::set_fj(dcomplex& f0, dcomplex& f1, dcomplex& f2,
                                  const double& u, const double& w)
{
  const double xi0   = func_xi0(w);
  const double u2    = u * u;
  const double w2    = w * w;
  const double cos_w = cos(w);

  const double cos_u = cos(u);
  const double sin_u = sin(u);

  const dcomplex emiu = cmplx(cos_u, -sin_u);
  const dcomplex e2iu = cmplx(cos_u * cos_u - sin_u * sin_u, 2.0 * sin_u * cos_u);

  const dcomplex h0 = e2iu * cmplx(u2 - w2, 0.0)
                      + emiu * cmplx(8.0 * u2 * cos_w, 2.0 * u * (3.0 * u2 + w2) * xi0);
  const dcomplex h1 = cmplx(2 * u, 0.0) * e2iu
                      - emiu * cmplx(2.0 * u * cos_w, -(3.0 * u2 - w2) * xi0);
  const dcomplex h2 = e2iu - emiu * cmplx(cos_w, 3.0 * u * xi0);

  const double fden = 1.0 / (9.0 * u2 - w2);

  //- output
  f0 = h0 * fden;
  f1 = h1 * fden;
  f2 = h2 * fden;
}


//====================================================================
void Projection_Stout_SU3::set_uw(double& u, double& w,
                                  const Mat_SU_N& iQ2, const Mat_SU_N& iQ3)
{
  const double c0    = -(iQ3.i(0, 0) + iQ3.i(1, 1) + iQ3.i(2, 2)) / 3.0;
  const double c1    = -0.5 * (iQ2.r(0, 0) + iQ2.r(1, 1) + iQ2.r(2, 2));
  const double c13r  = sqrt(c1 / 3.0);
  const double c0max = 2.0 * c13r * c13r * c13r;

  const double theta = acos(c0 / c0max);

  //- output
  u = c13r * cos(theta / 3.0);
  w = sqrt(c1) * sin(theta / 3.0);
}


//====================================================================
double Projection_Stout_SU3::func_xi0(const double w)
{
  if (w == 0.0) {
    return 1.0;
  } else {
    return sin(w) / w;
  }
}


//====================================================================
double Projection_Stout_SU3::func_xi1(const double w)
{
  if (w < 0.25) {
    const double        w2 = w * w;
    const static double c0 = -1.0 / 3.0;
    const static double c1 = 1.0 / 30.0;
    const static double c2 = -1.0 / 840.0;
    const static double c3 = 1.0 / 45360.0;
    const static double c4 = -1.0 / 3991680.0;

    return c0 + w2 * (c1 + w2 * (c2 + w2 * (c3 + w2 * c4)));
  } else {
    return (w * cos(w) - sin(w)) / (w * w * w);
  }
}


//====================================================================
void Projection_Stout_SU3::exp_iQ_bf(Field_G& e_iQ, const Field_G& iQ)
{
  // brute force version of exponentiation: for check

  const static int Nprec = 32;

  const int Nvol = iQ.nvol();
  const int Nex  = iQ.nex();

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = 0; site < Nvol; ++site) {
      Mat_SU_N u0(NC);
      u0.unit();

      Mat_SU_N u1(NC);
      u1.unit();

      Mat_SU_N h1 = iQ.mat(site, ex);

      for (int iprec = 0; iprec < Nprec; ++iprec) {
        double   exf = 1.0 / (Nprec - iprec);
        Mat_SU_N u2  = h1 * u1;

        u2 *= exf;
        u1  = u2;
        u1 += u0;
      }

      u1.reunit();
      e_iQ.set_mat(site, ex, u1);
    }
  }
}


//====================================================================
//============================================================END=====
