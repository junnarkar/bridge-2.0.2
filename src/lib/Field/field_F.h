/*!
        @file    field_F.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#ifndef FIELD_F_INCLUDED
#define FIELD_F_INCLUDED

#include "field_G.h"
#include "Tools/vec_SU_N.h"
#include "Tools/gammaMatrix.h"
#include "ResourceManager/threadManager.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Wilson-type fermion field.

/*!
   This class defines 4-spinor (in the case of Ndim=4) fermion
   field, which is mainly used by Wilson-type fermions.
   Original version of this class was written by J.Noaki.
   H.Matsufuru added several functions and modified intefaces
   of several functionality.
                                     [28 Dec 2011 H.Matsufuru]
   mult_GMproj2 is generalized for Wilson_General.
                                     [21 Mar 2015 Y.Namekawa]
 */
class Field_F : public Field
{
 private:
  int m_Nc;   // num of the color elements
  int m_Nc2;  // num of the double color elements
  int m_Nd;   // num of the spinor elements

  inline
  int myindex(const int c2, const int s, const int site, const int ex)
  const
  {
    return Field::myindex(c2 + m_Nc2 * s, site, ex);
  }

 public:

  explicit
  Field_F(const int Nvol = CommonParameters::Nvol(), const int Nex = 1) :
    Field(2 * CommonParameters::Nc() * CommonParameters::Nd(),
          Nvol, Nex, Element_type::COMPLEX),
    m_Nc(CommonParameters::Nc()),
    m_Nc2(2 * CommonParameters::Nc()),
    m_Nd(CommonParameters::Nd())
  {
    check();
  }

  Field_F clone() const
  {
    return Field_F(nvol(), nex());
  }

  // conversion from Field type

  Field_F(const Field& x) :
    Field(x),
    m_Nc(CommonParameters::Nc()),
    m_Nc2(2 * CommonParameters::Nc()),
    m_Nd(CommonParameters::Nd())
  {
    check();
  }

  void reset(int Nvol, int Nex)
  {
    Field::reset(m_Nc2 * m_Nd, Nvol, Nex);
  }

  // assignment
  //Field_F& operator=(const double a) { field = a; return *this; }
  Field_F& operator=(const Field_F& v) { copy(*this, v); return *this; }

  int nc() const { return m_Nc; }
  int nc2() const { return m_Nc2; }
  int nd() const { return m_Nd; }

  // accessors
  double cmp_r(const int cc, const int s, const int site, const int e = 0)
  const
  {
    return field[myindex(2 * cc, s, site, e)];
  }

  double cmp_i(const int cc, const int s, const int site, const int e = 0)
  const
  {
    return field[myindex(2 * cc + 1, s, site, e)];
  }

  void set_r(const int cc, const int s, const int site, const int e, const double re)
  {
    field[myindex(2 * cc, s, site, e)] = re;
  }

  void set_i(const int cc, const int s, const int site, const int e, const double im)
  {
    field[myindex(2 * cc + 1, s, site, e)] = im;
  }

  void set_ri(const int cc, const int s, const int site, const int e, const double re, const double im)
  {
    field[myindex(2 * cc, s, site, e)]     = re;
    field[myindex(2 * cc + 1, s, site, e)] = im;
  }

  Vec_SU_N vec(const int s, const int site, const int e = 0) const
  {
    Vec_SU_N Tmp;

    for (int cc = 0; cc < m_Nc; ++cc) {
      Tmp.set(cc,
              field[myindex(2 * cc, s, site, e)],
              field[myindex(2 * cc + 1, s, site, e)]);
    }
    return Tmp;
  }

  void set_vec(const int s, const int site, const int e, const Vec_SU_N& F)
  {
    for (int cc = 0; cc < m_Nc; ++cc) {
      field[myindex(2 * cc, s, site, e)]     = F.r(cc);
      field[myindex(2 * cc + 1, s, site, e)] = F.i(cc);
    }
  }

  void add_vec(const int s, const int site, const int e, const Vec_SU_N& F)
  {
    for (int cc = 0; cc < m_Nc; ++cc) {
      field[myindex(2 * cc, s, site, e)]     += F.r(cc);
      field[myindex(2 * cc + 1, s, site, e)] += F.i(cc);
    }
  }

  void clear_vec(const int s, const int site, const int e)
  {
    for (int cc = 0; cc < m_Nc2; ++cc) {
      field[myindex(cc, s, site, e)] = 0.0;
    }
  }

 private:

  //! check several assumptions for performance implementation.
  void check();
};


//- function style
void mult_Field_Gn(Field_F& y, const int ex,
                   const Field_G& u, int ex1,
                   const Field_F& x, int ex2);

void mult_Field_Gd(Field_F& y, const int ex,
                   const Field_G& u, int ex1,
                   const Field_F& x, int ex2);

void multadd_Field_Gn(Field_F& y, const int ex,
                      const Field_G& u, int ex1,
                      const Field_F& x, int ex2,
                      const double a);

void multadd_Field_Gd(Field_F& y, const int ex,
                      const Field_G& u, int ex1,
                      const Field_F& x, int ex2,
                      const double a);

//! gamma matrix multiplication
void mult_GM(Field_F& y,
             const GammaMatrix& gm,
             const Field_F& x);

//! gamma matrix multiplication (i is multiplied)
void mult_iGM(Field_F& y,
              const GammaMatrix& gm,
              const Field_F& x);

//! projection with gamma matrix: (1 \pm gamma)/2
void mult_GMproj(Field_F& y,
                 const int pm, const GammaMatrix& gm,
                 const Field_F& x);

//! projection with gamma matrix: (1 \pm gamma)
void mult_GMproj2(Field_F& y,
                  const int pm, const GammaMatrix& gm,
                  const Field_F& x);

//! projection with gamma matrix: (r_s \pm nu_s * gamma)
void mult_GMproj2(Field_F& y,
                  const double nu_s,
                  const int pm,
                  const double r_s,
                  const GammaMatrix& gm,
                  const Field_F& x);

#endif
