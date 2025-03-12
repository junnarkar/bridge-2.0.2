/*!
        @file    force_F_Staggered.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$
        @version $LastChangedRevision: 2492 $
*/

#ifndef FORCE_F_STAGGERED_INCLUDED
#define FORCE_F_STAGGERED_INCLUDED

#include "Force/Fermion/force_F.h"
#include "Fopr/fopr_Staggered.h"
#include "Field/field_F_1spinor.h"
#include "Field/shiftField_lex.h"

//! Force calculation for even-odd staggered fermion operator.

/*!
    This is a temporary implementation for Complex Langevin version.
                                [20 Sep 2015 H.Matsufuru]
 */

class Force_F_Staggered : public Force
{
 public:
  static const std::string class_name;

 private:

  double m_mq;                   //!< quark mass.
  std::vector<int> m_boundary;   //!< boundary conditions.
  Bridge::VerboseLevel m_vl;

  Fopr_Staggered *m_fopr_ks;
  Field_F_1spinor m_psf;
  ShiftField_lex m_shift;

 public:

  Force_F_Staggered() : Force() { init(); }

  ~Force_F_Staggered()
  {
    delete m_fopr_ks;
  }

  void set_parameters(const Parameters& params);

  void set_parameters(double mq, const std::vector<int> bc);

  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
    m_fopr_ks->set_config(U);
  }

  void force_udiv(Field& force, const Field& eta);
  void force_udiv1(Field&, const Field&, const Field&);
  void force_udiv1(Field_G& force, const Field_F_1spinor& zeta,
                   const Field_F_1spinor& eta);

 private:
  void init();
  void tidyup();
};
#endif
