/*!
        @file    fprop_Wilson_Shift.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef FPROP_WILSON_SHIFT_INCLUDED
#define FPROP_WILSON_SHIFT_INCLUDED

#include "Fopr/fopr_Wilson.h"
#include "Field/index_lex.h"
#include "Solver/shiftsolver_CG.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Get shifted quark propagators.

/*!
    This class is used to determine the shifted quark propagator
    with multishift solver.
    Present implementation is temporary to test the shiftsolver.
                                      [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.              [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                      [21 Mar 2015 Y.Namekawa]
 */

class Fprop_Wilson_Shift
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  Fopr_Wilson *m_fopr;
  Index_lex *m_index_lex;

  int m_Niter;
  double m_Stop_cond;

  int m_Nshift;
  std::vector<double> m_sigma;

 public:
  Fprop_Wilson_Shift(Fopr_Wilson *fopr, Index_lex *index)
    : m_vl(CommonParameters::Vlevel()),
    m_fopr(fopr), m_index_lex(index) {}

  Fprop_Wilson_Shift(Fopr_Wilson *fopr, const Parameters& params)
    : m_vl(CommonParameters::Vlevel()),
    m_fopr(fopr), m_index_lex(nullptr)
  {
    set_parameters(params);
  }

 private:
  // non-copyable
  Fprop_Wilson_Shift(const Fprop_Wilson_Shift&);
  Fprop_Wilson_Shift& operator=(const Fprop_Wilson_Shift&);

 public:
  void set_parameters(const Parameters& params);
  void set_parameters(const int Nshift, const std::vector<double> sigma,
                      const int Niter, const double Stop_cond);

  void get_parameters(Parameters& params) const;

  double invert_D(std::vector<Field_F> *, const Field_F&);
};
#endif
