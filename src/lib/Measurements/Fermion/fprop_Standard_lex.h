/*!
        @file    fprop_Standard_lex.h

        @brief

        @author  Satoru Ueda  (sueda)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#ifndef FPROP_STANDARD_LEX_INCLUDED
#define FPROP_STANDARD_LEX_INCLUDED

#include "fprop.h"
#include "Solver/solver.h"

//! Get quark propagator for Fopr with lexical site index.

/*!
    This is temporary implementation.
                                        [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    Redundant argument of fopr is deleted.
                                        [03 Mar 2013 Y.Namekawa]
    Introduce unique_ptr to avoid memory leaks.
                                        [21 Mar 2015 Y.Namekawa]
    Add flop_count.                     [ 8 Aug 2016 Y.Namekawa]
 */

class Fprop_Standard_lex : public Fprop
{
 public:
  static const std::string class_name;

 private:
  Solver *m_solver;

 public:
  Fprop_Standard_lex(Solver *solver)
    : Fprop(), m_solver(solver) {}

  ~Fprop_Standard_lex() {}

  void set_config(Field *);

  void invert_D(Field&, const Field&, int&, double&);
  void invert_DdagD(Field&, const Field&, int&, double&);

  double flop_count();

  void mult_performance(const std::string mode, const int Nrepeat);
};
#endif
