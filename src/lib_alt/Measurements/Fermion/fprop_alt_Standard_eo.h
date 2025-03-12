/*!
        @file    fprop_alt_Standard_eo.h
        @brief
        @author  Satoru Ueda  (sueda)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$
        @version $LastChangedRevision: 2492 $
*/

#ifndef FPROP_ALT_STANDARD_EO_INCLUDED
#define FPROP_ALT_STANDARD_EO_INCLUDED

#include "lib_alt/Measurements/Fermion/fprop_alt.h"

#include "lib/Tools/timer.h"

#include "lib/Fopr/afopr.h"
#include "lib/Smear/director_Smear.h"

#include "lib_alt/Solver/asolver.h"

//! Get quark propagator for Fopr with lexical site index: alternative version.

/*!
    This is temporary implementation.
                                        [30 May 2017 H.Matsufuru]
 */

template<typename AFIELD>
class Fprop_alt_Standard_eo : public Fprop_alt<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;
  using Fprop_alt<AFIELD>::m_vl;
  using Fprop_alt<AFIELD>::m_mode;

 private:
  AFopr<AFIELD> *m_kernel;
  AFopr<AFIELD> *m_fopr;
  ASolver<AFIELD> *m_solver;

  Timer m_timer;
  double m_flop_count;
  double m_elapsed_time;

  Director_Smear *m_dr_smear;

 public:
  Fprop_alt_Standard_eo(const Parameters& params_fopr,
                        const Parameters& params_solver)
    : Fprop_alt<AFIELD>()
  { init(params_fopr, params_solver); }

  Fprop_alt_Standard_eo(const Parameters& params_fopr,
                        const Parameters& params_solver,
                        Director_Smear *dr_smear)
    : Fprop_alt<AFIELD>()
  { init(params_fopr, params_solver, dr_smear); }

  ~Fprop_alt_Standard_eo()
  { tidyup(); }

  void set_config(Field *);

  void invert(Field&, const Field&, int&, double&);

  void invert_D(Field&, const Field&, int&, double&);

  void invert_DdagD(Field&, const Field&, int&, double&);

  // inverter with AFIELD
  void invert(AFIELD&, const AFIELD&, int&, double&);

  double flop_count();

  void reset_performance();

  void report_performance();

  void mult_performance(const std::string mode, const int Nrepeat);

 private:

  void init(const Parameters& params_fopr,
            const Parameters& params_solver);

  void init(const Parameters& params_fopr,
            const Parameters& params_solver,
            Director_Smear *dr_smear);

  void tidyup();

  void invert_D(AFIELD&, const AFIELD&, int&, double&);

  void invert_DdagD(AFIELD&, const AFIELD&, int&, double&);

  void invert_De(AFIELD&, AFIELD&, AFIELD&, AFIELD&, int&, double&);
};
#endif
