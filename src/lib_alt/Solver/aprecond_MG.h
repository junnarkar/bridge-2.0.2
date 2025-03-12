/*!
      @file    aprecond_MG.h
      @brief   Multigrid preconditionor (SIMD version)
      @author  KANAMORI Issaku (kanamori)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2023-02-28 16:09:41 +0900 (Tue, 28 Feb 2023) $
      @version $LastChangedRevision: 2492 $
*/

#ifndef APRECOND_MG_H
#define APRECOND_MG_H

#include <cstdio>
#include <cstdlib>

#include <string>
using std::string;

#include <vector>
using std::vector;

#include  "lib/Parameters/commonParameters.h"
#include  "lib/IO/bridgeIO.h"
using Bridge::vout;
#include "lib/Tools/timer.h"
#include "lib/Fopr/afopr.h"

#include "lib_alt/Solver/asolver_BiCGStab.h"
#include "lib_alt/Solver/asolver_FBiCGStab.h"
#include "lib_alt/Solver/asolver_BiCGStab_Cmplx.h"
#include "lib_alt/Solver/MultiGrid.h"
#include "lib_alt/Solver/aprecond.h"

/*!
    explanation to be added.
 */

template<typename AFIELD, typename AFIELD2>
class APrecond_MG : public APrecond<AFIELD>
{
 public:
  typedef typename AFIELD::real_t         real_t0; // assumes double
  typedef typename AFIELD2::real_t        real_t;  // assumes float
  typedef MultiGrid<AFIELD2, AFIELD2>     MultiGrid_t;
  typedef typename MultiGrid_t::Index_t   block_index_t;

  static const Impl IMPL  = AFIELD::IMPL;
  static const Impl IMPL2 = AFIELD2::IMPL;

  static const std::string class_name;

 private:
  ASolver<AFIELD2> *m_coarse_solver;
  ASolver<AFIELD2> *m_smoother;
  const MultiGrid_t *m_multigrid;

  AFopr<AFIELD> *m_afoprD;  // double
  AFopr<AFIELD2> *m_afoprF; // float

  AFIELD2 m_coarse_v, m_coarse_w;
  AFIELD2 m_fine_w, m_fine_v, m_fine_tmp, m_fine_r;
  AFIELD m_fine_vecD, m_fine_debugD, m_fine_rD;
  //  const std::vector<AFIELD2> *m_testvectors;

  double m_accum_flop_coarse;
  double m_accum_flop_smoother;
  double m_accum_flop_other;
  double m_accum_flop_float;
  double m_accum_flop_double;
  double m_accum_flop_restrict;
  double m_accum_flop_prolong;

  int m_coarse_ncol;

  // elapsed time for each ingredient
  double m_time_f2d;
  double m_time_d2f;
  double m_time_residual;
  double m_time_restriction;
  double m_time_coarse_solver;
  double m_time_smoother;
  double m_time_prolongation;
  double m_time_mult_single_total;
  double m_time_mult_total;
  double m_time_double;
  int m_num_mult_called;        //!< number of mult call
  int m_num_mult_single_called; //!< number of mult call

  //  unique_ptr<Timer> m_timer_f2d;  // float to double
  //  unique_ptr<Timer> m_timer_d2f;  // double to flaot
  //  unique_ptr<Timer> m_timer_residual;
  //  unique_ptr<Timer> m_timer_restrcition;
  //  unique_ptr<Timer> m_timer_coarse_solver;
  //  unique_ptr<Timer> m_timer_smoother;
  //  unique_ptr<Timer> m_timer_prolongation;
  //  unique_ptr<Timer> m_timer_mult_single_total;

 public:
  APrecond_MG() { init(); }

  void mult(AFIELD&, const AFIELD&);
  void mult_as_setup(AFIELD2&, const AFIELD2&);
  void reset_flop_count();
  void set_coarse_solver(ASolver<AFIELD2> *solver);
  void set_smoother(ASolver<AFIELD2> *solver);

  void set_multigrid(const MultiGrid_t *multigrid)
  {
    m_multigrid = multigrid;
  }

  void set_fopr(AFopr<AFIELD> *foprD, AFopr<AFIELD2> *foprF) { m_afoprD = foprD; m_afoprF = foprF; }
  void report_timer();

  double flop_count() { return m_accum_flop_float + m_accum_flop_double; }
  double flop_count_coarse() const { return m_accum_flop_coarse; }
  double flop_count_smoother() const { return m_accum_flop_smoother; }
  double flop_count_other() const { return m_accum_flop_other; }
  double flop_count_double() const { return m_accum_flop_double; }

 private:
  void init();
  void tidyup();
  void update_flop_count();
  void mult_single(AFIELD2&, const AFIELD2&);
};

#endif // include guard
//============================================================END=====
