// Time-stamp: <2018-10-22 17:47:39 kanamori>

/*!

        @file    test_MG_solver.h

        @brief   test for MultiGrid solver (alt_Simd version)

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: #$

        @version $LastChangedRevision: 2492 $
*/


#ifndef TEST_MG_SOLVER_H
#define TEST_MG_SOLVER_H

#include <vector>
#include <string>

#include "lib/Parameters/parameters.h"
#include "lib/Field/field_G.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//class Fopr;
//class Source;

class Test_MG_Solver
{
 public:
  static const std::string class_name;

 private:
  Parameters params_all;
  Bridge::VerboseLevel m_vl;
  unique_ptr<Field_G> U;

 public:

  //! constructor
  Test_MG_Solver()
    : m_vl(CommonParameters::Vlevel())
  { init(); }

  //! destructor
  ~Test_MG_Solver() {}

  int test();
  int test_jps();

 private:

  //! initial setup
  void init();
};
#endif //
