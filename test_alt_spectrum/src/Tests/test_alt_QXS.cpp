/*!
      @File    test_alt_QXS.cpp
      @brief
      @author  Hideo Matsufuru
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-02-28 16:09:41 #$
      @version $LastChangedRevision: 2492 $
*/

#ifdef USE_ALT_QXS

#include "lib_alt_QXS/bridge_alt_qxs.h"
#define IMPL    QXS


#include "spectrum_Wilson_alt-tmpl.h"
#include "spectrum_Domainwall_alt-tmpl.h"
#include "spectrum_Staggered_alt-tmpl.h"

namespace Test_alt_QXS {
  int test_all();
  int test_Spectrum_Wilson();
  int test_Spectrum_Clover();
  int test_Spectrum_Staggered();
  int test_Spectrum_Domainwall();
  int test_Spectrum_Domainwall_5din();
}

//====================================================================
int Test_alt_QXS::test_all()
{
  int result = 0;
  result += test_Spectrum_Wilson();
  result += test_Spectrum_Clover();
  result += test_Spectrum_Staggered();
  result += test_Spectrum_Domainwall_5din();
  result += test_Spectrum_Domainwall();

  vout.general("Test result for QXS = %d\n", result);
  return result;
}


//====================================================================
int Test_alt_QXS::test_Spectrum_Wilson()
{
  Spectrum_Wilson_alt<IMPL> test_wilson;
  string file_params = "test_alt_Spectrum_Wilson_Hadron2ptFunction.yaml";

  int result = 0;
  result += test_wilson.hadron_2ptFunction(file_params, "double");
  test_wilson.hadron_2ptFunction(file_params, "float");
  result += test_wilson.hadron_2ptFunction(file_params, "double_eo");
  test_wilson.hadron_2ptFunction(file_params, "float_eo");

  return result;
}


//====================================================================
int Test_alt_QXS::test_Spectrum_Clover()
{
  Spectrum_Wilson_alt<IMPL> test_wilson;
  string file_params = "test_alt_Spectrum_Clover_Hadron2ptFunction.yaml";

  int result = 0;
  result += test_wilson.hadron_2ptFunction(file_params, "double");
  test_wilson.hadron_2ptFunction(file_params, "float");
  result += test_wilson.hadron_2ptFunction(file_params, "double_eo");
  test_wilson.hadron_2ptFunction(file_params, "float_eo");

  return result;
}


//====================================================================
int Test_alt_QXS::test_Spectrum_Staggered()
{
  Spectrum_Staggered_alt<IMPL> test_staggered;

  std::string run_mode = "test";

  int result = 0;

  string file_cube =
    "test_alt_Spectrum_Staggered_Hadron2ptFunction_Wall_Cube.yaml";
  result += test_staggered.hadron_2ptFunction_Cube(file_cube,
                                                   "double", run_mode);

  test_staggered.hadron_2ptFunction_Cube(file_cube,
                                         "float", run_mode);

  result += test_staggered.hadron_2ptFunction_Cube(file_cube,
                                                   "double_eo", run_mode);

  test_staggered.hadron_2ptFunction_Cube(file_cube,
                                         "float_eo", run_mode);

  string file_evenodd =
    "test_alt_Spectrum_Staggered_Hadron2ptFunction_Wall_Evenodd.yaml";
  result += test_staggered.hadron_2ptFunction_Evenodd(file_evenodd,
                                                      "double", run_mode);

  test_staggered.hadron_2ptFunction_Evenodd(file_evenodd,
                                            "float", run_mode);

  result += test_staggered.hadron_2ptFunction_Evenodd(file_evenodd,
                                                      "double_eo", run_mode);

  test_staggered.hadron_2ptFunction_Evenodd(file_evenodd,
                                            "float_eo", run_mode);

  return result;
}


//====================================================================
int Test_alt_QXS::test_Spectrum_Domainwall_5din()
{
  Spectrum_Domainwall_alt<IMPL> test_dw;
  string test_file = "test_alt_Spectrum_Domainwall_5din_Hadron2ptFunction.yaml";

  int result = 0;
  result += test_dw.hadron_2ptFunction(test_file, "double");
  result += test_dw.hadron_2ptFunction(test_file, "float");
  result += test_dw.hadron_2ptFunction(test_file, "double_eo");
  result += test_dw.hadron_2ptFunction(test_file, "float_eo");

  return result;
}


//====================================================================
int Test_alt_QXS::test_Spectrum_Domainwall()
{
  Spectrum_Domainwall_alt<IMPL> test_dw;
  string test_file = "test_alt_Spectrum_Domainwall_Hadron2ptFunction.yaml";

  int result = 0;
  result += test_dw.hadron_2ptFunction(test_file, "double");
  result += test_dw.hadron_2ptFunction(test_file, "float");
  result += test_dw.hadron_2ptFunction(test_file, "double_eo");
  result += test_dw.hadron_2ptFunction(test_file, "float_eo");

  return result;
}


#endif
//============================================================END=====
