/*
        @file    run_test_MG.cpp

        @brief

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate: 2023-02-28 16:09:41 +0900 (Tue, 28 Feb 2023) $

        @version $LastChangedRevision: 2492 $
*/

#include <stdlib.h>
#include "testlist_MG.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//====================================================================
int run_test_MG()
{
  //test_coarse::test();
  //  test_restrict::test();
  //test_prolong::test();
  //test_afopr_coarse::test_parallel();
  //test_sap_mult::test();
  test_MG_solver::test();

  return EXIT_SUCCESS;
}
