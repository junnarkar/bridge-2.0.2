/*
        @file    run_testmanager.cpp

        @brief

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2013-04-12 13:28:01 #$

        @version $LastChangedRevision: 1928 $
*/

#include "testlist.h"
#include "testManager.h"

//====================================================================
void testmanager_usage(const char *prog)
{
  vout.general("\n");
  vout.general("usage: %s\n", prog);
  vout.general("  (without args) -- run tests interactively.\n");
  vout.general("  test_names     -- run specified tests.\n");
  vout.general("  -a | --all     -- run all tests.\n");
  vout.general("  -l | --list    -- list registered tests.\n");
  vout.general("  -h | --help    -- show this message.\n");
  vout.general("\n");
}


//====================================================================
int run_testmanager(int argc, char **argv)
{
  TestManager& testmanager = TestManager::Instance();  // find singleton

#ifdef USE_TESTMANAGER_AUTOREGISTER
  //- tests are registered automatically at start up.
#else
  //- register tests manually:
  //
  // testmanager.registerTest(
  //   "Gauge.Plaquette",
  //   Test_Gauge::plaquette
  // );
#endif

  if (argc > 1) {
    const std::string arg1 = argv[1];

    if ((arg1 == "-a") || (arg1 == "--all")) {         // run all tests
      testmanager.batch_recursive();
    } else if ((arg1 == "-l") || (arg1 == "--list")) { // list registered tests
      testmanager.list();
    } else if ((arg1 == "-h") || (arg1 == "--help")) { // show usage and exit
      testmanager_usage(argv[0]);
      return EXIT_SUCCESS;
    } else {
      // argv[] = { "./bridge.elf", test1, (test2, ...) }
      testmanager.batch_recursive(argc - 1, argv + 1);
    }
  } else {
    testmanager.interactive();
  }

  return EXIT_SUCCESS;
}


//====================================================================
//====================================================================
