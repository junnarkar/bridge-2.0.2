/*!
        @file    bridge_setup.h

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#ifndef BRIDGE_SETUP_INCLUDED
#define BRIDGE_SETUP_INCLUDED

#include <vector>
#include <string>
#include "Parameters/parameters.h"

#ifdef USE_GROUP_SU3
#define Nc    3
#else
#ifdef USE_GROUP_SU2
#define Nc    2
#else
#ifdef USE_GROUP_SU_N
#define Nc    3
#endif
#endif
#endif

int bridge_initialize(int *pargc, char ***pargv);

int bridge_initialize(int *pargc, char ***pargv,
                      const Parameters& params);

int bridge_initialize(int *pargc, char ***pargv,
                      const std::vector<int>& lattice_size,
                      const std::vector<int>& grid_size = std::vector<int>(),
                      const int number_of_threads = 1,
                      const int number_of_colors       = Nc,
                      const std::string& logfile       = "stdout",
                      const std::string& ildg_logfile  = "stdout",
                      const std::string& verbose_level = "General"
                      );

int bridge_finalize();

void bridge_setup(const Parameters& params);

void bridge_setup(const std::vector<int>& lattice_size,
                  const std::vector<int>& grid_size = std::vector<int>(),
                  const int number_of_threads = 1,
                  const int number_of_colors       = Nc,
                  const std::string& logfile       = "stdout",
                  const std::string& ildg_logfile  = "stdout",
                  const std::string& verbose_level = "General"
                  );

#ifdef Nc
#undef Nc
#endif

#endif /* BRIDGE_SETUP_INCLUDED */
