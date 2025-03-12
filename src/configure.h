/*!
        @file    configure.h

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#ifndef CONFIGURE_INCLUDED
#define CONFIGURE_INCLUDED

#define BRIDGE_VERSION    "2.0.0"

// some #define options are set as compiler options from Makefile

//#define USE_MPI
#define ENABLE_MULTI_INSTANCE

// smart pointer support
#include <memory>
using std::unique_ptr;

#define DEPRECATED    [[deprecated]]
// #define DEPRECATED __attribute__((deprecated))
// #define DEPRECATED

#endif /* CONFIGURE_INCLUDED */
