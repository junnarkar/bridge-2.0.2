/*!
        @file    fopr_CloverTerm_General.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#ifndef FOPR_CLOVERTERM_GENERAL_INCLUDED
#define FOPR_CLOVERTERM_GENERAL_INCLUDED

// implementations

// original unoptimised version
#ifdef USE_ORG
#include "Org/fopr_CloverTerm_General_impl.h"
#endif

// improved version
#ifdef USE_IMP
#include "Imp/fopr_CloverTerm_General_impl.h"
#endif

#if defined(USE_IMP)
typedef Imp::Fopr_CloverTerm_General   Fopr_CloverTerm_General;
#else
typedef Org::Fopr_CloverTerm_General   Fopr_CloverTerm_General;
#endif

#endif
