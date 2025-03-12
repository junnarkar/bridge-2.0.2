/*!
        @file    fopr_CloverTerm.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#ifndef FOPR_CLOVERTERM_INCLUDED
#define FOPR_CLOVERTERM_INCLUDED

//! Clover term operator.

/*!
    This class implements the clover term for the clover (improved
    Wilson) fermion operator.
    This part was separated from the Fopr_Clover class.
    The field strength is calculate when the function
    set_config() is called.
    The `mode' for setting fermion operator mode is now only
    defined to the case 'D'.
                [30 Sep 2012 H.Matsufuru,
                 original clover operator: 24 Dec 2011 H.M.]
    (Coding history will be recovered from trac.)
    YAML is implemented.            [14 Nov 2012 Y.Namekawa]
    Selector is implemented.        [03 Mar 2013 Y.Namekawa]
    (Selectors are replaced with factories by Aoyama-san)
    unique_ptr is introduced to avoid memory leaks
                                    [21 Mar 2015 Y.Namekawa]
    A mode 'F' is added.
    Note: mult with mode 'D' or 'F' multiplies

      csw kappa sigma_{mu nu} F_{mu nu}.

    (this is different from that in fopr_CloverTerm_eo,
     which multiplies 1 - csw kappa sigma_{mu nu} F_{mu nu}. )
                                    [12 Jan 2019 I.Kanamori]
 */

// implementations

#ifdef USE_ORG
#include "Org/fopr_CloverTerm_impl.h"
#endif

#ifdef USE_IMP
#include "Imp/fopr_CloverTerm_impl.h"
#endif

#if defined(USE_IMP)
typedef Imp::Fopr_CloverTerm   Fopr_CloverTerm;
#else
typedef Org::Fopr_CloverTerm   Fopr_CloverTerm;
#endif

#endif
