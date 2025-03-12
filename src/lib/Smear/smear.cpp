/*!
        @file    smear.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:01:33 #$

        @version $LastChangedRevision: 1928 $
*/

#include "smear.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
// setup factories for all subclasses

#include "smear_APE.h"
#include "smear_APE_spatial.h"
#include "smear_HYP.h"
#include "smear_APE_SF.h"
#include "smear_HYP_SF.h"

bool Smear::init_factory()
{
  bool result = true;

  result &= Smear_APE::register_factory();
  result &= Smear_APE_spatial::register_factory();
  result &= Smear_HYP::register_factory();
  result &= Smear_APE_SF::register_factory();
  result &= Smear_HYP_SF::register_factory();

  return result;
}


#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */
