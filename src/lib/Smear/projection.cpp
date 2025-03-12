/*!
        @file    projection.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:01:33 #$

        @version $LastChangedRevision: 1928 $
*/

#include "projection.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
// setup factories for all subclasses

#include "projection_Maximum_SU_N.h"
#include "projection_Stout_SU3.h"

bool Projection::init_factory()
{
  bool result = true;

  result &= Projection_Maximum_SU_N::register_factory();
  result &= Projection_Stout_SU3::register_factory();

  return result;
}


#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */
