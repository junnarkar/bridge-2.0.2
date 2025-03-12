/*!
        @file    fopr_Wilson_TwistedMass.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#include "Fopr/fopr_Wilson_TwistedMass.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Fopr_Wilson_TwistedMass::register_factory();
}
#endif

#include "Fopr/afopr_Wilson_TwistedMass-tmpl.h"

template<>
const std::string Fopr_Wilson_TwistedMass::class_name
  = "Fopr_Wilson_TwistedMass";

template class AFopr_Wilson_TwistedMass<Field>;

//============================================================END=====
