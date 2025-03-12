/*!
        @file    fopr_Clover_Chemical.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#include "Fopr/fopr_Clover_Chemical.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Fopr_Clover_Chemical::register_factory();
}
#endif

#include "Fopr/afopr_Clover_Chemical-tmpl.h"

template<>
const std::string Fopr_Clover_Chemical::class_name
  = "Fopr_Clover_Chemical";

template class AFopr_Clover_Chemical<Field>;

//============================================================END=====
