/*!
        @file    afopr_SAP_MINRES_float.cpp
        @brief   MINRES solver inside a SAP solver (QXS version)
        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate::  $
        @version $LastChangedRevision: 2492 $
 */

//====================================================================
#include "lib_alt/Solver/asolver_SAP_MINRES.h"
#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_QXS/inline/define_vlen.h"

#define  VLEN     VLENS
#define  VLENX    VLENXS
#define  VLENY    VLENYS

typedef float real_t;

#include "lib_alt_QXS/inline/vsimd_float-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_float-inc.h"
#include "lib_alt_QXS/Field/aindex_block_lex.h"
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/afield_dd-inc.h"
#include "lib_alt_QXS/Fopr/afopr_Clover_dd.h"

//#define DEBUG_MINRES

template<typename AFIELD>
using AFopr_dd_t = AFopr_Clover_dd<AFIELD>;


#include "lib_alt/Solver/asolver_SAP_MINRES-tmpl.h"

//====================================================================
// explicit instanciation for AField<float>.
template<>
const std::string ASolver_SAP_MINRES<AField<float, QXS> >::class_name
  = "ASolver_SAP_MINRES<AField<float,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_SAP_MINRES<AField<float, QXS> >::register_factory();
}
#endif

template class ASolver_SAP_MINRES<AField<float, QXS> >;

//============================================================END=====
