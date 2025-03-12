/*!
      @file    afopr_Clover_eo_float.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2492 $
*/

#include "lib_alt_QXS/Fopr/afopr_Clover_eo.h"

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"

#define  VLEN     VLENS
#define  VLENX    VLENXS
#define  VLENY    VLENYS

typedef float real_t;

#define CHIRAL_ROTATION    // chiral rotation in clover term

#include "lib_alt_QXS/inline/vsimd_float-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_float-inc.h"
#include "lib_alt_QXS/inline/vsimd_Wilson_SU3_float-inc.h"

#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/afield_Gauge-inc.h"

#include "lib_alt_QXS/Fopr/mult_common_th-inc.h"
#include "lib_alt_QXS/Fopr/mult_Wilson_qxs_parts-inc.h"

#include "lib_alt_QXS/BridgeQXS/bridgeQXS_Wilson.h"
#include "lib_alt_QXS/BridgeQXS/bridgeQXS_Clover.h"


// template definition
#include "lib_alt_QXS/Fopr/afopr_Clover_eo-tmpl.h"

template<>
const std::string AFopr_Clover_eo<AField<float, QXS> >::class_name
  = "AFopr_Clover_eo<AField<float,QXS> >";


// explicit instanciation.
template class AFopr_Clover_eo<AField<float, QXS> >;

//============================================================END=====
