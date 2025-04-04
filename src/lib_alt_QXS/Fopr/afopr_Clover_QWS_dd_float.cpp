/*!
      @file    afopr_Clover_QWS_dd_float.cpp
      @brief
      @author  $Author: Issaku Kanamori$
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2492 $
*/

#include "lib_alt_QXS/Fopr/afopr_Clover_QWS_dd.h"
#include "lib/ResourceManager/threadManager.h"
#include "lib_alt_QXS/inline/define_params.h"

#define  VLEN     VLENS
#define  VLENX    VLENXS
#define  VLENY    VLENYS

//#if defined(USE_QWSLIB) && VLEND != VLENXD
//#warning VLEN is not 1-dim, not using qws
//#undef USE_QWSLIB
//#endif

typedef float real_t;

#include "lib_alt_QXS/inline/vsimd_float-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_float-inc.h"
#include "lib_alt_QXS/inline/vsimd_Wilson_SU3_float-inc.h"

#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/aindex_lex_QWS_dd.h"
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"

#include "lib_alt_QXS/Fopr/mult_common_th-inc.h"

//#ifdef USE_QXSLIB
#include "lib_alt_QXS/BridgeQXS/bridgeQXS_Wilson.h"
#include "lib_alt_QXS/BridgeQXS/bridgeQXS_Clover.h"
//#endif

#include "lib_alt_QXS/Fopr/mult_Wilson_parts_qxs_org-inc.h"

// template definition
#include "lib_alt_QXS/Fopr/afopr_Clover_QWS_dd-tmpl.h"

template<>
const std::string AFopr_Clover_QWS_dd<AField<float, QXS> >::class_name
  = "AFopr_Clover_QWS_dd<AField<float,QXS> >";


// explicit instanciation.
template class AFopr_Clover_QWS_dd<AField<float, QXS> >;

//============================================================END=====
