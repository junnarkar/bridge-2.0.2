/*!
      @file    shiftAField_eo_double.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-02-28 16:09:41 #$
      @version $LastChangedRevision: 2492 $
*/


#include "lib_alt_QXS/Field/shiftAField_eo.h"

#include <string>

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"

#define  VLEN     VLEND
#define  VLENX    VLENXD
#define  VLENY    VLENYD

typedef double real_t;

#include "lib_alt_QXS/Field/aindex_eo.h"
#include "lib_alt_QXS/Field/afield.h"

// inline function files
#include "lib_alt_QXS/inline/vsimd_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_double-inc.h"
#include "lib_alt_QXS/inline/afield_th-inc.h"

// function template files
#include "lib_alt_QXS/Field/afield-inc.h"

#include "lib_alt_QXS/Field/shiftAField_eo-tmpl.h"

template<>
const std::string ShiftAField_eo<AField<double, QXS> >::class_name
  = "ShiftAField_eo<AField<double,QXS> >";

// explicit instanciation.
template class ShiftAField_eo<AField<double, QXS> >;

//============================================================END=====
