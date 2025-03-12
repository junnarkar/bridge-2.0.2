/*!
      @file    aindex_eo_float.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-02-28 16:09:41 #$
      @version $LastChangedRevision: 2492 $
*/

#include <assert.h>

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/Field/aindex_eo.h"

#define  VLEN     VLENS
#define  VLENX    VLENXS
#define  VLENY    VLENYS

#include "lib_alt_QXS/Field/afield.h"

#include "lib_alt_QXS/Field/aindex_eo-tmpl.h"

// explicit instanciation.
template class AIndex_eo<float, QXS>;
//============================================================END=====
