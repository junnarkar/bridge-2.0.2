/*!
      @file    aindex_eo_alt.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-02-28 16:09:41 #$
      @version $LastChangedRevision: 2492 $
*/

#include <assert.h>

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/Field/aindex_eo.h"

#define  VLEN     VLEND
#define  VLENX    VLENXD
#define  VLENY    VLENYD

#include "lib_alt_QXS/Field/afield.h"

#include "lib_alt_QXS/Field/aindex_eo-tmpl.h"

// explicit instanciation.
template class AIndex_eo<double, QXS>;
//============================================================END=====
