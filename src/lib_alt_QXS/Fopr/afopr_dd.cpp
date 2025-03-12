/*!
        @file    afopr.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$
        @version $LastChangedRevision: 2492 $
*/

#include "lib_alt/Fopr/afopr_dd.h"
#include "lib_alt_QXS/Field/afield.h"

// explicit instanciation.
template class AFopr_dd<AField<double, QXS> >;
template class AFopr_dd<AField<float, QXS> >;

//============================================================END=====
