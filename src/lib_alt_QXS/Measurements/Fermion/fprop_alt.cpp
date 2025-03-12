/*!
        @file    fprop_alt.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$
        @version $LastChangedRevision: 2492 $
*/

#include "lib_alt/Measurements/Fermion/fprop_alt.h"

#include "lib_alt_QXS/Field/afield.h"


// explicit instanciation for AField<double>.

template class Fprop_alt<AField<double, QXS> >;

template class Fprop_alt<AField<float, QXS> >;

//============================================================END=====
