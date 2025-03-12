/*!
        @file    asolver_SAP.cpp
        @brief   SAP solver (QXS version)
        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate: 2023-02-28 16:09:41 +0900 (Tue, 28 Feb 2023) $
        @version $LastChangedRevision: 2492 $
*/
//====================================================================
#include "lib_alt/Solver/asolver_SAP.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"

//#define DEBUG

#include "lib_alt/Solver/asolver_SAP-tmpl.h"


//====================================================================
// explicit instanciation
template<>
const std::string ASolver_SAP<AField<double, QXS> >::class_name
  = "ASolver_SAP<AField<double,QXS> >";

template class ASolver_SAP<AField<double, QXS> >;

//====================================================================
// explicit instanciation
template<>
const std::string ASolver_SAP<AField<float, QXS> >::class_name
  = "ASolver_SAP<AField<float,QXS> >";

template class ASolver_SAP<AField<float, QXS> >;

//============================================================END=====
