/*!
        @file    eigensolver_IRArnoldi.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#include "Eigen/aeigensolver_IRArnoldi.h"
#include "Eigen/aeigensolver_IRArnoldi-tmpl.h"
#include "Eigen/eigensolver_IRArnoldi.h"

//#include "Field/field.h"
//#include "Fopr/fopr.h"

#define COMPLEX    dcomplex

// explicit instanciation for AField<double>.
template<>
const std::string AEigensolver_IRArnoldi<Field, Fopr>::class_name
  = "Eigensolver_IRArnoldi";

template class AEigensolver_IRArnoldi<Field, Fopr>;


//============================================================END=====
