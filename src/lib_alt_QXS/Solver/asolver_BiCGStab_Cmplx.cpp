/*!

        @file    asolver_BiCGStab_Cmplx.cpp

        @brief   BiCGStab Solver

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate: 2023-02-28 16:09:41 +0900 (Tue, 28 Feb 2023) $

        @version $LastChangedRevision: 2492 $
*/
//====================================================================

#include "lib_alt/Solver/asolver_BiCGStab_Cmplx.h"
#include "lib/ResourceManager/threadManager.h"
#include "lib_alt_QXS/Field/afield.h"

// function template files
#include "lib_alt_QXS/Field/afield-inc.h"

// template functions
#include "lib_alt/Solver/asolver_BiCGStab_Cmplx-tmpl.h"


//====================================================================
// explicit instanciation for AField<double>.
template<>
const std::string ASolver_BiCGStab_Cmplx<AField<double, QXS> >::class_name
  = "ASolver_BiCGStab_Cmplx<Afield<double,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_BiCGStab_Cmplx<AField<double, QXS> >::register_factory();
}
#endif

template class ASolver_BiCGStab_Cmplx<AField<double, QXS> >;

//====================================================================
// explicit instanciation for AField<float>.
template<>
const std::string ASolver_BiCGStab_Cmplx<AField<float, QXS> >::class_name
  = "ASolver_BiCGStab_Cmplx<AField<float,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_BiCGStab_Cmplx<AField_dev<float, QXS> >::register_factory();
}
#endif

template class ASolver_BiCGStab_Cmplx<AField<float, QXS> >;

//============================================================END=====
