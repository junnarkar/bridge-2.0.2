/*!
        @file    bridge_init_factory_alt_QXS.cpp
        @brief
        @author  Hideo Matsufuru  (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$
        @version $LastChangedRevision: 2492 $
*/

#include "bridge_init_factory_alt_QXS.h"

#ifdef USE_FACTORY

// alt-code
#include "lib/Fopr/afopr.h"
//#include "lib_alt_QXS/Force/Fermion/aforce_F.h"
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt/Solver/asolver.h"
#include "lib/Eigen/aeigensolver.h"

// adding to corelib
//#include "lib_alt_QXS/Action/action.h"

/*
#include "Force/Gauge/force_G.h"
#include "Smear/forceSmear.h"
#include "Smear/smear.h"
#include "Smear/projection.h"
#include "Measurements/Gauge/gaugeFixing.h"
#include "Measurements/Gauge/staple.h"
#include "Measurements/Fermion/source.h"
#include "Tools/randomNumbers.h"
#include "Tools/gammaMatrixSet.h"
*/

void bridge_report_factory_alt_QXS();

#ifdef USE_FACTORY_AUTOREGISTER
#else

bool bridge_init_factory_alt_QXS()
{
  bool result = true;
  vout.general("hoge:  in bridge_init_factory_alt_QXS()\n");

  result &= AFopr<AField<double, QXS> >::init_factory();
  result &= ASolver<AField<double, QXS> >::init_factory();
  result &= AEigensolver<AField<double, QXS>, AFopr<AField<double, QXS> > >::init_factory();
  // result &= AForce_F<AField<double> >::init_factory();

  result &= AFopr<AField<float, QXS> >::init_factory();
  result &= ASolver<AField<float, QXS> >::init_factory();

  //  bridge_report_factory_alt_QXS();

  /*
  result &= Force_G::init_factory();
  result &= ForceSmear::init_factory();
  result &= Solver::init_factory();
  result &= Action::init_factory();
  result &= Projection::init_factory();
  result &= Smear::init_factory();
  result &= GaugeFixing::init_factory();
  result &= Source::init_factory();
  result &= Staple::init_factory();
  result &= RandomNumbers::init_factory();
  result &= GammaMatrixSet::init_factory();
#ifdef USE_FFTWLIB
  result &= FFT::init_factory();
#endif
  */

  return result;
}


#endif /* USE_FACTORY_AUTOREGISTER */

bool bridge_fin_factory_alt_QXS()
{
  bool result = true;
//  vout.general("hoge:  in bridge_fin_factory_alt_QXS()\n");
//  result &= AFopr<AField<double, QXS> >::Factory_params::clear();
//  result &= AFopr<AField<float, QXS> >::Factory_params::clear();

  return result;
}


//#ifdef DEBUG
void bridge_report_factory_alt_QXS()
{
  vout.general("------------------------------------------------\n");
  vout.general("Factory entries\n");
  vout.general("------------------------------------------------\n");

  typedef AFopr<AField<double, QXS> >   AFOPR_d;
  typedef AFopr<AField<float, QXS> >    AFOPR_f;

  AFOPR_d::Factory_params::print("AFopr<double, QXS>(Parameters&)");
  AFOPR_f::Factory_params::print("AFopr<float, QXS>(Parameters&)");

  /*
  Fopr::Factory_noarg::print("Fopr(void)");
  Fopr::Factory_fopr::print("Fopr(Fopr*)");
  Fopr::Factory_fopr_director::print("Fopr(Fopr*, Director*)");
  Fopr::Factory_string::print("Fopr(string)");

  Force_G::Factory::print("Force_G");

  ForceSmear::Factory::print("ForceSmear");

  Solver::Factory::print("Solver");

  Action::Factory::print("Action");

  Projection::Factory::print("Projection");
  Smear::Factory::print("Smear");

  GaugeFixing::Factory::print("GaugeFixing");

  Source::Factory::print("Source");

  Staple::Factory::print("Staple");

  RandomNumbers::Factory_int::print("RandomNumbers(int)");
  RandomNumbers::Factory_file::print("RandomNumbers(string)");

  GammaMatrixSet::Factory::print("GammaMatrixSet");
  */
  vout.general("------------------------------------------------\n");
}


//#endif

#endif /* USE_FACTORY */
