/*!
        @file    aforce_F_Smeared.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$
        @version $LastChangedRevision: 2492 $
*/

#ifndef AFORCE_F_SMEARED_INCLUDED
#define AFORCE_F_SMEARED_INCLUDED

#include "lib/Force/Fermion/aforce_F.h"
#include "lib/Smear/forceSmear.h"
#include "lib/Smear/director_Smear.h"
#include "lib/Fopr/afopr.h"


//! Force calculation for smeared fermion operators.

/*!
    This class determines the force of smeared fermion operator
    using smearing director (MultiSmear instance) and base
    fermion force instance.
                                      [28 Dec 2011 H.Matsufuru]
    Modified: set_mode() is added to incorporate non-hermitian H
                                      [21 Jan 2012 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                      [21 Mar 2015 Y.Namekawa]
*/

template<typename AFIELD>
class AForce_F_Smeared : public AForce_F<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  using AForce_F<AFIELD>::m_vl;
  using AForce_F<AFIELD>::m_U;
  using AForce_F<AFIELD>::m_Ucp;
  static const std::string class_name;

 private:
  AForce_F<AFIELD> *m_force;
  Director_Smear *m_director_smear;
  // Note that this template version always smear the force
  // voa Director_Smear.

 public:
  AForce_F_Smeared(AForce_F<AFIELD> *force,
                   Director_Smear *director_smear)
    : AForce_F<AFIELD>(), m_force(force), m_director_smear(director_smear)
  { init(); }

  void set_parameters(const Parameters&);

  void set_config(Field *U);

  void set_mode(const std::string& mode)
  { m_force->set_mode(mode); }

  void force_udiv(AFIELD& force, const AFIELD& eta);

  void force_udiv1(AFIELD& force, const AFIELD& zeta, const AFIELD& eta);

 private:
  //! initial setup.
  void init();

  void mult_jacobian(Field_G& force);

#ifdef USE_FACTORY
 private:
  static AForce_F<AFIELD> *create_object_with_force_director(
    AForce_F<AFIELD> *fopr, Director *director)
  { return new AForce_F_Smeared(fopr, (Director_Smear *)director); }

 public:
  static bool register_factory()
  {
    bool init1 = AForce_F<AFIELD>::Factory_force_director::Register(
      "Smeared", create_object_with_force_director);
    return init1;
  }
#endif
};
#endif
