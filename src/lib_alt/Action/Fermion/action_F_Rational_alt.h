/*!
      @file    action_F_Rational_alt.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-02-28 16:09:41 #$
      @version $LastChangedRevision: 2492 $
*/

#ifndef ACTION_F_RATIONAL_ALT_INCLUDED
#define ACTION_F_RATIONAL_ALT_INCLUDED

#include "lib/Action/action.h"

#include "lib/Fopr/afopr_Rational.h"
#include "lib/Force/Fermion/aforce_F_Rational.h"

#include "lib/IO/bridgeIO.h"
using Bridge::vout;

// include files in alt-code
#include "Field/afield.h"
#include "Measurements/Fermion/afprop.h"

//! action class for RHMC, with externally constructed AFopr_Rational.

/*!
    Action class for RHMC that is an alternative to Action_F_Rational
    in the core library.
                                            [05 Feb 2019 H.Matsufuru]
 */

template<typename AFIELD>
class Action_F_Rational_alt : public Action
{
 public:
  //typedef AField<double> AFIELD;
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  std::string m_label;    // label of action

  AFopr<AFIELD> *m_fopr_langev;
  AFopr<AFIELD> *m_fopr_H;
  AForce_F<AFIELD> *m_fopr_force_MD;

  Field *m_U;

  AFIELD m_psf;

 public:
  //! constructor.
  Action_F_Rational_alt(AFopr<AFIELD> *fopr_langev,
                        AFopr<AFIELD> *fopr_H,
                        AForce_F<AFIELD> *fopr_force_MD)
    : Action(),
    m_fopr_langev(fopr_langev), m_fopr_H(fopr_H),
    m_fopr_force_MD(fopr_force_MD)
  { init(); }

  //! destructor.
  ~Action_F_Rational_alt() { tidyup(); }

  void set_parameters(const Parameters& params);

  //! set the label of action.
  void set_label(const std::string label)
  {
    m_label = label;
    vout.detailed(m_vl, "  label: %s\n", m_label.c_str());
  }

  //! returns the label of action.
  std::string get_label()
  { return m_label; }

  //! setting gauge configuration.
  void set_config(Field *U);

  //! Langevin step called at the beginning of HMC.
  double langevin(RandomNumbers *);

  //! calculation of Hamiltonian.
  double calcH();

  //! returns the force for updating conjugate momentum.
  //const Field force();
  void force(Field&);

 private:
  void init();
  void tidyup();
};
#endif
