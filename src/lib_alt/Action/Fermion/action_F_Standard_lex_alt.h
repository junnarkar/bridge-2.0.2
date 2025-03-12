/*!
      @file    action_F_Standard_lex_alt.h
      @brief
      @author  Hideo Matsufuru
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-02-28 16:09:41 #$
      @version $LastChangedRevision: 2492 $
*/

#ifndef ACTION_F_STANDARD_LEX_ALT_INCLUDED
#define ACTION_F_STANDARD_LEX_ALT_INCLUDED

#include "lib/Action/action.h"

#include "lib/Force/Fermion/force_F.h"
#include "lib/Measurements/Fermion/fprop.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

// include files in alt-code
//#include "Field/afield.h"
#include "lib_alt/Measurements/Fermion/afprop.h"

//! Standard fermion action for HMC.

/*!
    Standard fermion action with alternative implementation.
    [04 Oct 2018 H.Matsufuru]
 */

template<typename AFIELD>
class Action_F_Standard_lex_alt : public Action
{
 public:
  //  typedef AField<double> AFIELD;
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  AFopr<AFIELD> *m_fopr;
  AForce_F<AFIELD> *m_fopr_force;
  AFIELD m_psf;
  std::string m_label;

  AFprop<AFIELD> *m_fprop_MD;
  AFprop<AFIELD> *m_fprop_H;

  Field *m_U;

 public:
  Action_F_Standard_lex_alt(
    AFopr<AFIELD> *fopr, AForce_F<AFIELD> *fopr_force,
    AFprop<AFIELD> *fprop_MD, AFprop<AFIELD> *fprop_H)
    : Action(),
    m_fopr(fopr), m_fopr_force(fopr_force),
    m_fprop_MD(fprop_MD), m_fprop_H(fprop_H)
  {}

  Action_F_Standard_lex_alt(
    unique_ptr<AFopr<AFIELD> >& fopr,
    unique_ptr<AForce_F<AFIELD> >& fopr_force,
    unique_ptr<AFprop<AFIELD> >& fprop_MD,
    unique_ptr<AFprop<AFIELD> >& fprop_H)
    : Action(),
    m_fopr(fopr.get()), m_fopr_force(fopr_force.get()),
    m_fprop_MD(fprop_MD.get()), m_fprop_H(fprop_H.get())
  { init(); }

  ~Action_F_Standard_lex_alt()
  { tidyup(); }

  void set_parameters(const Parameters&);

  void set_label(const std::string label)
  {
    m_label = label;
    vout.detailed(m_vl, "  label: %s\n", m_label.c_str());
  }

  std::string get_label()
  { return m_label; }

  void set_config(Field *U);

  double langevin(RandomNumbers *);

  double calcH();

  void force(Field&);

 private:
  void init();
  void tidyup();
};
#endif
