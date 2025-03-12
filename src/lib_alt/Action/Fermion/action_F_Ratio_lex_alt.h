/*!
      @file    action_F_Ratio_lex_alt.h
      @brief
      @author  Yusuke Namekawa (namekawa)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-02-28 16:09:41 #$
      @version $LastChangedRevision: 2492 $
*/

#ifndef ACTION_F_RATIO_LEX_ALT_INCLUDED
#define ACTION_F_RATIO_LEX_ALT_INCLUDED

#include "lib/Action/action.h"

#include "lib/Force/Fermion/force_F.h"
#include "lib/Measurements/Fermion/fprop.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

// include files in alt-code
//#include "Field/afield.h"
#include "lib_alt/Measurements/Fermion/afprop.h"

//! HMC action for Hasenbusch preconditioned fermions.

/*!
    Hasenbusch preconditioned fermion action in alternative
    implementation.
                                    [04 Feb 2019 H.Matsufuru]
 */

template<typename AFIELD>
class Action_F_Ratio_lex_alt : public Action
{
 public:
  //  typedef AField<double> AFIELD;
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  Field *m_U;

  AFopr<AFIELD> *m_fopr_prec;          // preconditioner
  AForce_F<AFIELD> *m_fopr_prec_force; // force of preconditioner
  AFopr<AFIELD> *m_fopr;               // dynamical fermion
  AForce_F<AFIELD> *m_fopr_force;      // force of dynamical fermion
  AFIELD m_psf;                        // pseudofermion field
  std::string m_label;                 // label of action

  AFprop<AFIELD> *m_fprop_H_prec;
  AFprop<AFIELD> *m_fprop_MD;
  AFprop<AFIELD> *m_fprop_H;

 public:
  Action_F_Ratio_lex_alt(
    AFopr<AFIELD> *fopr_prec, AForce_F<AFIELD> *fopr_prec_force,
    AFopr<AFIELD> *fopr, AForce_F<AFIELD> *fopr_force,
    AFprop<AFIELD> *fprop_H_prec,
    AFprop<AFIELD> *fprop_MD, AFprop<AFIELD> *fprop_H)
    : Action(),
    m_fopr_prec(fopr_prec), m_fopr_prec_force(fopr_prec_force),
    m_fopr(fopr), m_fopr_force(fopr_force),
    m_fprop_H_prec(fprop_H_prec),
    m_fprop_MD(fprop_MD), m_fprop_H(fprop_H)
  { init(); }

  ~Action_F_Ratio_lex_alt()
  { tidyup(); }

  void set_parameters(const Parameters&);
  void set_parameters();

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
