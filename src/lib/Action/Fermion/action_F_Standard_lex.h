/*!
        @file    action_F_Standard_lex.h

        @brief

        @author  UEDA, Satoru  (sueda)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef ACTION_F_STANDARD_LEX_INCLUDED
#define ACTION_F_STANDARD_LEX_INCLUDED

#include "Action/action.h"

#include "Force/Fermion/force_F.h"
#include "Measurements/Fermion/fprop.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Standard fermion action for HMC.

/*!
    This class is used to define an action used in HMC.
    Fermion and Force operators are given at the construction.
                                        [05 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    Fprop is introduced.                [03 Mar 2013 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                        [21 Mar 2015 Y.Namekawa]
 */

class Action_F_Standard_lex : public Action
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  Fopr *m_fopr;
  Force *m_fopr_force;
  Field m_psf;
  std::string m_label;

  Fprop *m_fprop_MD;
  Fprop *m_fprop_H;

  Field *m_U;


 public:
  Action_F_Standard_lex(
    Fopr *fopr, Force *fopr_force,
    Fprop *fprop_MD, Fprop *fprop_H)
    : m_vl(CommonParameters::Vlevel()),
    m_fopr(fopr), m_fopr_force(fopr_force),
    m_fprop_MD(fprop_MD), m_fprop_H(fprop_H)
  {}

  Action_F_Standard_lex(
    Fopr *fopr, Force *fopr_force,
    Fprop *fprop_MD, Fprop *fprop_H,
    const Parameters& params)
    : m_vl(CommonParameters::Vlevel()),
    m_fopr(fopr), m_fopr_force(fopr_force),
    m_fprop_MD(fprop_MD), m_fprop_H(fprop_H)
  {
    set_parameters(params);
  }

  ~Action_F_Standard_lex()
  {
    // delete m_fprop;
  }

  void set_parameters(const Parameters&);

  void get_parameters(Parameters&) const;

  void set_label(const std::string label)
  {
    m_label = label;
    vout.detailed(m_vl, "  label: %s\n", m_label.c_str());
  }

  std::string get_label()
  {
    return m_label;
  }

  void set_config(Field *U);

  double langevin(RandomNumbers *);

  double calcH();

  void force(Field&);
};
#endif
