/*!
        @file    force_F_Domainwall.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#ifndef FORCE_F_DOMAINWALL_INCLUDED
#define FORCE_F_DOMAINWALL_INCLUDED

#include "force_F_Wilson_Nf2.h"

#include "Fopr/fopr_Domainwall.h"
#include "Field/index_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Force calculation for domain-wall fermions.

/*!
    At present, only the standard domain-wall setting is
    available.
                                [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
 */


class Force_F_Domainwall : public Force
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  Index_lex *m_index;

  // parameters common to overlap fermion
  double m_mq;             // quark mass
  double m_M0;             // domain-wall height
  int m_Ns;                // size of fifth-dimension
  std::vector<int> m_boundary;
  std::vector<double> m_b; //!< coefficient b (array)
  std::vector<double> m_c; //!< coefficient c (array)

  Fopr_Wilson *m_fopr_w;
  Fopr_Domainwall *m_fopr_dw;
  Force_F_Wilson_Nf2 *m_force_w;

 public:
  Force_F_Domainwall()
    : m_vl(CommonParameters::Vlevel())
  {
    init();
  }

  Force_F_Domainwall(const Parameters& params)
    : m_vl(CommonParameters::Vlevel())
  {
    init();
    set_parameters(params);
  }

  ~Force_F_Domainwall()
  {
    tidyup();
  }

 private:
  void init()
  {
    m_fopr_w  = new Fopr_Wilson;
    m_fopr_dw = new Fopr_Domainwall(m_fopr_w);
    m_force_w = new Force_F_Wilson_Nf2;
    m_boundary.resize(CommonParameters::Ndim());
  }

  void tidyup()
  {
    delete m_force_w;
    delete m_fopr_dw;
    delete m_fopr_w;
  }

 public:

  void set_parameters(const Parameters& params);

  void set_parameters(const double mq,
                      const double M0,
                      const int Ns,
                      const std::vector<int> bc,
                      const double b,
                      const double c);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
    m_fopr_w->set_config(U);
    m_fopr_dw->set_config(U);
    m_force_w->set_config(U);
  }

  void force_udiv(Field& force, const Field& eta);

  void force_core1(Field& force, const Field& zeta, const Field& eta); // override default
  void force_udiv1(Field& force, const Field& zeta, const Field& eta);
};
#endif
