/*!
        @file    builder_Integrator.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef BUILDER_INTEGRATOR_INCLUDED
#define BUILDER_INTEGRATOR_INCLUDED

#include <cassert>

#include "action_list.h"
#include "integrator.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Builder of MD integrator for HMC.

/*!
    This class builds a molecular dynamics integrator for hybrid Monte
    Carlo update algorithm by composing subclasses of Integrator.
    At present standard leapfrog and omelyan integrators can be
    composed of.
                                     [25 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.
                                     [03 Mar 2013 Y.Namekawa]
 */


class Builder_Integrator
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  std::vector<int> m_Nstep;            //!< Number of steps at each level

  int m_Nprec;                         //!< precision parameter of exponentiation
  double m_lambda_Omelyan;

  std::string m_str_integrator_type;

  ActionList m_actions;
  std::vector<Director *> m_director;


  std::vector<Integrator *> m_integs;  //!< Integrator to be constructed

 public:
  //! constructor with ActionList
  Builder_Integrator(const ActionList& action_list,
                     std::vector<Director *> director = std::vector<Director *>());

  Builder_Integrator(const ActionList& action_list,
                     std::vector<Director *> director,
                     const Parameters& params);

  Builder_Integrator(const ActionList& action_list,
                     const Parameters& params);

  //! destructor
  ~Builder_Integrator()
  {
    tidyup();
  }

 private:
  // non-copyable
  Builder_Integrator(const Builder_Integrator&);
  Builder_Integrator& operator=(const Builder_Integrator&);

 public:
  void set_parameters(const Parameters& params);
  void set_parameters(const std::string str_integrator_type,
                      const std::vector<int>& Nstep,
                      const int Nprec,
                      const double lambda_Omelyan);

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  void get_parameters(Parameters& params) const;

  Integrator *build();

  Integrator *build_leapfrog();
  Integrator *build_omelyan();

  void tidyup();
};
#endif
