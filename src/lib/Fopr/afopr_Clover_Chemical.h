/*!
        @file    afopr_Clover_Chemical.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-12-07 21:08:08 #$

        @version $LastChangedRevision: 2558 $
*/

#ifndef AFOPR_CLOVER_CHEMICAL_INCLUDED
#define AFOPR_CLOVER_CHEMICAL_INCLUDED

#include "lib/Fopr/afopr.h"

#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! Clover fermion operator with chemical potential.

/*!
    This class implements the Clover fermion operator at
    finite chemical potential.
    Former Fopr_Clover_Isochemical class was converted to
    a template class and renamed.
    The first version was implemented by H.M. [22 Aug 2012]
    and has been improved by Y.Namekawa and H.M.
                                  [27 Feb 2022 H.Matsufuru]
 */

template<typename AFIELD>
class AFopr_Clover_Chemical : public AFopr<AFIELD>
{
 public:
  static const std::string class_name;
  typedef typename AFIELD::real_t real_t;

 private:
  real_t m_kappa;               //!< hopping parameter
  real_t m_cSW;                 //!< clover coefficient
  real_t m_mu;                  //!< isospin chemical potential
  std::vector<int> m_boundary;  //!< boundary conditions
  Bridge::VerboseLevel m_vl;    //!< verbose level
  std::string m_repr;           //!< gamma-matrix representation

  std::string m_kernel_type;    //!< kernel type

  std::string m_mode;

  // internal data members
  int m_Nvol, m_Ndim;
  int m_Nc, m_Nd, m_NinF;
  real_t m_exp_mu;               //!< exp(mu)

  AFopr<AFIELD> *m_fopr_w;
  AFopr<AFIELD> *m_fopr_csw;

  const Field_G *m_U;

  AFIELD m_v1, m_v2;      //!< working field.

  bool m_is_initial_step; //!< to avoid redundant setup

 public:

  DEPRECATED
  AFopr_Clover_Chemical() { init("Dirac"); }

  DEPRECATED
  AFopr_Clover_Chemical(const std::string repr) { init(repr); }

  AFopr_Clover_Chemical(const Parameters& params) { init(params); }

  ~AFopr_Clover_Chemical() { tidyup(); }

  void set_parameters(const Parameters& params);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U);

  void set_mode(const std::string mode);

  void mult(AFIELD& v, const AFIELD& w);

  void mult_dag(AFIELD& v, const AFIELD& w);

  std::string get_mode() const { return m_mode; }

  void mult_gm5p(const int mu, AFIELD& v, const AFIELD& w);

  void mult_gm5(AFIELD&, const AFIELD&);
  void D(AFIELD&, const AFIELD&);
  void Dminmu(AFIELD&, const AFIELD&);

  void DdagD(AFIELD&, const AFIELD&);
  void Ddag(AFIELD&, const AFIELD&);
  void H(AFIELD&, const AFIELD&);
  void Hdag(AFIELD&, const AFIELD&);

  void mult_up(const int mu, AFIELD& v, const AFIELD& w)
  { m_fopr_w->mult_up(mu, v, w); }

  void mult_dn(const int mu, AFIELD& v, const AFIELD& w)
  { m_fopr_w->mult_dn(mu, v, w); }

  //! returns true if additional field conversion is needed.
  virtual bool needs_convert()
  { return m_fopr_w->needs_convert(); }

  //! converts a Field object into other format if necessary.
  virtual void convert(AFIELD& v, const Field& w)
  { m_fopr_w->convert(v, w); }

  //! reverses to a Field object from other format if necessary.
  virtual void reverse(Field& v, const AFIELD& w)
  { m_fopr_w->reverse(v, w); }

  int field_nin()  { return m_NinF; }
  int field_nvol() { return m_Nvol; }
  int field_nex()  { return 1; }

  //! this returns the number of floating point operations.
  double flop_count();

 private:
  void init(const std::string repr);

  void init(const Parameters& params);

  void tidyup();

  //! sets parameters given as values: private for composite operator.
  void set_parameters_impl(const real_t kappa,
                           const real_t cSW,
                           const real_t mu,
                           const std::vector<int> bc);

#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object()
  { return new AFopr_Clover_Chemical<AFIELD>(); }

  static AFopr<AFIELD> *create_object_with_repr(const std::string& repr)
  { return new AFopr_Clover_Chemical<AFIELD>(repr); }

  static AFopr<AFIELD> *create_object_with_params(const Parameters& params)
  { return new AFopr_Clover_Chemical<AFIELD>(params); }

 public:
  static bool register_factory()
  {
    bool                init = true;
    init &= AFopr<AFIELD>::Factory_noarg::Register("Clover_Chemical",
                                                   create_object);
    init &= AFopr<AFIELD>::Factory_string::Register("Clover_Chemical",
                                                    create_object_with_repr);
    init &= AFopr<AFIELD>::Factory_params::Register("Clover_Chemical",
                                                    create_object_with_params);
    return init;
  }
#endif
};
#endif
