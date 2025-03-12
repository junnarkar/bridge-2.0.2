/*!
        @file    fopr_Wilson_impl.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#ifndef FOPR_WILSON_IMPL_IMP_INCLUDED
#define FOPR_WILSON_IMPL_IMP_INCLUDED

//#include "Fopr/fopr_Wilson.h"
#include "Fopr/fopr.h"

#include "Field/shiftField_lex.h"
#include "Tools/gammaMatrixSet.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Wilson fermion operator.

/*!
    This fermion operator defines the standard Wilson fermion.
    The gamma matrix representation is given as a control string
    "Dirac"(default) or "Chiral" at the construction, which is
    used to construct the Fopr_Wilson instance.
    The `mode', which of D, Ddag, H, DdagD are multiplied, is
    controlled by setting the pointers to these functions,
    m_mult and m_mult_dag.
    At the beginning, they are set to point mult_undef() which
    just represent the mode has not been set.
    set_mode(string) must be called before mult() is called.
                                    [24 Dec 2011 H.Matsufuru]
 */

namespace Imp {
  class Fopr_Wilson : public Fopr
  {
   public:
    static const std::string class_name;

   private:
    // input parameters
    double m_kappa;              //!< hopping parameter
    std::vector<int> m_boundary; //!< boundary condition
    std::string m_repr;          //!< gamma-matrix representation
    Bridge::VerboseLevel m_vl;   //!< verbose level

    std::string m_mode;          //!< mult mode

    // internal data members
    int m_Nc, m_Nd, m_Nvc, m_Ndf;
    int m_Nx, m_Ny, m_Nz, m_Nt;
    int m_Nvol, m_Ndim;

    const Field_G *m_U;                       //!< gauge configuration.

    std::vector<double> m_boundary_each_node; //!< b.c. on each node.

    //! arrays for communication buffer.
    double *vcp1_xp, *vcp2_xp, *vcp1_xm, *vcp2_xm;
    double *vcp1_yp, *vcp2_yp, *vcp1_ym, *vcp2_ym;
    double *vcp1_zp, *vcp2_zp, *vcp1_zm, *vcp2_zm;
    double *vcp1_tp, *vcp2_tp, *vcp1_tm, *vcp2_tm;

    Field m_w1, m_w2;   //!< working fields

   public:
    //! standard constructor.
    Fopr_Wilson(const Parameters& params) { init(params); }

    DEPRECATED
    Fopr_Wilson() { init("Dirac"); }

    DEPRECATED
    Fopr_Wilson(const std::string repr) { init(repr); }

    ~Fopr_Wilson() { tidyup(); }

    void set_parameters(const Parameters& params);

    void set_parameters(const double kappa, const std::vector<int> bc);

    void get_parameters(Parameters& params) const;

    void set_config(Field *U);

    void set_mode(const std::string mode);

    std::string get_mode() const { return m_mode; }

    void mult(Field& v, const Field& w);

    void mult_dag(Field& v, const Field& w);

    void mult(Field& v, const Field& w, const std::string mod);

    void mult_dag(Field& v, const Field& w, const std::string mode);

    void mult_gm5(Field& v, const Field& w);

    void D(Field& v, const Field& w);

    void Ddag(Field& v, const Field& w);

    void DdagD(Field& v, const Field& w);

    void DDdag(Field& v, const Field& w);

    void H(Field& v, const Field& w);

    void D_ex(Field& v, const int ex1,
              const Field& f, const int ex2);

    void mult_gm5p(const int mu, Field&, const Field&);

    void proj_chiral(Field& w, const int ex1,
                     const Field& v, const int ex2, const int ipm);

    void mult_up(const int mu, Field&, const Field&);
    void mult_dn(const int mu, Field&, const Field&);

    int field_nvol() { return m_Nvol; }
    int field_nin()  { return m_Nvc * m_Nd; }
    int field_nex()  { return 1; }

    double flop_count();

   private:
    // prohibit copy
    Fopr_Wilson(const Fopr_Wilson&) {}
    Fopr_Wilson& operator=(const Fopr_Wilson&);

    //! to be discarded.
    void init(const std::string repr);

    //! to be discarded.
    void init();

    //! standard initial setup.
    void init(const Parameters& params);

    //! initial setup main.
    void setup();

    //! final clean-up.
    void tidyup();

    void mult_gm5_chiral(Field&, const Field&);
    void mult_gm5_dirac(Field&, const Field&);

    void D_ex_chiral(Field&, const int ex1, const Field&, const int ex2);
    void D_ex_dirac(Field&, const int ex1, const Field&, const int ex2);

    void D_ex_chiral_alt(Field&, const int ex1, const Field&, const int ex2);
    void D_ex_dirac_alt(Field&, const int ex1, const Field&, const int ex2);

    void mult_xp(Field&, const Field&);
    void mult_xm(Field&, const Field&);
    void mult_yp(Field&, const Field&);
    void mult_ym(Field&, const Field&);
    void mult_zp(Field&, const Field&);
    void mult_zm(Field&, const Field&);

    void mult_tp_dirac(Field&, const Field&);
    void mult_tm_dirac(Field&, const Field&);
    void mult_tp_chiral(Field&, const Field&);
    void mult_tm_chiral(Field&, const Field&);

    void daypx(Field&, const double, const Field&);
    void clear(Field&);


#ifdef USE_FACTORY
   private:
    static Fopr *create_object() { return new Fopr_Wilson(); }

    static Fopr *create_object_with_repr(const std::string& repr)
    { return new Fopr_Wilson(repr); }

    static Fopr *create_object_with_params(const Parameters& params)
    { return new Fopr_Wilson(params); }

   public:
    static bool register_factory()
    {
      bool init = true;
      init &= Fopr::Factory_noarg::Register("Wilson/Imp", create_object);
      init &= Fopr::Factory_string::Register("Wilson/Imp",
                                             create_object_with_repr);
      init &= Fopr::Factory_params::Register("Wilson/Imp",
                                             create_object_with_params);
      return init;
    }
#endif
  };
}
#endif /* FOPR_WILSON_IMPL_IMP_INCLUDED */
