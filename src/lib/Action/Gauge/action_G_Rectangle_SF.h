/*!
        @file    action_G_Rectangle_SF.h

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-01-26 21:43:53 #$

        @version $LastChangedRevision: 2457 $
*/


#ifndef ACTION_G_RECTANGLE_SF_INCLUDED
#define ACTION_G_RECTANGLE_SF_INCLUDED

#include "Action/action.h"
#include "Force/Gauge/force_G_Rectangle_SF.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! HMC action class for rectangular gauge action with the SF BC.

/*!
  Gauge action with plaquette and rectangular Wilson loops.
  Iwasaki, Luscher-Weisz, DBW2 are examples of this type
  of action.

  SF boundary condition is intrduced accrdong to the following policy.
  <ul>
    <li>The total size of temporal links is Nt with the Dirichlet BC.
    <li>We need t=0, ..., Nt sites for Nt links.
    <li>The boundary spatial link at t=0 is dummy.
    <ul>
      <li>Calculation of the Hamiltonian and force is overridden to given the correct result without using the spatial link at t=0 stored in the Field_G object.
      <li>The spatial link at t=0 is updated by a randomly given conjugate momentum but is not used for the Hamiltoniand nor the force.
      <li>The corresponding force for the spatial conjugate momentum at t=0 is set to zero. This conjugate momentum is not updated at the boundary.
    </ul>
    <li>The boundary spatial sites at t=Nt does not exist in the Field_G object.
    <ul>
      <li>Calculation of the Hamiltonian and force is overridden to given the correct result using an approriate matrix for the spatial link at t=Nt.
    </ul>
  </ul>
  The rectangle improved gauge action is given by
\f[
S[U]=
\frac{\beta}{N_c}\sum_{{\cal C}\in{\cal S}_0}W_0({\cal C})
{\rm Re}{\rm Tr}\left(-P({\cal C})\right)
+\frac{\beta}{N_c}\sum_{{\cal C}\in{\cal S}_1}W_1({\cal C})
{\rm Re}{\rm Tr}\left(-R({\cal C})\right)
\f]
where \f${\cal C}\f$ is an oriented plaquette or rectangle.
One needs to choose the weight factors appropriately
to achieve the O(a) improvement.
\f[
W_0({\cal C}) =
      \left\{
       \begin{array}{ll}
         c_0 c^P_{\rm{t}}(g^2_0)
       & \mbox{for } {\cal C} \in P_{\rm{t}} :
          \mbox{Set of temporal plaquettes that just touch}
       \\
       &
          \mbox{\hspace{19mm} one of the boundaries, }
       \\
         c_0
       & \mbox{for } {\cal C} \in P_{\rm{other}} :
          \mbox{otherwise, }
       \end{array}
            \right.
\f]
\f[
W_1({\cal C})
      =
      \left\{
       \begin{array}{ll}
         c_1 c^R_{\rm{t}}(g^2_0)
       & \mbox{for } {\cal C} \in R^2_{\rm{t}} :
          \mbox{Set of temporal rectangles that have exactly}
       \\
       &
          \mbox{\hspace{19mm} two links on a boundary, }
       \\
         c_1
       & \mbox{for } {\cal C} \in R_{\rm{other}} :
          \mbox{otherwise, }
       \end{array}
            \right.
\f]
  <ul>
  <li>A major difference from Action_G_Rectangle class is to override a calculation of the Hamltonian and the force.
  <li>Boundary condition can be accessed with m_phi and m_phipr.
  <li>Boundary improvement factor is stored in m_ct, m_ctr.
  <li> [03 Feb. 2012 Y.Taniguchi]
  </ul>
    (Coding history will be recovered from trac.)
    YAML is implemented.           [14 Nov 2012 Y.Namekawa]
 */


class Action_G_Rectangle_SF : public Action
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  int m_Nc;

  double m_beta;
  double m_c_plaq;
  double m_c_rect;

  std::string m_label;

  Field_G *m_U;
  Staple_SF m_staple;
  ShiftField_lex m_shift;
  Force_G *m_force_G;

  //! SF boundary condition
  Mat_SU_N m_wk, m_wkpr;
  //  double *m_phi;
  //  double *m_phipr;
  double m_ct;  //!< SF boundary improvement coefficient for the plaquatte action
  double m_ctr; //!< SF boundary improvement coefficient for the rectangle action

  // raw parameter values
  std::vector<double> m_phi;
  std::vector<double> m_phipr;

 public:
  Action_G_Rectangle_SF()
    : m_vl(CommonParameters::Vlevel()),
    m_Nc(CommonParameters::Nc()), m_wk(m_Nc), m_wkpr(m_Nc)
  {
    m_force_G = Force_G::New("Force_G_Rectangle_SF");
  }

  Action_G_Rectangle_SF(const Parameters& params)
    : m_vl(CommonParameters::Vlevel()),
    m_Nc(CommonParameters::Nc()), m_wk(m_Nc), m_wkpr(m_Nc)
  {
    m_force_G = new Force_G_Rectangle_SF(params);
    set_parameters(params);
  }

  ~Action_G_Rectangle_SF()
  {
    delete m_force_G;
  }

  void set_parameters(const Parameters& params);

  void set_parameters(const double beta, const double c_plaq, const double c_rect,
                      double *phi, double *phipr,
                      const double ct, double ctr);

  void get_parameters(Parameters& params) const;

  void set_label(const std::string label)
  {
    m_label = label;
    vout.detailed(m_vl, "  label: %s\n", m_label.c_str());
  }

  std::string get_label()
  {
    return m_label;
  }

  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
  }

  double langevin(RandomNumbers *);

  double calcH();

  void force(Field&);

#ifdef USE_FACTORY
 private:
  static Action *create_object()
  {
    return new Action_G_Rectangle_SF();
  }

  static Action *create_object_with_params(const Parameters& params)
  {
    return new Action_G_Rectangle_SF(params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Action::Factory::Register("Action_G_Rectangle_SF", create_object);
    init &= Action::Factory_params::Register("Action_G_Rectangle_SF", create_object_with_params);
    return init;
  }
#endif
};

#endif
