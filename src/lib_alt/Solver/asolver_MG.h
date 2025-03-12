/*!

        @file    asolver_MG.h

        @brief   Multigrid solver

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate: 2023-02-28 16:09:41 +0900 (Tue, 28 Feb 2023) $

        @version $LastChangedRevision: 2492 $
*/
//====================================================================
#ifndef ASOLVER_MG_H
#define ASOLVER_MG_H

#include <cstdio>
#include <cstdlib>
#include <vector>
using std::vector;
#include <string>
using std::string;

#include "lib_alt/Solver/asolver.h"
#include "lib/Fopr/afopr.h"
#include "lib_alt/Fopr/afopr_dd.h"
#include "lib_alt/Solver/aprecond_MG.h"
#include "lib_alt/Field/afield_base.h"
#include "lib_alt/Solver/MultiGrid.h"

//====================================================================
template<typename AFIELD>
class ASolver_MG : public ASolver<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  using AFIELD_f = AField<float, AFIELD::IMPL>;
  // fopr types: should be given in the .cpp file
  //using FoprD_t = AFopr_Clover<AFIELD>;
  //using FoprF_t = AFopr_Clover_dd<AFIELD_f>;
  //using FoprCoarse_t = AFopr_Clover_coarse<AFIELD_f >;

  // solver types: should be given in the .cpp file
  //using OuterSolver_t  = ASolver_FBiCGStab<AFIELD>;
  //using CoarseSolver_t = ASolver_BiCGStab_Cmplx< AFIELD_f>;
  //using Smoother_t     = ASolver_SAP<AFIELD_f>;


  // index
  //  using MultiGrid_t = MultiGrid_Clover<AFIELD_f, AFIELD_f>;

  using ASolver<AFIELD>::m_vl;
  static const std::string class_name;

 protected:

  // preconditinor: I am the owener
  unique_ptr<APrecond_MG<AFIELD, AFIELD_f> > m_prec_mg;
  std::vector<int> m_sap_block_size;

  // fine grid fermion operators: only references
  AFopr<AFIELD> *m_afopr_fineD;         // assumes double prec.
  AFopr_dd<AFIELD_f> *m_afopr_fineF;    //  a block version is required
  AFopr_dd<AFIELD_f> *m_afopr_smoother; //  a block version is required

  // coarse grid fermion operator: I am the owner
  unique_ptr<AFopr<AFIELD_f> > m_afopr_coarse;

  // solvers: I am the owner
  unique_ptr<ASolver<AFIELD_f> > m_asolver_coarse;
  unique_ptr<ASolver<AFIELD_f> > m_asolver_smoother;
  unique_ptr<ASolver<AFIELD> > m_outer_solver;

  // multigrid converter: I am the owner
  unique_ptr<MultiGrid<AFIELD_f, AFIELD_f> > m_multigrid;

  int m_Niter;           //!< maximum iteration number.
  real_t m_Stop_cond;    //!< stopping criterion (squared).
  int m_Nconv;           //!< iteratoin number to calculate flop

  int m_nvec;            //!< number of testvectors
  int m_nsetup;          //!< setup iterations

  static constexpr int m_min_res_iter = 6;
  int m_smoother_niter;
  int m_smoother_stop_cond;

  Parameters m_params_asolver_coarse;  //!< parameters for coarse grid solver
  Parameters m_params_asolver_outer;   //!< parameters for outer solver
  //  Parameters m_params_fopr_dd;         //!< parameters for fopr, used to generated dd op.

  std::string m_mode; //! set "D" for colver (default), "DdagD" for domainwall

  //! to remember convergence iteration to provide flop count.
  int m_nconv;

  //! calling constructor without fermion operator is forbidden.
  //ASolver_MG(){ };

  //! working vectors.
  AFIELD m_vec_work;


 public:
  //! constructor.
  ASolver_MG()
    : m_Niter(0), m_Stop_cond(0.0L)
  {
    this->init();
  }

  //! destructor.
  ~ASolver_MG() { this->tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters by a Parameter object.
  void set_parameters_level0(const Parameters& params);

  //! setting parameters by a Parameter object.
  void set_parameters_level1(const Parameters& params);

  //! setting parameters.
  void set_parameters(const int Niter, const real_t Stop_cond, const std::vector<int>& sap_block_size, const int nvec);

  void set_parameters(const int Niter, const real_t Stop_cond, const std::string& outer_vlevel,
                      const std::vector<int>& sap_block_size, const int nvec, const int nsetup,
                      const int coarse_niter, const real_t coarse_stop_cond, const std::string& coarse_vlevel,
                      const int smoother_niter, const real_t smoother_stop_cond);

  //! setup
  void run_setup();
  void run_setup(std::vector<AFIELD_f>& testvec_work);
  void run_setup_initial(std::vector<AFIELD_f>& testvec_work);
  void run_setup_iterative(int niter, std::vector<AFIELD_f>& testvec_work);


  //! solver main.
  void solve(AFIELD& xq, const AFIELD& b, int& nconv, real_t& diff);

  //! returns the pointer to the fermion operator.
  AFopr<AFIELD> *get_fopr() { return m_afopr_fineD; }
  AFopr<AFIELD> *get_foprD() { return m_afopr_fineD; }
  AFopr_dd<AFIELD_f> *get_foprF() { return m_afopr_fineF; }
  AFopr_dd<AFIELD_f> *get_fopr_smoother() { return m_afopr_smoother; }

  AFopr_dd<AFIELD_f> *new_fopr_smoother(const Parameters& param) const;

  AFopr<AFIELD_f> *get_fopr_coarse() { return m_afopr_coarse.get(); }

  // ! set fopr
  void set_foprD(AFopr<AFIELD> *foprD);

  void set_foprF(AFopr_dd<AFIELD_f> *foprF);

  void set_fopr_smoother(AFopr_dd<AFIELD_f> *foprF);

  void set_lattice(const vector<int>& sap_block_size);

  // ! initialize internal solvers
  void init_solver(std::string mode);

  // ! initialize internal solvers
  void init_solver()
  {
    init_solver("D");
  }

  //! returns the floating point operation count.
  double flop_count();

  //! returns the floating point operation count [setup].
  double flop_count_setup();

  void reset_flop_count();

  bool use_fopr_for_smoother() const;

 protected:

  void init_coarse_grid();

  void init(void);

  void tidyup(void);

  void solve_MG_init(real_t& rrp, real_t& rr);

  void solve_MG_step(real_t& rrp, real_t& rr);

  std::vector<int> m_coarse_lattice;

  unique_ptr<Timer> m_timer_gramschmidt;
  unique_ptr<Timer> m_timer_generate_coarse_op;
};

#endif // include gurad
//============================================================END=====
