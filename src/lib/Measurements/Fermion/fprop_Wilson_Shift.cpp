/*!
        @file    fprop_Wilson_Shift.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#include "fprop_Wilson_Shift.h"

const std::string Fprop_Wilson_Shift::class_name = "Fprop_Wilson_Shift";

//====================================================================
void Fprop_Wilson_Shift::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  int                 Nshift;
  std::vector<double> sigma;
  int                 Niter;
  double              Stop_cond;

  int err = 0;
  err += params.fetch_int("number_of_shifts", Nshift);
  err += params.fetch_double_vector("shifted_mass_difference", sigma);
  err += params.fetch_int("maximum_number_of_iteration", Niter);
  err += params.fetch_double("convergence_criterion_squared", Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(Nshift, sigma, Niter, Stop_cond);
}


//====================================================================
void Fprop_Wilson_Shift::get_parameters(Parameters& params) const
{
  params.set_int("number_of_shifts", m_Nshift);
  params.set_double_vector("shifted_mass_difference", m_sigma);
  params.set_int("maximum_number_of_iteration", m_Niter);
  params.set_double("convergence_criterion_squared", m_Stop_cond);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Fprop_Wilson_Shift::set_parameters(const int Nshift, const std::vector<double> sigma,
                                        const int Niter, const double Stop_cond)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Nshift    = %d\n", Nshift);
  for (int i = 0; i < Nshift; ++i) {
    vout.general(m_vl, "  sigma[%d] = %16.8e\n", i, sigma[i]);
  }
  vout.general(m_vl, "  Niter     = %d\n", Niter);
  vout.general(m_vl, "  Stop_cond = %8.2e\n", Stop_cond);

  //- range check
  // NB. Nshift,sigma == 0 is allowed.
  int err = 0;
  err += ParameterCheck::non_zero(Niter);
  err += ParameterCheck::square_non_zero(Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_Nshift = Nshift;
  m_sigma.resize(Nshift);
  m_sigma = sigma;

  m_Niter     = Niter;
  m_Stop_cond = Stop_cond;
}


//====================================================================
double Fprop_Wilson_Shift::invert_D(std::vector<Field_F> *xq2, const Field_F& b)
{
  const int Nin  = b.nin();
  const int Nvol = b.nvol();
  const int Nex  = b.nex();

  const int Nshift = m_Nshift;

  std::vector<double> sigma = m_sigma;

  std::vector<Field> xq(Nshift);

  for (int i = 0; i < Nshift; ++i) {
    xq[i].reset(Nin, Nvol, Nex);
  }

  int    Nconv = 0;
  double diff  = 1.0;


  //  vout.general(m_vl, "size of xq = %d\n", xq->size());
  //  vout.general(m_vl, "size of xq[0] = %d\n", xq[0].nvol());

  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Number of shift values = %d\n", sigma.size());
  assert(xq2->size() == sigma.size());

  m_fopr->set_mode("DdagD");

  //  std::vector<Field_F>* xq2 = new std::vector<Field_F>;

  /*
  std::vector<Field_F>* xq2;
  xq2->resize(xq->size());
  for(int i=0; i<xq->size(); ++i){
    xq2[i] = (Field*) xq[i];
  }
  */

  Shiftsolver_CG solver(m_fopr, m_Niter, m_Stop_cond);
  solver.solve(xq, sigma, (Field)b, Nconv, diff);

  vout.general(m_vl, "  residues of solutions(2):\n");

  // Field version: works
  Field s((Field)b);
  Field x((Field)b);
  Field t((Field)b);

  double diff1 = 1.0;  // superficial initialization
  for (int i = 0; i < Nshift; ++i) {
    x = xq[i];
    double sg = sigma[i];
    m_fopr->mult(s, x);
    axpy(s, sg, x);
    axpy(s, -1.0, t);
    diff1 = dot(s, s);

    vout.general(m_vl, "i_shift,diff = %6d  %22.15e\n", i, diff1);
  }

  for (int i = 0; i < Nshift; ++i) {
    (*xq2)[i] = (Field_F)xq[i];
  }

  const double result = diff1;

  return result;
}
