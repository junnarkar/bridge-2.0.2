/*!
        @file    corr2pt_4spinor.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#include "corr2pt_4spinor.h"

#include "Tools/epsilonTensor.h"

const std::string Corr2pt_4spinor::class_name = "Corr2pt_4spinor";

//====================================================================
void Corr2pt_4spinor::set_parameters(const Parameters& params)
{
  m_filename_output = params.get_string("filename_output");
  if (m_filename_output.empty()) {
    m_filename_output = "stdout";
  }

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }
}


//====================================================================
void Corr2pt_4spinor::get_parameters(Parameters& params) const
{
  params.set_string("filename_output", m_filename_output);
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Corr2pt_4spinor::init()
{
  assert(CommonParameters::Nc() == 3);

  m_filename_output = "stdout";
}


//====================================================================
double Corr2pt_4spinor::pion(const std::vector<Field_F>& sq1, 
                             const std::vector<Field_F>& sq2)
{
  const int Lt = CommonParameters::Lt();
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  std::vector<dcomplex> corr(Lt);

  GammaMatrix unity  = m_gmset->get_GM(m_gmset->UNITY);
  GammaMatrix gm5  = m_gmset->get_GM(m_gmset->GAMMA5);
  GammaMatrix gm_src  = m_gmset->get_GM(m_gmset->GAMMA5);
  GammaMatrix gm_sink = m_gmset->get_GM(m_gmset->GAMMA5);
  GammaMatrix parity = gm5.add(unity); 
  
  vout.general(m_vl, "PS <-- PS correlator:\n");
  
  //meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  //pion_correlator(corr, gm_sink, gm_src, sq1, sq2);
  //pion_test(corr, gm_sink, gm_src, sq1, sq2);
  //nucleon_correlator(corr, sq1);
  nucleon_test(corr, sq1);
  //nucleon_test1(corr, sq1);
  //nucleon_test2(corr, sq1);
  //proton_correlator(corr, parity, sq1, sq2);
  //pion_modsq(corr, gm_sink, gm_src, sq1, sq2);
  //pion_correlator_modsq(corr, gm_sink, gm_src, sq1, sq2);
  
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }
  
  const double result = real(corr[0]);
  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }

  return result;
}

//====================================================================

void Corr2pt_4spinor::pion_correlator(std::vector<dcomplex>& corr_global,
                                       const GammaMatrix& gm_sink,
                                       const GammaMatrix& gm_src,
                                       const std::vector<Field_F>& sq1,
                                       const std::vector<Field_F>& sq2)
{
  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);

  const int Nvol   = sq1[0].nvol();
  const int Nvol_s = Nvol / CommonParameters::Nt();

  const GammaMatrix gm5         = m_gmset->get_GM(m_gmset->GAMMA5);
  const GammaMatrix gm_gm5_src  = gm_src.mult(gm5);
  const GammaMatrix gm5_gm_sink = gm5.mult(gm_sink);

  std::vector<dcomplex> corr_local(Nt, cmplx(0.0, 0.0));

  int id1[Nd];
  int id2[Nd];
  int id_src[Nd];

  int NC2 = 6 ;
  int NCD2 = 24;

  for (int id = 0; id < Nd; ++id) {
    id_src[id] = gm_gm5_src.index(id);
    id1[id] = id * NC2;
    id2[id] = gm5_gm_sink.index(id) * NC2;
  }
  
  dcomplex corr_t;

  for (int t = 0; t < Nt; ++t) {
    corr_t = 0 ;
    for (int ss = 0; ss < Nvol_s; ++ss) {

    int site = NCD2 * (ss + t * Nvol_s);
    
    for(int c0=0; c0 < Nc; ++c0){
      for(int d0=0; d0 < Nd; ++d0){
        
        const double *u_quark, *d_quark;

        u_quark = sq1[c0 + Nc * id_src[d0]].ptr(0);
        d_quark = sq2[c0 + Nc * d0].ptr(0);

        for (int c1 = 0; c1 < Nc; ++c1) {
          for (int d1 = 0; d1 < Nd; ++d1) {

            dcomplex u_prop, d_prop;

            int ic1_r = 2 * c1 + id1[d1] + site;
            int ic1_i = 2 * c1 + 1 + id1[d1] + site;
            u_prop = cmplx(u_quark[ic1_r], u_quark[ic1_i]);

            int ic2_r = 2 * c1 + id2[d1] + site;
            int ic2_i = 2 * c1 + 1 + id2[d1] + site;
            d_prop = cmplx(d_quark[ic2_r],d_quark[ic2_i]);

            corr_t += gm_gm5_src.value(d0) * gm5_gm_sink.value(d1) * conj(u_prop) * d_prop ;

            }
          }
        }
      }
    }
  corr_local[t] = corr_t;
  }

  global_corr_t(corr_global, corr_local);

}

void Corr2pt_4spinor::pion_corr(std::vector<dcomplex>& corr_global,
                                       const GammaMatrix& gm_sink,
                                       const GammaMatrix& gm_src,
                                       const std::vector<Field_F>& sq1,
                                       const std::vector<Field_F>& sq2)
{
  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);

  const int Nvol   = sq1[0].nvol();
  const int Nvol_s = Nvol / CommonParameters::Nt();
  
  const GammaMatrix gm5         = m_gmset->get_GM(m_gmset->GAMMA5);
  const GammaMatrix gm_gm5_src  = gm_src.mult(gm5);
  const GammaMatrix gm5_gm_sink = gm5.mult(gm_sink);

  std::vector<dcomplex> corr_local(Nt, cmplx(0.0, 0.0));
  
  int id1[Nd];
  int id2[Nd];
  int id_src[Nd];

  int NC2 = 6 ;
  int NCD2 = 24;

  for (int id = 0; id < Nd; ++id) {
    id_src[id] = gm_gm5_src.index(id);
    id1[id] = id * NC2;
    id2[id] = gm5_gm_sink.index(id) * NC2;
  }

  dcomplex corr;

  for (int c0 = 0; c0 < Nc; ++c0) {
    for (int d0 = 0; d0 < Nd; ++d0) {
      
      const double *w1 = sq1[c0 + Nc * d0].ptr(0);
      const double *w2 = sq2[c0 + Nc * id_src[d0]].ptr(0);
      
      for (int t = 0; t < Nt; ++t) {
        double c_r[Nd], c_i[Nd];
        
        for (int id = 0; id < Nd; ++id) {
          c_r[id] = 0.0;
          c_i[id] = 0.0;
        }
        
        for (int ss = 0; ss < Nvol_s; ++ss) {
          int site = NCD2 * (ss + t * Nvol_s);

          for (int c1 = 0; c1 < Nc; ++c1) {
            for (int d1 = 0; d1 < Nd; ++d1) {
              int ic1_r = 2 * c1 + id1[d1] + site;
              int ic2_r = 2 * c1 + id2[d1] + site;

              int ic1_i = 2 * c1 + 1 + id1[d1] + site;
              int ic2_i = 2 * c1 + 1 + id2[d1] + site;

              c_r[d1] += w1[ic2_r] * w2[ic1_r] + w1[ic2_i] * w2[ic1_i];

              c_i[d1] += -w1[ic2_r] * w2[ic1_i] + w1[ic2_i] * w2[ic1_r];

              //corr_t += gm_gm5_src.value(d0) * gm5_gm_sink.value(d1) * cmplx(c_r, c_i);
            }
          }
        }

      //corr = cmplx(0.0, 0.0);
      //for (int id = 0; id < Nd; ++id) {
      //corr += gm5_gm_sink.value(id) * cmplx(c_r[id], c_i[id]);
     // }
     
      corr_local[t] += gm_gm5_src.value(d0) * corr;
      
      }
    }
  }
  global_corr_t(corr_global, corr_local);
}


void Corr2pt_4spinor::pion_modsq(std::vector<dcomplex>& corr_global,
                                       const GammaMatrix& gm_sink,
                                       const GammaMatrix& gm_src,
                                       const std::vector<Field_F>& sq1,
                                       const std::vector<Field_F>& sq2)
{
  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);

  const int Nvol   = sq1[0].nvol();
  const int Nvol_s = Nvol / CommonParameters::Nt();
  

  std::vector<dcomplex> corr_local(Nt, cmplx(0.0, 0.0));
  

  int NC2 = 6 ;
  int NCD2 = 24;
      
  for (int t = 0; t < Nt; ++t) {
    double c_r=0.0;
    double c_i=0.0;
    for (int ss = 0; ss < Nvol_s; ++ss) {
      int site = NCD2 * (ss + t * Nvol_s);
      
      for (int c0 = 0; c0 < Nc; ++c0) {
        for (int d0 = 0; d0 < Nd; ++d0) {
      
          const double *w1 = sq1[c0 + Nc * d0].ptr(0);

          for (int c1 = 0; c1 < Nc; ++c1) {
            for (int d1 = 0; d1 < Nd; ++d1) {
              int ic1_r = 2 * c1 + (NC2*d1) + site;
              int ic1_i = ic1_r + 1;

              c_r += (w1[ic1_r] * w1[ic1_r]) + (w1[ic1_i] * w1[ic1_i]);

              c_i = 0.0; 

            }
          }
        }
      }
    }
    corr_local[t] += cmplx(c_r,c_i);
  }
  global_corr_t(corr_global, corr_local);
}

//====================================================================
void Corr2pt_4spinor::nucleon_correlator(std::vector<dcomplex>& corr_global,
                                       const std::vector<Field_F>& sq1)
{
  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);

  const int Nvol   = sq1[0].nvol();
  const int Nvol_s = Nvol / CommonParameters::Nt();

  const GammaMatrix unity       = m_gmset->get_GM(m_gmset->UNITY);
  const GammaMatrix gm5       = m_gmset->get_GM(m_gmset->GAMMA5);
  const GammaMatrix C         = m_gmset->get_GM(m_gmset->CHARGECONJG);
  const GammaMatrix Cg5_src  = C.mult(gm5);
  const GammaMatrix Cg5_snk  = C.mult(gm5);
  //const GammaMatrix parity_src   = gm5.add(unity);
  //const GammaMatrix parity_snk   = gm5.add(unity);
  const GammaMatrix parity_src   = unity;
  const GammaMatrix parity_snk   = unity;

  EpsilonTensor e_src, e_snk;

  std::vector<dcomplex> corr_local(Nt, cmplx(0.0, 0.0));

  int NC2 = 6 ;
  int NCD2 = 24;
  int Nfaces = 6;
  
  int id1_src[Nd],id2_src[Nd],id3_src[Nd], id1_snk[Nd],id2_snk[Nd],id3_snk[Nd];

  for (int id = 0; id < Nd; ++id) {
    id1_src[id] = id ;
    id2_src[id] = parity_src.index(id) ;
    id3_src[id] = Cg5_src.index(id) ;
    id1_snk[id] = id * NC2;
    id2_snk[id] = parity_snk.index(id) * NC2;
    id3_snk[id] = Cg5_snk.index(id) * NC2;
  }
  
  const double *u1_quark, *d_quark, *u2_quark;

  for (int beta = 0; beta < Nd; ++beta){
    for (int mu = 0; mu < Nd; ++mu){
      for (int nu = 0; nu < Nd; ++nu){
        for (int c_src = 0; c_src < Nfaces; ++c_src){

          int a  = e_src.epsilon_3_index(c_src, 0);
          int b  = e_src.epsilon_3_index(c_src, 1);
          int c  = e_src.epsilon_3_index(c_src, 2);

          dcomplex src_factor = parity_src.value(beta) * Cg5_src.value(mu) * static_cast<double>(e_src.epsilon_3_value(c_src));
          
          u1_quark = sq1[a + Nc * id2_src[beta]].ptr(0) ;
          u2_quark = sq1[b + Nc * id1_src[mu]].ptr(0) ;
          d_quark  = sq1[c + Nc * id3_src[nu]].ptr(0) ;

          for (int t = 0; t < Nt; ++t) {
            dcomplex corr_t = cmplx(0.0,0.0) ;
            for (int ss = 0; ss < Nvol_s; ++ss) {
              int site = NCD2 * (ss + t * Nvol_s);
              
              for (int beta_p = 0; beta_p < Nd; ++beta_p){
                for (int mu_p = 0; mu_p < Nd; ++mu_p){
                  for (int nu_p = 0; nu_p < Nd; ++nu_p){
                    for (int c_snk = 0; c_snk < Nfaces; ++c_snk){

                      int a_p  = e_snk.epsilon_3_index(c_snk, 0);
                      int b_p  = e_snk.epsilon_3_index(c_snk, 1);
                      int c_p  = e_snk.epsilon_3_index(c_snk, 2);

                      dcomplex u1_prop, u2_prop, u3_prop, u4_prop, d_prop;

                      dcomplex snk_factor = parity_snk.value(beta_p) * Cg5_snk.value(mu_p) * static_cast<double>(e_snk.epsilon_3_value(c_snk));
                      
                      int ic1_r = 2 * a_p + id2_snk[beta_p] + site;
                      int ic1_i = ic1_r + 1 ;
                      u1_prop = cmplx(u1_quark[ic1_r], u1_quark[ic1_i]);
                      
                      int ic2_r = 2 * b_p + id1_snk[mu_p] + site;
                      int ic2_i = ic2_r + 1 ;
                      u2_prop = cmplx(u2_quark[ic2_r], u2_quark[ic2_i]);

                      u3_prop = cmplx(u1_quark[ic2_r], u1_quark[ic2_i]);
                      
                      u4_prop = cmplx(u2_quark[ic1_r], u2_quark[ic1_i]);
                      
                      int ic3_r = 2 * c_p + id3_snk[nu_p] + site;
                      int ic3_i = ic3_r + 1 ;
                      d_prop = cmplx(d_quark[ic3_r], d_quark[ic3_i]);

                      //corr_t += src_factor * snk_factor * (((u1_prop * u2_prop) - (u3_prop * u4_prop)) * d_prop);
                      corr_t += src_factor * snk_factor * u1_prop * u2_prop * d_prop;
                    }
                  }
                }
              }
            }
          corr_local[t] += corr_t;
          }
        }
      }
    }
  }
  global_corr_t(corr_global, corr_local);
}
//====================================================================
void Corr2pt_4spinor::nucleon_test(std::vector<dcomplex>& corr_global, // This version has parity projection at the source and sink. Agrees with the imported code. 
                                       const std::vector<Field_F>& sq1) // The proton correlator code which was native just provides a gm multiplication to the source indices
{                                                                       // and therefore not very useful.
  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);

  const int Nvol   = sq1[0].nvol();
  const int Nvol_s = Nvol / CommonParameters::Nt();

  const GammaMatrix unity       = m_gmset->get_GM(m_gmset->UNITY);
  const GammaMatrix gm5       = m_gmset->get_GM(m_gmset->GAMMA5);
  const GammaMatrix C         = m_gmset->get_GM(m_gmset->CHARGECONJG);
  const GammaMatrix Cg5_src  = C.mult(gm5);
  const GammaMatrix Cg5_snk  = C.mult(gm5);
  const GammaMatrix parity_src   = gm5.add(unity);
  const GammaMatrix parity_snk   = gm5.add(unity);
  //const GammaMatrix parity_src   = unity;
  //const GammaMatrix parity_snk   = unity;

  EpsilonTensor e_src, e_snk;

  std::vector<dcomplex> corr_local(Nt, cmplx(0.0, 0.0));

  int NC2 = 6 ;
  int NCD2 = 24;
  int Nfaces = 6;
  
  int id1_src[Nd],id2_src[Nd],id3_src[Nd], id1_snk[Nd],id2_snk[Nd],id3_snk[Nd];

  for (int id = 0; id < Nd; ++id) {
    id1_src[id] = parity_src.index(id) ;
    id3_src[id] = Cg5_src.index(id) ;
    id1_snk[id] = parity_snk.index(id) * NC2;
    id2_snk[id] = id * NC2;
    id3_snk[id] = Cg5_snk.index(id) * NC2;
  }
  
  const double *u1_quark, *d_quark, *u2_quark;

  for (int alpha = 0; alpha < Nd; ++alpha){
      for (int mu = 0; mu < Nd; ++mu){
          for (int c_src = 0; c_src < Nfaces; ++c_src){

            int a  = e_src.epsilon_3_index(c_src, 0);
            int b  = e_src.epsilon_3_index(c_src, 1);
            int c  = e_src.epsilon_3_index(c_src, 2);

            dcomplex src_factor = parity_src.value(alpha) * Cg5_src.value(mu) * static_cast<double>(e_src.epsilon_3_value(c_src));
            
            u1_quark = sq1[a + Nc * id1_src[alpha]].ptr(0) ;
            u2_quark = sq1[b + Nc * mu].ptr(0) ;
            d_quark  = sq1[c + Nc * id3_src[mu]].ptr(0) ;

            for (int t = 0; t < Nt; ++t) {
              dcomplex corr_t = cmplx(0.0,0.0) ;
              for (int ss = 0; ss < Nvol_s; ++ss) {
                int site = NCD2 * (ss + t * Nvol_s);
                 
                  for (int mu_p = 0; mu_p < Nd; ++mu_p){
                      for (int c_snk = 0; c_snk < Nfaces; ++c_snk){

                        int a_p  = e_snk.epsilon_3_index(c_snk, 0);
                        int b_p  = e_snk.epsilon_3_index(c_snk, 1);
                        int c_p  = e_snk.epsilon_3_index(c_snk, 2);

                        dcomplex u1_prop, u2_prop, u3_prop, u4_prop, d_prop;

                        dcomplex snk_factor = parity_snk.value(alpha) * Cg5_snk.value(mu_p) * static_cast<double>(e_snk.epsilon_3_value(c_snk));
                        
                        int ic1_r = 2 * a_p + id1_snk[alpha] + site;
                        int ic1_i = ic1_r + 1 ;
                        u1_prop = cmplx(u1_quark[ic1_r], u1_quark[ic1_i]);
                        
                        int ic2_r = 2 * b_p + id2_snk[mu_p] + site;
                        int ic2_i = ic2_r + 1 ;
                        u2_prop = cmplx(u2_quark[ic2_r], u2_quark[ic2_i]);

                        int ic3_r = 2 * c_p + id3_snk[mu_p] + site;
                        int ic3_i = ic3_r + 1 ;
                        d_prop = cmplx(d_quark[ic3_r], d_quark[ic3_i]);
                        
                        u3_prop = cmplx(u1_quark[ic2_r], u1_quark[ic2_i]);
                        
                        u4_prop = cmplx(u2_quark[ic1_r], u2_quark[ic1_i]);
                        
                        corr_t += src_factor * snk_factor * (((u1_prop * u2_prop) - (u3_prop * u4_prop)) * d_prop);
                        //corr_t +=  src_factor * snk_factor * u1_prop * u2_prop * d_prop;
                        }
                      }
                    }
                corr_local[t] += corr_t;
              }
            }
          }
        }
      
  global_corr_t(corr_global, corr_local);
}


void Corr2pt_4spinor::convert_gamma( const GammaMatrix gm, dcomplex gamma[4][4])
{
  const int Nd = CommonParameters::Nd();
  
  int id1[Nd];
  dcomplex val[Nd];

  for (int id = 0; id < Nd; ++id) {
    id1[id] = gm.index(id) ; // 
    val[id] = gm.value(id) ;
  }
  
  for (int alpha = 0; alpha < Nd; ++alpha ){
    for (int beta = 0; beta < Nd; ++beta ){
      gamma[alpha][beta] = cmplx(0.0,0.0); 
    }
  }

  for (int alpha = 0; alpha < Nd; ++alpha){
    gamma[alpha][id1[alpha]] = val[alpha];
  }

}

//====================================================================
void Corr2pt_4spinor::convert_prop( const std::vector<Field_F>& quark_prop, dcomplex *****propagator)
{
  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);

  const int Nvol   = quark_prop[0].nvol();
  const int Nvol_s = Nvol / CommonParameters::Nt();

  int NC2 = 6 ;
  int NCD2 = 24;
  
  int id1[Nd];

  for (int id = 0; id < Nd; ++id) {
    id1[id] = id * NC2;
  }
  
  const double *quark; 

  for (int alpha = 0; alpha < Nd; ++alpha){
    for (int a = 0; a < Nc; ++a){
      quark       = quark_prop[a + Nc * alpha].ptr(0) ;
      
        for (int ss = 0; ss < Nvol; ++ss) {
          int site = NCD2 * ss;
          
          for (int alpha_p = 0; alpha_p < Nd; ++alpha_p){
            for (int a_p = 0; a_p < Nc; ++a_p){
              int ic1_r = 2 * a_p + id1[alpha_p] + site;
              int ic1_i = ic1_r + 1 ;
              propagator[a][alpha][a_p][alpha_p][ss] = cmplx(quark[ic1_r],quark[ic1_i]);
          }
        }
      }
    }
  }
}

//====================================================================
void Corr2pt_4spinor::pion_test(std::vector<dcomplex>& corr_global,
                                       const GammaMatrix& gm_sink,
                                       const GammaMatrix& gm_src,
                                       const std::vector<Field_F>& sq1,
                                       const std::vector<Field_F>& sq2)
{
  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);

  const int Nvol   = sq1[0].nvol();
  const int Nvol_s = Nvol / CommonParameters::Nt();

  const GammaMatrix gm5         = m_gmset->get_GM(m_gmset->GAMMA5);
  const GammaMatrix gm_gm5_src  = gm_src.mult(gm5);
  const GammaMatrix gm5_gm_sink = gm5.mult(gm_sink);

  std::vector<dcomplex> corr_local(Nt, cmplx(0.0, 0.0));
  
  dcomplex *****quark_prop;

  quark_prop=(dcomplex*****)malloc(Nc*sizeof(dcomplex****));
  
  for(int a = 0 ; a < Nc ; ++a){ 
    quark_prop[a]=(dcomplex****)malloc(Nd*sizeof(dcomplex***));
    for(int alpha = 0; alpha < Nd; ++alpha){ 
      quark_prop[a][alpha]=(dcomplex***)malloc(Nc*sizeof(dcomplex**));
      for(int a_p = 0; a_p < Nc; ++a_p){  
        quark_prop[a][alpha][a_p]=(dcomplex**)malloc(Nd*sizeof(dcomplex*));
          for(int alpha_p = 0; alpha_p < Nd; ++alpha_p){   
             quark_prop[a][alpha][a_p][alpha_p]=(dcomplex*)malloc(Nvol*sizeof(dcomplex));
           }   
         }   
       }   
     }   
  
  convert_prop(sq1,quark_prop);
  
  dcomplex corr_t;

  for (int t = 0; t < Nt; ++t) {
    corr_t = 0 ;
    for (int ss = 0; ss < Nvol_s; ++ss) {
      int site = ss + t * Nvol_s;
      for(int a = 0; a < Nc; ++a){
        for(int alpha = 0; alpha < Nd; ++alpha){
          for(int a_p = 0; a_p < Nc; ++a_p){
            for(int alpha_p = 0; alpha_p < Nd; ++alpha_p){
            corr_t += quark_prop[a][alpha][a_p][alpha_p][site] * conj(quark_prop[a][alpha][a_p][alpha_p][site]);
            }
          }
        }
      }
    }
  corr_local[t] = +corr_t;
  }

  global_corr_t(corr_global, corr_local);

}

//====================================================================
dcomplex Corr2pt_4spinor::Sign(int a,int b,int c)
{  if(((a==0)&&(b==1)&&(c==2))|| ((a==1)&&(b==2)&&(c==0)) || ((a==2)&&(b==0)&&(c==1)))  { return ( cmplx(1.0,0.0) );} 
   else if (((a==0)&&(b==2)&&(c==1))|| ((a==1)&&(b==0)&&(c==2)) || ((a==2)&&(b==1)&&(c==0))) {return ( cmplx(-1.0,0.0) );}
   else {return (cmplx(0.0,0.0));}
}
//====================================================================

void Corr2pt_4spinor::nucleon_test1(std::vector<dcomplex>& corr_global,
                                       const std::vector<Field_F>& sq1)
{
  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);

  const int Nvol   = sq1[0].nvol();
  const int Nvol_s = Nvol / CommonParameters::Nt();

  const GammaMatrix unity       = m_gmset->get_GM(m_gmset->UNITY);
  const GammaMatrix gm5       = m_gmset->get_GM(m_gmset->GAMMA5);
  const GammaMatrix C         = m_gmset->get_GM(m_gmset->CHARGECONJG);
  const GammaMatrix Cg5_b  = C.mult(gm5);
  const GammaMatrix parity_b   = gm5.add(unity);
  const GammaMatrix polarisation_b   = unity;

  EpsilonTensor e_src, e_snk;
  
  dcomplex Cg5[4][4], parity[4][4], polarisation[4][4], spin_product, sign;

  convert_gamma(Cg5_b, Cg5);
  convert_gamma(parity_b, parity);
  //convert_gamma(unity, parity);
  convert_gamma(unity, polarisation);
  
  /*
  for (int alpha = 0; alpha < Nd; ++alpha ){
    for (int beta = 0; beta < Nd; ++beta ){
      vout.general(m_vl, "  %4d  %4d   %20.12e  %20.12e\n",
                 alpha, beta, real(polarisation[alpha][beta]), imag(polarisation[alpha][beta]));
      }
    }
  */
  std::vector<dcomplex> corr_local(Nt, cmplx(0.0, 0.0));
  
  dcomplex *****quark_prop;

  quark_prop=(dcomplex*****)malloc(Nc*sizeof(dcomplex****));
  
  for(int a = 0 ; a < Nc ; ++a){ 
    quark_prop[a]=(dcomplex****)malloc(Nd*sizeof(dcomplex***));
    for(int alpha = 0; alpha < Nd; ++alpha){ 
      quark_prop[a][alpha]=(dcomplex***)malloc(Nc*sizeof(dcomplex**));
      for(int a_p = 0; a_p < Nc; ++a_p){  
        quark_prop[a][alpha][a_p]=(dcomplex**)malloc(Nd*sizeof(dcomplex*));
          for(int alpha_p = 0; alpha_p < Nd; ++alpha_p){   
             quark_prop[a][alpha][a_p][alpha_p]=(dcomplex*)malloc(Nvol*sizeof(dcomplex));
           }   
         }   
       }   
     }   
  
  convert_prop(sq1,quark_prop);
  
  dcomplex corr_t;
  for(int alpha = 0; alpha < 4; alpha++){
    for(int beta = 0; beta < 4; beta++){
      for(int kappa = 0; kappa < 4; kappa++){
        for(int rho = 0; rho < 4; rho++){
          for(int alpha_p = 0; alpha_p < 4; alpha_p++){
            for(int beta_p = 0; beta_p < 4; beta_p++){
              for(int kappa_p = 0; kappa_p < 4; kappa_p++){
                for(int rho_p = 0; rho_p < 4; rho_p++){
                  spin_product = Cg5[kappa][rho] * Cg5[rho_p][kappa_p] * parity[alpha][beta] * parity[beta_p][alpha_p] * polarisation[alpha_p][alpha];
                  if( spin_product == cmplx(0.0, 0.0)){continue;}
                  for (int t = 0; t < Nt; ++t) {
                    corr_t = cmplx(0.0,0.0) ;
                    for (int ss = 0; ss < Nvol_s; ++ss) {
                      int site = ss + t * Nvol_s;
                      for(int c_src = 0; c_src < 6; c_src++){
                        int a  = e_src.epsilon_3_index(c_src, 0);
                        int b  = e_src.epsilon_3_index(c_src, 1);
                        int c  = e_src.epsilon_3_index(c_src, 2);
                        for(int c_snk = 0; c_snk < 6; c_snk++){
                          int a_p  = e_snk.epsilon_3_index(c_snk, 0);
                          int b_p  = e_snk.epsilon_3_index(c_snk, 1);
                          int c_p  = e_snk.epsilon_3_index(c_snk, 2);
                          sign = spin_product *   static_cast<double>(e_snk.epsilon_3_value(c_src)) * static_cast<double>(e_snk.epsilon_3_value(c_snk)); 
                          corr_t += sign * ((quark_prop[a][beta][a_p][beta_p][site] * quark_prop[b][kappa][b_p][kappa_p][site] * quark_prop[c][rho][c_p][rho_p][site] * cmplx(-1.0,0.0) ) + (quark_prop[a][beta][b_p][kappa_p][site] * quark_prop[b][kappa][a_p][beta_p][site] * quark_prop[c][rho][c_p][rho_p][site] * cmplx(1.0,0.0))) ;
                          //corr_t += sign * quark_prop[a][beta][a_p][beta_p][site] * quark_prop[b][kappa][b_p][kappa_p][site] * quark_prop[c][rho][c_p][rho_p][site] * cmplx(-1.0,0.0); 
                        }
                      }
                    }
                  corr_local[t] += corr_t;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  global_corr_t(corr_global, corr_local);

}


//====================================================================
//====================================================================

void Corr2pt_4spinor::nucleon_test2(std::vector<dcomplex>& corr_global,
                                       const std::vector<Field_F>& sq1)
{
  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);

  const int Nvol   = sq1[0].nvol();
  const int Nvol_s = Nvol / CommonParameters::Nt();

  const GammaMatrix unity       = m_gmset->get_GM(m_gmset->UNITY);
  const GammaMatrix gm5       = m_gmset->get_GM(m_gmset->GAMMA5);
  const GammaMatrix C         = m_gmset->get_GM(m_gmset->CHARGECONJG);
  const GammaMatrix Cg5_b  = C.mult(gm5);
  const GammaMatrix parity_b   = unity;
  const GammaMatrix polarisation_b   = unity;

  EpsilonTensor e_src, e_snk;
  
  dcomplex Cg5[4][4], parity[4][4], polarisation[4][4], spin_product, sign;

  convert_gamma(Cg5_b, Cg5);
  convert_gamma(unity, parity);
  convert_gamma(unity, polarisation);
  
  /*
  for (int alpha = 0; alpha < Nd; ++alpha ){
    for (int beta = 0; beta < Nd; ++beta ){
      vout.general(m_vl, "  %4d  %4d   %20.12e  %20.12e\n",
                 alpha, beta, real(polarisation[alpha][beta]), imag(polarisation[alpha][beta]));
      }
    }
  */
  std::vector<dcomplex> corr_local(Nt, cmplx(0.0, 0.0));
  
  dcomplex *****quark_prop;

  quark_prop=(dcomplex*****)malloc(Nc*sizeof(dcomplex****));
  
  for(int a = 0 ; a < Nc ; ++a){ 
    quark_prop[a]=(dcomplex****)malloc(Nd*sizeof(dcomplex***));
    for(int alpha = 0; alpha < Nd; ++alpha){ 
      quark_prop[a][alpha]=(dcomplex***)malloc(Nc*sizeof(dcomplex**));
      for(int a_p = 0; a_p < Nc; ++a_p){  
        quark_prop[a][alpha][a_p]=(dcomplex**)malloc(Nd*sizeof(dcomplex*));
          for(int alpha_p = 0; alpha_p < Nd; ++alpha_p){   
             quark_prop[a][alpha][a_p][alpha_p]=(dcomplex*)malloc(Nvol*sizeof(dcomplex));
           }   
         }   
       }   
     }   
  
  convert_prop(sq1,quark_prop);
  
  dcomplex corr_t;
    for(int beta = 0; beta < 4; beta++){
      for(int kappa = 0; kappa < 4; kappa++){
        for(int rho = 0; rho < 4; rho++){
            for(int beta_p = 0; beta_p < 4; beta_p++){
              for(int kappa_p = 0; kappa_p < 4; kappa_p++){
                for(int rho_p = 0; rho_p < 4; rho_p++){
                  spin_product = Cg5[kappa][rho] * Cg5[rho_p][kappa_p] * parity[beta_p][beta]; 
                  if( spin_product == cmplx(0.0, 0.0)){continue;}
                  for (int t = 0; t < Nt; ++t) {
                    corr_t = cmplx(0.0,0.0) ;
                    for (int ss = 0; ss < Nvol_s; ++ss) {
                      int site = ss + t * Nvol_s;
                      for(int c_src = 0; c_src < 6; c_src++){
                        int a  = e_src.epsilon_3_index(c_src, 0);
                        int b  = e_src.epsilon_3_index(c_src, 1);
                        int c  = e_src.epsilon_3_index(c_src, 2);
                        for(int c_snk = 0; c_snk < 6; c_snk++){
                          int a_p  = e_snk.epsilon_3_index(c_snk, 0);
                          int b_p  = e_snk.epsilon_3_index(c_snk, 1);
                          int c_p  = e_snk.epsilon_3_index(c_snk, 2);
                          sign = spin_product *   static_cast<double>(e_snk.epsilon_3_value(c_src)) * static_cast<double>(e_snk.epsilon_3_value(c_snk)); 
                          //corr_t += sign * ((quark_prop[a][beta][a_p][beta_p][site] * quark_prop[b][kappa][b_p][kappa_p][site] * quark_prop[c][rho][c_p][rho_p][site] * cmplx(-1.0,0.0) ) + (quark_prop[a][beta][b_p][kappa_p][site] * quark_prop[b][kappa][a_p][beta_p][site] * quark_prop[c][rho][c_p][rho_p][site] * cmplx(1.0,0.0))) ;
                          corr_t += sign * quark_prop[a][beta][a_p][beta_p][site] * quark_prop[b][kappa][b_p][kappa_p][site] * quark_prop[c][rho][c_p][rho_p][site]; 
                        }
                      }
                    }
                  corr_local[t] += corr_t;
                  }
                }
              }
            }
          }
        }
      }
  global_corr_t(corr_global, corr_local);

}


//====================================================================
//====================================================================

void Corr2pt_4spinor::get_antiquark(const std::vector<Field_F>& quark_prop, std::vector<Field_F>& anti_quark_prop)
{
  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);

  const int Nvol   = quark_prop[0].nvol();
  const int Nvol_s = Nvol / CommonParameters::Nt();

  const GammaMatrix unity       = m_gmset->get_GM(m_gmset->UNITY);
  const GammaMatrix gm5       = m_gmset->get_GM(m_gmset->GAMMA5);

  int NC2 = 6 ;
  int NCD2 = 24;
  
  int id1[Nd], id2[Nd], id3[Nd], id4[Nd], id5[Nd];

  for (int id = 0; id < Nd; ++id) {
    id1[id] = id * NC2;
    id2[id] = gm5.index(id) * NC2;
  }
  
  const double *quark; 
  double *anti_quark;

  for (int alpha = 0; alpha < Nd; ++alpha){
    for (int a = 0; a < Nc; ++a){
      quark       = quark_prop[a + Nc * id2[alpha]].ptr(0) ;
      anti_quark  = anti_quark_prop[a + Nc * id1[alpha]].ptr(0) ;
      
      for (int t = 0; t < Nt; ++t) {
        for (int ss = 0; ss < Nvol_s; ++ss) {
          int site = NCD2 * (ss + t * Nvol_s);
          
          for (int alpha_p = 0; alpha_p < Nd; ++alpha_p){
            for (int a_p = 0; a_p < Nc; ++a_p){
              dcomplex prop, prop1;

              int ic1_r = 2 * a_p + id2[alpha_p] + site;
              int ic1_i = ic1_r + 1 ;

              int ic2_r = 2 * a_p + id1[alpha_p] + site;
              int ic2_i = ic2_r + 1 ;
              
              prop = cmplx(quark[ic1_r],quark[ic1_i]);
              prop1 = gm5.value(alpha) *   conj(prop) * gm5.value(alpha_p);

              anti_quark[ic2_r] += real(prop1); 
              anti_quark[ic2_i] += imag(prop1);
            }
          }
        }
      }
    }
  }
}

void Corr2pt_4spinor::Tbb_Tbb(std::vector<dcomplex>& corr_global,
                                       const GammaMatrix& gm_snk,
                                       const GammaMatrix& gm_src,
                                       const double overall_factor,
                                       const std::vector<Field_F>& sq1,
                                       const std::vector<Field_F>& sq2)
{
  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);

  const int Nvol   = sq1[0].nvol();
  const int Nvol_s = Nvol / CommonParameters::Nt();

  const GammaMatrix gm5      = m_gmset->get_GM(m_gmset->GAMMA5);
  const GammaMatrix C        = m_gmset->get_GM(m_gmset->CHARGECONJG);
  const GammaMatrix Cg5_src  = C.mult(gm5);
  const GammaMatrix Cg5_snk  = C.mult(gm5);
  
  int NC2 = 6 ;
  int NCD2 = 24;
  
  int id1[Nd], id2[Nd], id3[Nd], id4[Nd], id5[Nd];

  for (int id = 0; id < Nd; ++id) {
    id1[id] = id * NC2;
    id2[id] = gm_src.index(id) * NC2;
    id3[id] = Cg5_src.index(id) * NC2;
    id4[id] = gm_snk.index(id) * NC2;
    id5[id] = Cg5_snk.index(id) * NC2;
  }
  
  const double *u_quark, *d_quark, *b1_quark, *b2_quark;

  std::vector<dcomplex> corr_local(Nt, cmplx(0.0, 0.0));
  
  dcomplex src_factor, snk_factor;

  for (int alpha = 0; alpha < Nd; ++alpha){
    for (int beta = 0; beta < Nd; ++beta){
      for (int rho = 0; rho < Nd; ++rho){
        for (int kappa = 0; kappa < Nd; ++kappa){
          for (int a = 0; a < Nc; ++a){
            for (int b = 0; b < Nc; ++b){
              u_quark  = sq1[a + Nc * id1[alpha]].ptr(0) ;
              d_quark  = sq1[b + Nc * id3[beta]].ptr(0) ;
              b1_quark = sq2[b + Nc * id2[rho]].ptr(0);
              b2_quark = sq2[a + Nc * id1[kappa]].ptr(0);
              
              src_factor = overall_factor * Cg5_src.value(alpha) * gm_src.value(kappa);
              
              for (int t = 0; t < Nt; ++t) {
                dcomplex corr_t = cmplx(0.0,0.0) ;
                
                for (int ss = 0; ss < Nvol_s; ++ss) {
                  int site = NCD2 * (ss + t * Nvol_s);
              
                  for (int alpha_p = 0; alpha_p < Nd; ++alpha_p){
                    for (int beta_p = 0; beta_p < Nd; ++beta_p){
                      for (int rho_p = 0; rho_p < Nd; ++rho_p){
                        for (int kappa_p = 0; kappa_p < Nd; ++kappa_p){
                          for (int a_p = 0; a_p < Nc; ++a_p){
                            for (int b_p = 0; b_p < Nc; ++b_p){
                              
                              dcomplex u_prop, d_prop, b1_prop, b2_prop, b3_prop, b4_prop;
                              
                              snk_factor = Cg5_snk.value(beta_p) * gm_snk.value(rho_p);
                              
                              int ic1_r = 2 * a_p + id5[alpha_p] + site;
                              int ic1_i = ic1_r + 1 ;
                              u_prop = cmplx(u_quark[ic1_r], u_quark[ic1_i]);

                              int ic2_r = 2 * b_p + id1[beta_p] + site;
                              int ic2_i = ic2_r + 1 ;
                              d_prop = cmplx(d_quark[ic2_r], d_quark[ic2_i]);
                              
                              int ic3_r = 2 * b_p + id1[rho_p] + site;
                              int ic3_i = ic3_r + 1 ;
                              b1_prop = cmplx(b1_quark[ic3_r], b1_quark[ic3_i]);
                              
                              int ic4_r = 2 * a_p + id4[kappa_p] + site;
                              int ic4_i = ic4_r + 1 ;
                              b2_prop = cmplx(b2_quark[ic4_r], b2_quark[ic4_i]);

                              b3_prop = cmplx(b1_quark[ic4_r], b1_quark[ic4_i]);
                              b4_prop = cmplx(b2_quark[ic3_r], b2_quark[ic3_i]);

                              corr_t += src_factor * snk_factor * u_prop * d_prop * ( (conj(b1_prop) * conj(b2_prop)) - (conj(b3_prop) * conj(b4_prop)) );

                            }
                          }
                        }
                      }
                    }
                  }
                }
              corr_local[t] += corr_t;
              }
            }
          }
        }
      }
    }
  }
  global_corr_t(corr_global, corr_local);
}

//====================================================================

void Corr2pt_4spinor::M_Tbb(std::vector<dcomplex>& corr_global,
                                       const GammaMatrix& gm_src,
                                       const GammaMatrix& gm_snk,
                                       const double overall_factor,
                                       const std::vector<Field_F>& sq1,
                                       const std::vector<Field_F>& sq2)
{
  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);

  const int Nvol   = sq1[0].nvol();
  const int Nvol_s = Nvol / CommonParameters::Nt();

  const GammaMatrix gm5_src = m_gmset->get_GM(m_gmset->GAMMA5);
  const GammaMatrix C       = m_gmset->get_GM(m_gmset->CHARGECONJG);
  const GammaMatrix Cg5_snk = C.mult(gm5_src);
  
  int NC2 = 6 ;
  int NCD2 = 24;
  
  int id1[Nd], id2[Nd], id3[Nd], id4[Nd], id5[Nd];

  for (int id = 0; id < Nd; ++id) {
    id1[id] = id * NC2;
    id2[id] = gm_src.index(id) * NC2; // this is gamma^i_{\kappa \rho}
    id3[id] = gm5_src.index(id) * NC2;
    id4[id] = gm_snk.index(id) * NC2; // this one (C gamma_i)_{\rho^\prime \kappa^\prime}
    id5[id] = Cg5_snk.index(id) * NC2;
  }
  
  std::vector<Field_F> anti_quark_prop(Nd*Nc);

  get_antiquark(sq2 , anti_quark_prop);

  const double *u_quark, *d_quark, *b1_quark, *b2_quark;

  std::vector<dcomplex> corr_local(Nt, cmplx(0.0, 0.0));
  
  dcomplex src_factor, snk_factor;

  for (int alpha = 0; alpha < Nd; ++alpha){
    for (int beta = 0; beta < Nd; ++beta){
      for (int rho = 0; rho < Nd; ++rho){
        for (int kappa = 0; kappa < Nd; ++kappa){
          for (int a = 0; a < Nc; ++a){
            for (int b = 0; b < Nc; ++b){
              d_quark  = sq1[b + Nc * id2[rho]].ptr(0) ;
              u_quark  = sq1[a + Nc * id3[beta]].ptr(0) ;
              
              src_factor = overall_factor * gm5_src.value(alpha) * gm_src.value(kappa);

              for (int t = 0; t < Nt; ++t) {
                dcomplex corr_t = cmplx(0.0,0.0) ;
                
                for (int ss = 0; ss < Nvol_s; ++ss) {
                  int site = NCD2 * (ss + t * Nvol_s);
              
                  for (int alpha_p = 0; alpha_p < Nd; ++alpha_p){
                    for (int beta_p = 0; beta_p < Nd; ++beta_p){
                      for (int rho_p = 0; rho_p < Nd; ++rho_p){
                        for (int kappa_p = 0; kappa_p < Nd; ++kappa_p){
                          for (int a_p = 0; a_p < Nc; ++a_p){
                            for (int b_p = 0; b_p < Nc; ++b_p){
                              
                              b1_quark = anti_quark_prop[b_p + Nc * id1[rho_p]].ptr(0);
                              b2_quark = anti_quark_prop[a_p + Nc * id4[kappa_p]].ptr(0);
                              
                              dcomplex u1_prop, u2_prop, d1_prop, d2_prop, b1_prop, b2_prop, b3_prop, b4_prop;
                              
                              snk_factor = Cg5_snk.value(beta_p) * gm_snk.value(rho_p);
                              
                              int ic1_r = 2 * a_p + id5[alpha_p] + site;
                              int ic1_i = ic1_r + 1 ;
                              u1_prop = cmplx(u_quark[ic1_r], u_quark[ic1_i]);

                              int ic2_r = 2 * b_p + id1[beta_p] + site;
                              int ic2_i = ic2_r + 1 ;
                              d1_prop = cmplx(d_quark[ic2_r], d_quark[ic2_i]);
                              
                              u2_prop = cmplx(u_quark[ic2_r], u_quark[ic2_i]);
                              d2_prop = cmplx(d_quark[ic1_r], d_quark[ic1_i]);
                              
                              int ic3_r = 2 * b + id1[kappa] + site;
                              int ic3_i = ic3_r + 1 ;
                              b1_prop = cmplx(b1_quark[ic3_r], b1_quark[ic3_i]);
                              
                              int ic4_r = 2 * a + id1[alpha] + site;
                              int ic4_i = ic4_r + 1 ;
                              b2_prop = cmplx(b2_quark[ic4_r], b2_quark[ic4_i]);

                              b3_prop = cmplx(b1_quark[ic4_r], b1_quark[ic4_i]);
                              b4_prop = cmplx(b2_quark[ic3_r], b2_quark[ic3_i]);

                              corr_t += src_factor * snk_factor * ( (u1_prop * d1_prop) +  (u2_prop * d2_prop) ) * ( (b1_prop * b2_prop) - (b3_prop * b4_prop) );

                            }
                          }
                        }
                      }
                    }
                  }
                }
              corr_local[t] += corr_t;
              }
            }
          }
        }
      }
    }
  }
  global_corr_t(corr_global, corr_local);
}

//====================================================================

void Corr2pt_4spinor::M_M(std::vector<dcomplex>& corr_global,
                                       const GammaMatrix& gm_src,
                                       const GammaMatrix& gm_snk,
                                       const double overall_factor,
                                       const std::vector<Field_F>& sq1,
                                       const std::vector<Field_F>& sq2)
{
  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);

  const int Nvol   = sq1[0].nvol();
  const int Nvol_s = Nvol / CommonParameters::Nt();

  const GammaMatrix gm5_src = m_gmset->get_GM(m_gmset->GAMMA5);
  const GammaMatrix gm5_snk = m_gmset->get_GM(m_gmset->GAMMA5);
  
  int NC2 = 6 ;
  int NCD2 = 24;
  
  int id1[Nd], id2[Nd], id3[Nd], id4[Nd], id5[Nd];

  for (int id = 0; id < Nd; ++id) {
    id1[id] = id * NC2;
    id2[id] = gm_src.index(id) * NC2; // 
    id3[id] = gm5_src.index(id) * NC2;
    id4[id] = gm_snk.index(id) * NC2; // 
    id5[id] = gm5_snk.index(id) * NC2;
  }
  
  std::vector<Field_F> anti_quark_prop(Nd*Nc);

  get_antiquark(sq2 , anti_quark_prop);

  const double *u_quark, *d_quark, *b1_quark, *b2_quark;

  std::vector<dcomplex> corr_local(Nt, cmplx(0.0, 0.0));
  
  dcomplex src_factor, snk_factor;

  for (int alpha = 0; alpha < Nd; ++alpha){
    for (int beta = 0; beta < Nd; ++beta){
      for (int rho = 0; rho < Nd; ++rho){
        for (int kappa = 0; kappa < Nd; ++kappa){
          for (int a = 0; a < Nc; ++a){
            for (int b = 0; b < Nc; ++b){
              d_quark  = sq1[b + Nc * id2[rho]].ptr(0) ;
              u_quark  = sq1[a + Nc * id3[beta]].ptr(0) ;
              
              src_factor = overall_factor * gm5_src.value(alpha) * gm_src.value(kappa);

              for (int t = 0; t < Nt; ++t) {
                dcomplex corr_t = cmplx(0.0,0.0) ;
                
                for (int ss = 0; ss < Nvol_s; ++ss) {
                  int site = NCD2 * (ss + t * Nvol_s);
              
                  for (int alpha_p = 0; alpha_p < Nd; ++alpha_p){
                    for (int beta_p = 0; beta_p < Nd; ++beta_p){
                      for (int rho_p = 0; rho_p < Nd; ++rho_p){
                        for (int kappa_p = 0; kappa_p < Nd; ++kappa_p){
                          for (int a_p = 0; a_p < Nc; ++a_p){
                            for (int b_p = 0; b_p < Nc; ++b_p){
                              
                              b1_quark = anti_quark_prop[b_p + Nc * id4[kappa_p]].ptr(0);
                              b2_quark = anti_quark_prop[a_p + Nc * id5[alpha_p]].ptr(0);
                              
                              dcomplex u1_prop, u2_prop, d1_prop, d2_prop, b1_prop, b2_prop, b3_prop, b4_prop;
                              
                              snk_factor = gm5_snk.value(beta_p) * gm_snk.value(rho_p);
                              
                              int ic1_r = 2 * a_p + id1[beta_p] + site;
                              int ic1_i = ic1_r + 1 ;
                              u1_prop = cmplx(u_quark[ic1_r], u_quark[ic1_i]);

                              int ic2_r = 2 * b_p + id1[rho_p] + site;
                              int ic2_i = ic2_r + 1 ;
                              d1_prop = cmplx(d_quark[ic2_r], d_quark[ic2_i]);
                              
                              u2_prop = cmplx(u_quark[ic2_r], u_quark[ic2_i]);
                              d2_prop = cmplx(d_quark[ic1_r], d_quark[ic1_i]);
                              
                              int ic3_r = 2 * b + id1[kappa] + site;
                              int ic3_i = ic3_r + 1 ;
                              b1_prop = cmplx(b1_quark[ic3_r], b1_quark[ic3_i]);
                              
                              int ic4_r = 2 * a + id1[alpha] + site;
                              int ic4_i = ic4_r + 1 ;
                              b2_prop = cmplx(b2_quark[ic4_r], b2_quark[ic4_i]);

                              b3_prop = cmplx(b1_quark[ic4_r], b1_quark[ic4_i]);
                              b4_prop = cmplx(b2_quark[ic3_r], b2_quark[ic3_i]);

                              corr_t += cmplx(2.0,0.0) * src_factor * snk_factor * ( (u1_prop * d1_prop) +  (u2_prop * d2_prop) ) * ( (b1_prop * b2_prop) - (b3_prop * b4_prop) );

                            }
                          }
                        }
                      }
                    }
                  }
                }
              corr_local[t] += corr_t;
              }
            }
          }
        }
      }
    }
  }
  global_corr_t(corr_global, corr_local);
}

//====================================================================

//====================================================================
double Corr2pt_4spinor::meson_all(const std::vector<Field_F>& sq1,
                                  const std::vector<Field_F>& sq2)
{
  const int Lt = CommonParameters::Lt();

  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  std::vector<dcomplex> corr(Lt);

  GammaMatrix gm_src  = m_gmset->get_GM(m_gmset->GAMMA5);
  GammaMatrix gm_sink = m_gmset->get_GM(m_gmset->GAMMA5);
  vout.general(m_vl, "PS <-- PS correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }
  const double result = real(corr[0]);

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA1);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA1);
  vout.general(m_vl, "V1 <-- V1 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA2);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA2);
  vout.general(m_vl, "V2 <-- V2 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA3);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA3);
  vout.general(m_vl, "V3 <-- V3 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA5);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA45);
  vout.general(m_vl, "A4 <-- PS correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA54);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA5);
  vout.general(m_vl, "PS <-- A4 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  gm_src  = m_gmset->get_GM(m_gmset->UNITY);
  gm_sink = m_gmset->get_GM(m_gmset->UNITY);
  vout.general(m_vl, "S <-- S correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA51);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA51);
  vout.general(m_vl, "GAMMA51 <-- GAMMA51 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA52);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA52);
  vout.general(m_vl, "GAMMA52 <-- GAMMA52 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA53);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA53);
  vout.general(m_vl, "GAMMA53 <-- GAMMA53 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  gm_src  = m_gmset->get_GM(m_gmset->SIGMA12);
  gm_sink = m_gmset->get_GM(m_gmset->SIGMA12);
  vout.general(m_vl, "SIGMA12 <-- SIGMA12 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  gm_src  = m_gmset->get_GM(m_gmset->SIGMA23);
  gm_sink = m_gmset->get_GM(m_gmset->SIGMA23);
  vout.general(m_vl, "SIGMA23 <-- SIGMA23 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  gm_src  = m_gmset->get_GM(m_gmset->SIGMA31);
  gm_sink = m_gmset->get_GM(m_gmset->SIGMA31);
  vout.general(m_vl, "SIGMA31 <-- SIGMA31 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }


  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }

  return result;
}


//====================================================================
void Corr2pt_4spinor::meson_correlator(std::vector<dcomplex>& corr_global,
                                       const GammaMatrix& gm_sink,
                                       const GammaMatrix& gm_src,
                                       const std::vector<Field_F>& sq1,
                                       const std::vector<Field_F>& sq2)
{
  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);

  const GammaMatrix gm5         = m_gmset->get_GM(m_gmset->GAMMA5);
  const GammaMatrix gm_gm5_src  = gm_src.mult(gm5);
  const GammaMatrix gm5_gm_sink = gm5.mult(gm_sink);

  std::vector<dcomplex> corr_local(Nt, cmplx(0.0, 0.0));
  for (int c0 = 0; c0 < Nc; ++c0) {
    for (int d0 = 0; d0 < Nd; ++d0) {
      int d1 = gm_gm5_src.index(d0);

      for (int t = 0; t < Nt; ++t) {
        dcomplex corr_t;

        contract_at_t(corr_t, gm5_gm_sink,
                      sq1[c0 + Nc * d0], sq2[c0 + Nc * d1], t);

        corr_local[t] += gm_gm5_src.value(d0) * corr_t;
      }
    }
  }
  global_corr_t(corr_global, corr_local);
}

//====================================================================
double Corr2pt_4spinor::meson_momentum_all(const std::vector<Field_F>& sq1,
                                           const std::vector<Field_F>& sq2,
                                           const std::vector<int>& source_position)
{
  const int Ndim = CommonParameters::Ndim();
  const int Lt   = CommonParameters::Lt();

  const int N_momentum = 10;

  typedef std::vector<int> MomentumSet;
  std::vector<MomentumSet> momentum_sink(N_momentum);
  for (int i_momentum = 0; i_momentum < N_momentum; i_momentum++) {
    momentum_sink[i_momentum].resize(Ndim - 1);
  }

  //- momentum_sink[0] = (1,0,0)
  int i_momentum = 0;
  momentum_sink[i_momentum][0] = 1;
  momentum_sink[i_momentum][1] = 0;
  momentum_sink[i_momentum][2] = 0;

  //- momentum_sink[1] = (0,1,0)
  i_momentum = 1;
  momentum_sink[i_momentum][0] = 0;
  momentum_sink[i_momentum][1] = 1;
  momentum_sink[i_momentum][2] = 0;

  //- momentum_sink[2] = (0,0,1)
  i_momentum = 2;
  momentum_sink[i_momentum][0] = 0;
  momentum_sink[i_momentum][1] = 0;
  momentum_sink[i_momentum][2] = 1;

  //- momentum_sink[3] = (1,1,0)
  i_momentum = 3;
  momentum_sink[i_momentum][0] = 1;
  momentum_sink[i_momentum][1] = 1;
  momentum_sink[i_momentum][2] = 0;

  //- momentum_sink[4] = (0,1,1)
  i_momentum = 4;
  momentum_sink[i_momentum][0] = 0;
  momentum_sink[i_momentum][1] = 1;
  momentum_sink[i_momentum][2] = 1;

  //- momentum_sink[5] = (1,0,1)
  i_momentum = 5;
  momentum_sink[i_momentum][0] = 1;
  momentum_sink[i_momentum][1] = 0;
  momentum_sink[i_momentum][2] = 1;

  //- momentum_sink[6] = (1,1,1)
  i_momentum = 6;
  momentum_sink[i_momentum][0] = 1;
  momentum_sink[i_momentum][1] = 1;
  momentum_sink[i_momentum][2] = 1;

  //- momentum_sink[7] = (2,0,0)
  i_momentum = 7;
  momentum_sink[i_momentum][0] = 2;
  momentum_sink[i_momentum][1] = 0;
  momentum_sink[i_momentum][2] = 0;

  //- momentum_sink[8] = (0,2,0)
  i_momentum = 8;
  momentum_sink[i_momentum][0] = 0;
  momentum_sink[i_momentum][1] = 2;
  momentum_sink[i_momentum][2] = 0;

  //- momentum_sink[9] = (0,0,2)
  i_momentum = 9;
  momentum_sink[i_momentum][0] = 0;
  momentum_sink[i_momentum][1] = 0;
  momentum_sink[i_momentum][2] = 2;


  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  std::vector<dcomplex> corr(Lt);

  GammaMatrix gm_src  = m_gmset->get_GM(m_gmset->GAMMA5);
  GammaMatrix gm_sink = m_gmset->get_GM(m_gmset->GAMMA5);
  for (int i_momentum = 0; i_momentum < N_momentum; i_momentum++) {
    vout.general(m_vl, "PS_momentum(%d %d %d) <-- PS correlator:\n",
                 momentum_sink[i_momentum][0],
                 momentum_sink[i_momentum][1],
                 momentum_sink[i_momentum][2]);
    meson_momentum_correlator(corr, momentum_sink[i_momentum], gm_sink, gm_src,
                              sq1, sq2, source_position);
    for (int t = 0; t < corr.size(); ++t) {
      vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                   t, real(corr[t]), imag(corr[t]));
    }
  }

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA1);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA1);
  for (int i_momentum = 0; i_momentum < N_momentum; i_momentum++) {
    vout.general(m_vl, "V1_momentum(%d %d %d) <-- V1 correlator:\n",
                 momentum_sink[i_momentum][0],
                 momentum_sink[i_momentum][1],
                 momentum_sink[i_momentum][2]);
    meson_momentum_correlator(corr, momentum_sink[i_momentum], gm_sink, gm_src,
                              sq1, sq2, source_position);
    for (int t = 0; t < corr.size(); ++t) {
      vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                   t, real(corr[t]), imag(corr[t]));
    }
  }

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA2);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA2);
  for (int i_momentum = 0; i_momentum < N_momentum; i_momentum++) {
    vout.general(m_vl, "V2_momentum(%d %d %d) <-- V2 correlator:\n",
                 momentum_sink[i_momentum][0],
                 momentum_sink[i_momentum][1],
                 momentum_sink[i_momentum][2]);
    meson_momentum_correlator(corr, momentum_sink[i_momentum], gm_sink, gm_src,
                              sq1, sq2, source_position);
    for (int t = 0; t < corr.size(); ++t) {
      vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                   t, real(corr[t]), imag(corr[t]));
    }
  }

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA3);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA3);
  for (int i_momentum = 0; i_momentum < N_momentum; i_momentum++) {
    vout.general(m_vl, "V3_momentum(%d %d %d) <-- V3 correlator:\n",
                 momentum_sink[i_momentum][0],
                 momentum_sink[i_momentum][1],
                 momentum_sink[i_momentum][2]);
    meson_momentum_correlator(corr, momentum_sink[i_momentum], gm_sink, gm_src,
                              sq1, sq2, source_position);
    for (int t = 0; t < corr.size(); ++t) {
      vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                   t, real(corr[t]), imag(corr[t]));
    }
  }


  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }

  return EXIT_SUCCESS;
}


//====================================================================
void Corr2pt_4spinor::meson_momentum_correlator(std::vector<dcomplex>& corr_global,
                                                const std::vector<int>& momentum_sink,
                                                const GammaMatrix& gm_sink,
                                                const GammaMatrix& gm_src,
                                                const std::vector<Field_F>& sq1,
                                                const std::vector<Field_F>& sq2,
                                                const std::vector<int>& source_position)
{
  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);

  const GammaMatrix gm5         = m_gmset->get_GM(m_gmset->GAMMA5);
  const GammaMatrix gm_gm5_src  = gm_src.mult(gm5);
  const GammaMatrix gm5_gm_sink = gm5.mult(gm_sink);

  std::vector<dcomplex> corr_local(Nt, cmplx(0.0, 0.0));
  for (int c0 = 0; c0 < Nc; ++c0) {
    for (int d0 = 0; d0 < Nd; ++d0) {
      int d1 = gm_gm5_src.index(d0);

      for (int t = 0; t < Nt; ++t) {
        dcomplex corr_t;

        contract_at_t(corr_t, momentum_sink, gm5_gm_sink, source_position,
                      sq1[c0 + Nc * d0], sq2[c0 + Nc * d1], t);

        corr_local[t] += gm_gm5_src.value(d0) * corr_t;
      }
    }
  }
  global_corr_t(corr_global, corr_local);
}


//====================================================================
double Corr2pt_4spinor::proton_test(const std::vector<Field_F>& sq_u,
                                    const std::vector<Field_F>& sq_d)
{
  const int Lt = CommonParameters::Lt();

  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }


  vout.general(m_vl, "proton <-- proton correlator(UNITY):\n");

  const GammaMatrix gm_unit   = m_gmset->get_GM(m_gmset->UNITY);
  const GammaMatrix gm_gamma0 = m_gmset->get_GM(m_gmset->GAMMA4);

  std::vector<dcomplex> p_corr_unity(Lt);
  proton_correlator(p_corr_unity, gm_unit, sq_u, sq_d);

  for (int t = 0; t < p_corr_unity.size(); t++) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(p_corr_unity[t]), imag(p_corr_unity[t]));
  }

  vout.general(m_vl, "proton <-- proton correlator(UPPER):\n");

  std::vector<dcomplex> p_corr_gamma0(Lt);
  proton_correlator(p_corr_gamma0, gm_gamma0, sq_u, sq_d);

  std::vector<dcomplex> p_corr_upper(Lt);
  for (int t = 0; t < p_corr_upper.size(); t++) {
    p_corr_upper[t] = (p_corr_unity[t] + p_corr_gamma0[t]) * 0.5;
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(p_corr_upper[t]), imag(p_corr_upper[t]));
  }

  vout.general(m_vl, "proton <-- proton correlator(GAMMA0):\n");

  for (int t = 0; t < p_corr_gamma0.size(); t++) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(p_corr_gamma0[t]), imag(p_corr_gamma0[t]));
  }


  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }

  const double result = real(p_corr_gamma0[0]);

  return result;
}


//====================================================================
void Corr2pt_4spinor::proton_correlator(std::vector<dcomplex>& corr_global,
                                        const GammaMatrix& gm,
                                        const std::vector<Field_F>& sq_u,
                                        const std::vector<Field_F>& sq_d)
{
  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(Nc == 3);
  assert(corr_global.size() == Lt);

  const GammaMatrix gm5 = m_gmset->get_GM(m_gmset->GAMMA5);
  const GammaMatrix c   = m_gmset->get_GM(m_gmset->CHARGECONJG);
  const GammaMatrix cg5 = c.mult(gm5);

#ifdef DEBUG
  vout.general(m_vl, "i:\tgm5\t\t\t\tc\t\t\t\tcg5\t\t\t\tgm\n");
  for (int i = 0; i < Nd; i++) {
    vout.general(m_vl, "%d:\t %d %e %e \t %d  %e %e \t %d  %e %e \t %d  %e %e \n",
                 i,
                 gm5.index(i), real(gm5.value(i)), imag(gm5.value(i)),
                 c.index(i), real(c.value(i)), imag(c.value(i)),
                 cg5.index(i), real(cg5.value(i)), imag(cg5.value(i)),
                 gm.index(i), real(gm.value(i)), imag(gm.value(i))
                 );
  }
#endif

  const int FactNc = 6;
  // This is valid only when Nc =3, which was already asserted.

  std::vector<dcomplex> corr_local(Nt);

  for (int t = 0; t < Nt; t++) {
    vout.paranoiac(m_vl, "# t= %d\n", t);

    dcomplex sum = cmplx(0.0, 0.0);
    for (int i_alpha = 0; i_alpha < Nd; i_alpha++) {
      int i_alphaP  = gm.index(i_alpha);
      int i_alpha3  = i_alpha;
      int i_alpha3P = i_alphaP;

      for (int i_alpha1P = 0; i_alpha1P < Nd; i_alpha1P++) {
        int i_alpha2P = cg5.index(i_alpha1P);

        for (int ic123P = 0; ic123P < FactNc; ic123P++) {
          EpsilonTensor epsilon_tensor;

          int      ic1P   = epsilon_tensor.epsilon_3_index(ic123P, 0);
          int      ic2P   = epsilon_tensor.epsilon_3_index(ic123P, 1);
          int      ic3P   = epsilon_tensor.epsilon_3_index(ic123P, 2);
          dcomplex factor = gm.value(i_alpha)
                            * cg5.value(i_alpha1P)
                            * static_cast<double>(epsilon_tensor.epsilon_3_value(ic123P));

          dcomplex sum1;
          contract_at_t(sum1, cg5, i_alpha3,
                        sq_u[ic1P + Nc * i_alpha1P],
                        sq_d[ic2P + Nc * i_alpha2P],
                        sq_u[ic3P + Nc * i_alpha3P], t);

          dcomplex sum2;
          contract_at_t(sum2, cg5, i_alpha3,
                        sq_u[ic3P + Nc * i_alpha3P],
                        sq_d[ic2P + Nc * i_alpha2P],
                        sq_u[ic1P + Nc * i_alpha1P], t);

          //sum += factor * (sum1 - sum2);
          sum += factor * (sum1 );
        }
      }
    }

    corr_local[t] = sum;
  } // t loop end.

  global_corr_t(corr_global, corr_local);
}


//====================================================================
// meson =tr(gm5.qn_sink.sq1.qn_src.gm5.(sq2)^\dagger);
void Corr2pt_4spinor::meson_correlator_x(std::vector<dcomplex>& meson,
                                         const GammaMatrix& qn_sink,
                                         const GammaMatrix& qn_src,
                                         const std::vector<Field_F>& sq1,
                                         const std::vector<Field_F>& sq2)
{
  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
  int Lx = CommonParameters::Lx();
  int Nx = CommonParameters::Nx();

  assert(meson.size() == Lx);

  GammaMatrix gm_src, gm_sink, gm5;
  gm5     = m_gmset->get_GM(m_gmset->GAMMA5);
  gm_src  = qn_src.mult(gm5);
  gm_sink = gm5.mult(qn_sink);

  std::vector<dcomplex> corr_local(Nx);
  for (int i = 0; i < Nx; ++i) {
    corr_local[i] = 0.0;
  }

  for (int c0 = 0; c0 < Nc; ++c0) {
    for (int d0 = 0; d0 < Nd; ++d0) {
      int d1 = gm_src.index(d0);

      for (int x = 0; x < Nx; ++x) {
        dcomplex corr_x;
        contract_at_x(corr_x, gm_sink,
                      sq1[c0 + Nc * d0], sq2[c0 + Nc * d1], x);

        corr_local[x] += gm_src.value(d0) * corr_x;
      }
    }
  }

  global_corr_x(meson, corr_local);
}


//====================================================================
void Corr2pt_4spinor::meson_momentum_correlator_x(std::vector<dcomplex>& corr_global,
                                                  const std::vector<int>& momentum_sink,
                                                  const GammaMatrix& gm_sink,
                                                  const GammaMatrix& gm_src,
                                                  const std::vector<Field_F>& sq1,
                                                  const std::vector<Field_F>& sq2,
                                                  const std::vector<int>& source_position)
{
  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
  int Lx = CommonParameters::Lx();
  int Nx = CommonParameters::Nx();

  assert(corr_global.size() == Lx);

  GammaMatrix gm_gm5_src, gm5_gm_sink, gm5;
  gm5         = m_gmset->get_GM(m_gmset->GAMMA5);
  gm_gm5_src  = gm_src.mult(gm5);
  gm5_gm_sink = gm5.mult(gm_sink);

  std::vector<dcomplex> corr_local(Nx);

  for (int c0 = 0; c0 < Nc; ++c0) {
    for (int d0 = 0; d0 < Nd; ++d0) {
      int d1 = gm_gm5_src.index(d0);

      for (int x = 0; x < Nx; ++x) {
        dcomplex corr_x;

        contract_at_x(corr_x, momentum_sink, gm5_gm_sink, source_position,
                      sq1[c0 + Nc * d0], sq2[c0 + Nc * d1], x);

        corr_local[x] += gm_gm5_src.value(d0) * corr_x;
      }
    }
  }

  global_corr_x(corr_global, corr_local);
}


//====================================================================
void Corr2pt_4spinor::proton_correlator_x(std::vector<dcomplex>& proton,
                                          const GammaMatrix& gm,
                                          const std::vector<Field_F>& squ,
                                          const std::vector<Field_F>& sqd)
{
  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
  int Lx = CommonParameters::Lx();
  int Nx = CommonParameters::Nx();

  assert(Nc == 3);
  assert(proton.size() == Lx);

  EpsilonTensor epsilon_tensor;

  GammaMatrix cg5, c, gm5;
  gm5 = m_gmset->get_GM(m_gmset->GAMMA5);
  c   = m_gmset->get_GM(m_gmset->CHARGECONJG);
  cg5 = c.mult(gm5);

  int FactNc = 6;
  // This is valid only when Nc =3, which was already asserted.

  std::vector<dcomplex> corr_local(Nx);

  for (int ix = 0; ix < Nx; ix++) {
    dcomplex sum = 0.0;
    dcomplex sum1, sum2;

    for (int ialph = 0; ialph < Nd; ialph++) {
      int ialphP  = gm.index(ialph);
      int ialph3  = ialph;
      int ialph3P = ialphP;

      for (int ialph1P = 0; ialph1P < Nd; ialph1P++) {
        int ialph2P = cg5.index(ialph1P);

        for (int ic123P = 0; ic123P < FactNc; ic123P++) {
          int      ic1P   = epsilon_tensor.epsilon_3_index(ic123P, 0);
          int      ic2P   = epsilon_tensor.epsilon_3_index(ic123P, 1);
          int      ic3P   = epsilon_tensor.epsilon_3_index(ic123P, 2);
          dcomplex factor = gm.value(ialph)
                            * cg5.value(ialph1P)
                            * static_cast<double>(epsilon_tensor.epsilon_3_value(ic123P));

          contract_at_x(sum1, cg5, ialph3,
                        squ[ic1P + Nc * ialph1P],
                        sqd[ic2P + Nc * ialph2P],
                        squ[ic3P + Nc * ialph3P], ix);
          contract_at_x(sum2, cg5, ialph3,
                        squ[ic3P + Nc * ialph3P],
                        sqd[ic2P + Nc * ialph2P],
                        squ[ic1P + Nc * ialph1P], ix);
          sum += factor * (sum1 - sum2);
        }
      }
    }

    corr_local[ix] = sum;
  } // it loop end.

  global_corr_x(proton, corr_local);
}


/* moved by tanigchi
//====================================================================
void Corr2pt_4spinor::global_corr_t(std::vector<dcomplex>& corr_global,
                                    const std::vector<dcomplex>& corr_local)
{
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);
  assert(corr_local.size() == Nt);

  const int ipe_t = Communicator::ipe(3);

  std::vector<dcomplex> corr_tmp(Lt, cmplx(0.0, 0.0));

  for (int t = 0; t < Nt; ++t) {
    int t_global = t + ipe_t * Nt;
    corr_tmp[t_global] = corr_local[t];
  }

  for (int t_global = 0; t_global < Lt; ++t_global) {
    double cr_r = Communicator::reduce_sum(real(corr_tmp[t_global]));
    double cr_i = Communicator::reduce_sum(imag(corr_tmp[t_global]));

    corr_global[t_global] = cmplx(cr_r, cr_i);
  }
}
*/

//====================================================================
//============================================================END=====
