/*!
        @file    tensorProd.cpp

        @brief

        @author
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2013-03-21 15:28:34 #$

        @version $LastChangedRevision: 1928 $
*/

#include "Force/Fermion/tensorProd.h"
#include <cassert>

//====================================================================
void tensorProd_Field_F(Field_G& u, const Field_F& v1, const Field_F& v2)
{
  return tensorProd_Field_F(u, 0, v1, v2);
}


//====================================================================
void tensorProd_Field_F(Field_G& u, const int ex, const Field_F& v1, const Field_F& v2)
{
  const int Nvol = u.nvol();
  const int Nc   = CommonParameters::Nc();
  const int Nd   = CommonParameters::Nd();

  assert(u.nvol() == v1.nvol());
  assert(u.nvol() == v2.nvol());
  assert(ex < u.nex());
  assert(v1.nex() == 1);
  assert(v2.nex() == 1);

  Mat_SU_N ut(v1.nc());
  double   ut_r, ut_i;

  for (int site = 0; site < Nvol; ++site) {
    for (int c1 = 0; c1 < Nc; ++c1) {
      for (int c2 = 0; c2 < Nc; ++c2) {
        ut_r = 0.0;
        ut_i = 0.0;
        for (int s = 0; s < Nd; ++s) {
          ut_r += v1.cmp_r(c2, s, site) * v2.cmp_r(c1, s, site)
                  + v1.cmp_i(c2, s, site) * v2.cmp_i(c1, s, site);
          ut_i += v1.cmp_r(c2, s, site) * v2.cmp_i(c1, s, site)
                  - v1.cmp_i(c2, s, site) * v2.cmp_r(c1, s, site);
        }
        ut.set(c1, c2, ut_r, ut_i);
      }
    }
    u.set_mat(site, ex, ut);
  }
}


//====================================================================
//============================================================END=====
