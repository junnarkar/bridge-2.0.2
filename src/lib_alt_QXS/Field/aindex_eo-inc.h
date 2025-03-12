/*!
      @file    aindex_eo-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-02-28 16:09:41 #$
      @version $LastChangedRevision: 2492 $
*/

#include <assert.h>

#include "lib_alt_QXS/Field/aindex_eo.h"
#include "lib_alt_QXS/Field/afield.h"

//====================================================================
template<typename REALTYPE>
template<typename AFIELD>
void AIndex_eo<REALTYPE, QXS>::split(AFIELD& field_e,
                                     AFIELD& field_o,
                                     const AFIELD& field_lex)
{
  int Nin   = field_lex.nin();
  int Nex   = field_lex.nex();
  int Nvol  = field_lex.nvol();
  int Nvol2 = Nvol / 2;

  assert(field_e.check_size(Nin, Nvol2, Nex));
  assert(field_o.check_size(Nin, Nvol2, Nex));

  AIndex_lex<REALTYPE, QXS> index_lex;

  for (int ex = 0; ex < Nex; ++ex) {
    for (int ist = 0; ist < Nvol; ++ist) {
      int ix   = ist % Nx;
      int iyzt = ist / Nx;
      int ist2 = ist / 2;
      int ieo  = (ix + Leo[iyzt]) % 2;
      if (ieo == 0) {
        for (int in = 0; in < Nin; ++in) {
          int index1 = index_lex.idx(in, Nin, ist, ex);
          int index2 = idxh(in, Nin, ist2, ex);
          field_e.set(index2, field_lex.cmp(index1));
        }
      } else {
        for (int in = 0; in < Nin; ++in) {
          int index1 = index_lex.idx(in, Nin, ist, ex);
          int index2 = idxh(in, Nin, ist2, ex);
          field_o.set(index2, field_lex.cmp(index1));
        }
      }
    }
  }
}


//====================================================================
template<typename REALTYPE>
template<typename AFIELD>
void AIndex_eo<REALTYPE, QXS>::merge(AFIELD& field_lex,
                                     const AFIELD& field_e,
                                     const AFIELD& field_o)
{
  int Nin   = field_lex.nin();
  int Nex   = field_lex.nex();
  int Nvol  = field_lex.nvol();
  int Nvol2 = Nvol / 2;

  assert(field_e.check_size(Nin, Nvol2, Nex));
  assert(field_o.check_size(Nin, Nvol2, Nex));

  AIndex_lex<REALTYPE, QXS> index_lex;

  for (int ex = 0; ex < Nex; ++ex) {
    for (int ist = 0; ist < Nvol; ++ist) {
      int ix   = ist % Nx;
      int iyzt = ist / Nx;
      int ist2 = ist / 2;
      int ieo  = (ix + Leo[iyzt]) % 2;
      if (ieo == 0) {
        for (int in = 0; in < Nin; ++in) {
          int index1 = index_lex.idx(in, Nin, ist, ex);
          int index2 = idxh(in, Nin, ist2, ex);
          field_lex.set(index1, field_e.cmp(index2));
        }
      } else {
        for (int in = 0; in < Nin; ++in) {
          int index1 = index_lex.idx(in, Nin, ist, ex);
          int index2 = idxh(in, Nin, ist2, ex);
          field_lex.set(index1, field_o.cmp(index2));
        }
      }
    }
  }
}


//============================================================END=====
