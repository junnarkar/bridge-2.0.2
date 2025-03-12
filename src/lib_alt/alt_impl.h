/*!
        @file    alt_impl.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$
        @version $LastChangedRevision: 2492 $
*/

#ifndef ALT_IMPL_INCLUDED
#define ALT_IMPL_INCLUDED

enum Impl
{
  CORELIB, SIMD, SIMD2, VECTOR, QXS, ACCEL, OPENACC
};

// alignment for each IMPL
template<Impl IMPL>
constexpr int alignment_size();

#endif
