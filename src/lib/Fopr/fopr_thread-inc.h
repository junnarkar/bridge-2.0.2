/*!
        @file    fopr_thread-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$
        @version $LastChangedRevision: 2492 $
*/

#ifndef FOPR_THREAD_INC_INCLUDED
#define FOPR_THREAD_INC_INCLUDED

#include "ResourceManager/threadManager.h"

namespace {
  // Equal tasks for all threads
  inline void set_threadtask(int& ith, int& nth, int& is, int& ns,
                             const int size)
  {
    nth = ThreadManager::get_num_threads();
    ith = ThreadManager::get_thread_id();

    size_t is2 = size_t(size) * size_t(ith) / nth;
    size_t ns2 = size_t(size) * size_t(ith + 1) / nth;
    is = int(is2);
    ns = int(ns2);
  }
}

#endif
//============================================================END=====
