/*!
        @file    communicator.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#ifndef COMMUNICATOR_INCLUDED
#define COMMUNICATOR_INCLUDED

#include "configure.h"
#include "bridge_complex.h"
#include "bridge_defs.h"

#include <cstdio>
#include <cstdlib>
#include <cstddef>
#include <vector>
#include <string>
using std::string;

class Channel;

//! Communication library which wraps MPI.

/**
  This class provides a communication library which wraps
  MPI (Message Passing Interface) if the implementation
  file communicator_mpi.cpp is bound.
  In the single processor environment with no MPI library,
  communicator_single.cpp should be bound instead.
  [28 Dec 2011 H.Matsufuru]

  This class defines interface of inter-node communication routines.
  All methods are static, i.e. class-methods (like global).
  The explicit definitions are hidden in the implementation classes.

  Add complex args                [08 Aug 2020 Y.Namekawa]
*/

// forward declaration
//class Channel;

class Communicator {
 public:
  //! initialize communicator

  /**
     initialize communicator.
     @param pargc pointer to argc passed from main function.
     @param pargv pointer to argv passed from main function.

     they would further be passed to communication library
     for hints to machine-dependent configurations.
  */
  static int init(int *pargc, char ***pargv);

  //! finalize communicator

  /**
     finalize communicator.
     terminates communication environment.
  */
  static int finalize();

  //! terminate communicator

  /**
     terminate communicator immediately
     at some erroneous situations.
  */
  static void abort();

  //! setup communicator

  /**
     setup communicator environment such as logical layout.
     called after the parameters are obtained.
     @param ninstance specifies multiplicity of trivial parallelism.
     (expects 1 at present).
   */
  static int setup(int ninstance = 1);

  /**
     setup communicator environment with lattice size and
     (hint of) grid size provided.
     @param grid_size may be incomplete, in that case it is
     determined from the total number of processes and the
     lattice size and overwritten.
   */

  static int setup(const std::vector<int>& lattice_size,
                   std::vector<int>& grid_size,
                   int ninstance = 1);

// info about rank
  static bool is_primary();        //!< check if the present node is primary in small communicator.
  static bool is_primary_master(); //!< check if the present node is primary in global communicator.

  static int self();               //!< rank within small world.

  static int nodeid() { return self(); }  //!< alternative name for self().
  static int size();   //!< size of small world.

#ifdef ENABLE_MULTI_INSTANCE
  static int self_global(); //!< rank within global communicator.
  static int world_id();    //!< id for the current small world.
#endif

// layout
  static int ipe(const int dir);                          //!< logical coordinate of current proc.
  static int npe(const int dir);                          //!< logical grid extent

  static int grid_rank(int *rank, const int *grid_coord); //!< find rank number from grid coordinate.
  static int grid_coord(int *grid_coord, const int rank); //!< find grid coordinate from rank number.
  static int grid_dims(int *grid_dims);                   //!< find grid dimensions.

// synchronize
  static int sync();   //!< synchronize within small world.

// synchronize (w/o busy wait)
  static int sync_usleep();   //!< synchronize within small world. (slow but no busy wait)

#ifdef ENABLE_MULTI_INSTANCE
  static int sync_global();   //!< synchronize all processes.
#endif

// data transfer
  static int broadcast(int count, dcomplex *data, int sender);                                            //!< broadcast array of dcomplex from sender.
  static int broadcast(int count, double *data, int sender);                                              //!< broadcast array of double from sender.
  static int broadcast(int count, float *data, int sender);                                               //!< broadcast array of float from sender.
  static int broadcast(int count, int *data, int sender);                                                 //!< broadcast array of integer from sender.
  static int broadcast(int count, string& data, int sender);                                              //!< broadcast a string from sender. count is insignificant.

  static int exchange(int count, dcomplex *recv_buf, dcomplex *send_buf, int idir, int ipm, int tag);     //!< receive array of dcomplex from upstream specified by idir and ipm, and send array to downstream.
  static int exchange(int count, double *recv_buf, double *send_buf, int idir, int ipm, int tag);         //!< receive array of double from upstream specified by idir and ipm, and send array to downstream.
  static int exchange(int count, float *recv_buf, float *send_buf, int idir, int ipm, int tag);           //!< receive array of float from upstream specified by idir and ipm, and send array to downstream.
  static int exchange(int count, int *recv_buf, int *send_buf, int idir, int ipm, int tag);               //!< receive array of int from upstream specified by idir and ipm, and send array to downstream.

  static int send_1to1(int count, dcomplex *recv_buf, dcomplex *send_buf, int p_to, int p_from, int tag); //!< send array of dcomplex from rank p_from to rank p_to. communication distinguished by tag.
  static int send_1to1(int count, double *recv_buf, double *send_buf, int p_to, int p_from, int tag);     //!< send array of double from rank p_from to rank p_to. communication distinguished by tag.
  static int send_1to1(int count, float *recv_buf, float *send_buf, int p_to, int p_from, int tag);       //!< send array of float from rank p_from to rank p_to. communication distinguished by tag.
  static int send_1to1(int count, int *recv_buf, int *send_buf, int p_to, int p_from, int tag);           //!< send array of int from rank p_from to rank p_to. communication distinguished by tag.

  static int reduce_sum(int count, dcomplex *recv_buf, dcomplex *send_buf, int pattern = 0);              //!< make a global sum of an array of dcomplex over the communicator. pattern specifies the dimensions to be reduced.
  static int reduce_sum(int count, double *recv_buf, double *send_buf, int pattern     = 0);              //!< make a global sum of an array of double over the communicator. pattern specifies the dimensions to be reduced.
  static int reduce_sum(int count, float *recv_buf, float *send_buf, int pattern       = 0);              //!< make a global sum of an array of float over the communicator. pattern specifies the dimensions to be reduced.
  static int reduce_sum(int count, int *recv_buf, int *send_buf, int pattern           = 0);              //!< make a global sum of an array of int over the communicator. pattern specifies the dimensions to be reduced.

  static int reduce_max(int count, double *recv_buf, double *send_buf, int pattern = 0);                  //!< find a global maximum of an array of double over the communicator. pattern specifies the dimensions to be reduced.
  static int reduce_max(int count, float *recv_buf, float *send_buf, int pattern   = 0);                  //!< find a global maximum of an array of float over the communicator. pattern specifies the dimensions to be reduced.
  static int reduce_max(int count, int *recv_buf, int *send_buf, int pattern       = 0);                  //!< find a global maximum of an array of int over the communicator. pattern specifies the dimensions to be reduced.

  static int reduce_min(int count, double *recv_buf, double *send_buf, int pattern = 0);                  //!< find a global minimum of an array of double over the communicator. pattern specifies the dimensions to be reduced.
  static int reduce_min(int count, float *recv_buf, float *send_buf, int pattern   = 0);                  //!< find a global minimum of an array of float over the communicator. pattern specifies the dimensions to be reduced.
  static int reduce_min(int count, int *recv_buf, int *send_buf, int pattern       = 0);                  //!< find a global minimum of an array of int over the communicator. pattern specifies the dimensions to be reduced.

  static dcomplex reduce_sum(dcomplex);                                                                   //!< alternative interface to reduce_sum(). returns the global sum of a dcomplex over the whole communicator.
  //- NB. no reduce_{max,min} for dcomplex
  // static dcomplex reduce_max(dcomplex);
  // static dcomplex reduce_min(dcomplex);

  static double reduce_sum(double);                                                                   //!< alternative interface to reduce_sum(). returns the global sum of a double over the whole communicator.
  static double reduce_max(double);                                                                   //!< alternative interface to reduce_max(). returns the global sum of a double over the whole communicator.
  static double reduce_min(double);                                                                   //!< alternative interface to reduce_min(). returns the global sum of a double over the whole communicator.

  static float reduce_sum(float);                                                                     //!< alternative interface to reduce_sum(). returns the global sum of a float over the whole communicator.
  static float reduce_max(float);                                                                     //!< alternative interface to reduce_max(). returns the global sum of a float over the whole communicator.
  static float reduce_min(float);                                                                     //!< alternative interface to reduce_min(). returns the global sum of a float over the whole communicator.

  static double get_time();                                                                           //!< obtain a wall-clock time.

  // async communication
  static Channel *send_init(int count, int idir, int ipm);
  static Channel *recv_init(int count, int idir, int ipm);

  // async communication with given buffer [2017.09.02 H.Matsufuru]
  static Channel *send_init(int count, int idir, int ipm, void *buf);
  static Channel *recv_init(int count, int idir, int ipm, void *buf);

// debug
  static int status();

//! base case

  /**
     this class defines base case for data exchange that are specified
     by plain streams of bytes of size.
   */
  class Base {
   public:
    static int broadcast(size_t size, void *data, int sender);
    static int exchange(size_t size, void *recv_buf, void *send_buf, int idir, int ipm, int tag);

    static int send_1to1(size_t size, void *recv_buf, void *send_buf, int send_to, int recv_from, int tag);
  };

 private:

//! no instance at all

/**
   communicator class is not indented to be instantiated.
   constructor, copy constroctor, assignment operator, and destructor
   are defined as private.
*/
  Communicator() {}
  Communicator(const Communicator&) {}
  Communicator& operator=(const Communicator&);

  ~Communicator() {}
};
#endif /* COMMUNICATOR_INCLUDED */
