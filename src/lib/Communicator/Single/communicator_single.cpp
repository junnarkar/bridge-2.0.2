/*!
        @file    communicator_single.cpp

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate: 2013-03-21 15:28:34 #$

        @version $LastChangedRevision: 2492 $
 */

#include "Communicator/communicator.h"

#include <cstdarg>
#include <cstring>
#include <cassert>

#include <sys/time.h>
#include <time.h>

#include "Parameters/commonParameters.h"

//====================================================================
int Communicator::init(int *pargc, char ***pargv)
{
  return EXIT_SUCCESS;
}


//====================================================================
int Communicator::finalize()
{
  return EXIT_SUCCESS;
}


//====================================================================
void Communicator::abort()
{
  ::abort();  // call system function.
}


//====================================================================
int Communicator::setup(int Ninstance)
{
  assert(Ninstance == 1);

  // check consistency
  if (!((CommonParameters::NPE() == 1) &&
        (CommonParameters::NPEx() == 1) &&
        (CommonParameters::NPEy() == 1) &&
        (CommonParameters::NPEz() == 1) &&
        (CommonParameters::NPEt() == 1)))
  {
    fprintf(stderr, "Error at Communicator::setup(): inappropriate grid_size.\n");
    exit(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}


//====================================================================
int Communicator::setup(
  const std::vector<int>& lattice_size,
  std::vector<int>& grid_size,
  int Ninstance)
{
  int Ndim = lattice_size.size();

  if (grid_size.size() != Ndim) {
    grid_size.resize(Ndim, 1);
  }

  int err = 0;
  for (int i = 0; i < Ndim; ++i) {
    if (grid_size[i] != 1) {
      ++err;
    }
  }

  if (err > 0) {
    printf("Error at Communicator: %s: unexpected grid_size.\n", __func__);
    exit(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}


//====================================================================
bool Communicator::is_primary()
{
  return true;
}


//====================================================================
bool Communicator::is_primary_master()
{
  return true;
}


//====================================================================
int Communicator::self()
{
  return 0;
}


//====================================================================
int Communicator::size()
{
  return 1;
}


//====================================================================
int Communicator::self_global()
{
  return 0;
}


//====================================================================
int Communicator::ipe(const int dir)
{
  return 0;
}


//====================================================================
int Communicator::npe(const int dir)
{
  return 1;
}


//====================================================================
int Communicator::grid_rank(int *rank, const int *g_coord)
{
  *rank = 0;
  return EXIT_SUCCESS;
}


//====================================================================
int Communicator::grid_coord(int *g_coord, const int rank)
{
  int Ndim = CommonParameters::Ndim();

  for (int i = 0; i < Ndim; ++i) {
    g_coord[i] = 0;
  }
  return EXIT_SUCCESS;
}


//====================================================================
int Communicator::grid_dims(int *g_dims)
{
  int Ndim = CommonParameters::Ndim();

  for (int i = 0; i < Ndim; ++i) {
    g_dims[i] = 1;
  }
  return EXIT_SUCCESS;
}


//====================================================================
int Communicator::sync()
{
  return EXIT_SUCCESS;
}


//====================================================================
// sync w/o (possible) busy wait
int Communicator::sync_usleep()
{
  return EXIT_SUCCESS;
}


//====================================================================
int Communicator::sync_global()
{
  return EXIT_SUCCESS;
}


//====================================================================
int Communicator::Base::broadcast(size_t size, void *data, int sender)
{
  // stay intact
  return EXIT_SUCCESS;
}


int Communicator::broadcast(int count, dcomplex *data, int sender)
{
  // stay intact
  return EXIT_SUCCESS;
}


int Communicator::broadcast(int count, double *data, int sender)
{
  // stay intact
  return EXIT_SUCCESS;
}


int Communicator::broadcast(int count, float *data, int sender)
{
  // stay intact
  return EXIT_SUCCESS;
}


int Communicator::broadcast(int count, int *data, int sender)
{
  // stay intact
  return EXIT_SUCCESS;
}


int Communicator::broadcast(int count, string& data, int sender)
{
  // stay intact
  return EXIT_SUCCESS;
}


//====================================================================
int Communicator::Base::exchange(size_t size, void *recv_buf, void *send_buf, int idir, int ipm, int tag)
{
  memcpy(recv_buf, send_buf, size);
  return EXIT_SUCCESS;
}


int Communicator::exchange(int count, dcomplex *recv_buf, dcomplex *send_buf, int idir, int ipm, int tag)
{
  memcpy(recv_buf, send_buf, sizeof(dcomplex) * count);
  return EXIT_SUCCESS;
}


int Communicator::exchange(int count, double *recv_buf, double *send_buf, int idir, int ipm, int tag)
{
  memcpy(recv_buf, send_buf, sizeof(double) * count);
  return EXIT_SUCCESS;
}


int Communicator::exchange(int count, float *recv_buf, float *send_buf, int idir, int ipm, int tag)
{
  memcpy(recv_buf, send_buf, sizeof(float) * count);
  return EXIT_SUCCESS;
}


int Communicator::exchange(int count, int *recv_buf, int *send_buf, int idir, int ipm, int tag)
{
  memcpy(recv_buf, send_buf, sizeof(int) * count);
  return EXIT_SUCCESS;
}


//====================================================================
int Communicator::Base::send_1to1(size_t size, void *recv_buf, void *send_buf, int send_to, int recv_from, int tag)
{
  memcpy(recv_buf, send_buf, size);
  return EXIT_SUCCESS;
}


int Communicator::send_1to1(int count, dcomplex *recv_buf, dcomplex *send_buf, int send_to, int recv_from, int tag)
{
  memcpy(recv_buf, send_buf, sizeof(dcomplex) * count);
  return EXIT_SUCCESS;
}


int Communicator::send_1to1(int count, double *recv_buf, double *send_buf, int send_to, int recv_from, int tag)
{
  memcpy(recv_buf, send_buf, sizeof(double) * count);
  return EXIT_SUCCESS;
}


int Communicator::send_1to1(int count, float *recv_buf, float *send_buf, int send_to, int recv_from, int tag)
{
  memcpy(recv_buf, send_buf, sizeof(float) * count);
  return EXIT_SUCCESS;
}


int Communicator::send_1to1(int count, int *recv_buf, int *send_buf, int send_to, int recv_from, int tag)
{
  memcpy(recv_buf, send_buf, sizeof(int) * count);
  return EXIT_SUCCESS;
}


//====================================================================
int Communicator::reduce_sum(int count, dcomplex *recv_buf, dcomplex *send_buf, int pattern)
{
  memcpy(recv_buf, send_buf, sizeof(dcomplex) * count);
  return EXIT_SUCCESS;
}


int Communicator::reduce_sum(int count, double *recv_buf, double *send_buf, int pattern)
{
  memcpy(recv_buf, send_buf, sizeof(double) * count);
  return EXIT_SUCCESS;
}


int Communicator::reduce_sum(int count, float *recv_buf, float *send_buf, int pattern)
{
  memcpy(recv_buf, send_buf, sizeof(float) * count);
  return EXIT_SUCCESS;
}


int Communicator::reduce_sum(int count, int *recv_buf, int *send_buf, int pattern)
{
  memcpy(recv_buf, send_buf, sizeof(int) * count);
  return EXIT_SUCCESS;
}


//====================================================================
int Communicator::reduce_max(int count, double *recv_buf, double *send_buf, int pattern)
{
  memcpy(recv_buf, send_buf, sizeof(double) * count);
  return EXIT_SUCCESS;
}


int Communicator::reduce_max(int count, float *recv_buf, float *send_buf, int pattern)
{
  memcpy(recv_buf, send_buf, sizeof(float) * count);
  return EXIT_SUCCESS;
}


int Communicator::reduce_max(int count, int *recv_buf, int *send_buf, int pattern)
{
  memcpy(recv_buf, send_buf, sizeof(int) * count);
  return EXIT_SUCCESS;
}


//====================================================================
dcomplex Communicator::reduce_sum(dcomplex v)
{
  return v;
}


//- NB. no reduce_{max,min} for dcomplex
// dcomplex Communicator::reduce_max(dcomplex v);
// dcomplex Communicator::reduce_min(dcomplex v);


//====================================================================
double Communicator::reduce_sum(double v)
{
  return v;
}


double Communicator::reduce_max(double v)
{
  return v;
}


double Communicator::reduce_min(double v)
{
  return v;
}


//====================================================================
float Communicator::reduce_sum(float v)
{
  return v;
}


float Communicator::reduce_max(float v)
{
  return v;
}


float Communicator::reduce_min(float v)
{
  return v;
}


//====================================================================
int Communicator::status()
{
#ifdef DEBUG
  printf("Communicator Single\n");
#endif
  return EXIT_SUCCESS;
}


//====================================================================
double Communicator::get_time()
{
  struct timeval now;

  if (gettimeofday(&now, (struct timezone *)0) != 0) {
    return double();
  }

  double sec = (double)now.tv_sec + ((double)now.tv_usec) * 1.0e-6;

  return sec;
}


//====================================================================
//============================================================END=====
