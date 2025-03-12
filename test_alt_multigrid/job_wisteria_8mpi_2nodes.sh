#!/bin/bash
#PJM -L rscgrp=regular-o
#####PJM -L rscgrp=debug-o
#PJM -L node=2
#PJM -L elapse=0:30:00
#PJM -g wo22i060
#PJM --mpi proc=8
#PJM --omp thread=12
#PJM -j
#PJM -S

#export OMP_NUM_THREADS=12
export FLIB_BARRIER=HARD

mpiexec build/bridge_test_alt

