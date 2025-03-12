#!/bin/bash
#PJM -L rscgrp=regular-o
#####PJM -L rscgrp=debug-o
#PJM -L node=1
#PJM -L elapse=0:30:00
#PJM -g wo22i060
#PJM --mpi proc=4
#PJM --omp thread=12
#PJM -j

#export OMP_NUM_THREADS=12
export FLIB_BARRIER=HARD

mpiexec ./bridge.elf

