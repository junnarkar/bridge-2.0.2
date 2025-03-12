#!/bin/bash
#PJM -L rscunit=fx
#PJM -L rscgrp=fx-debug
#PJM -L node=1
#PJM -L elapse=1:00:00
#PJM --mpi proc=4
#PJM -j

export OMP_NUM_THREADS=12
#export XOS_MMM_L_PAGING_POLICY=demand:demand:demand
export FLIB_BARRIER=HARD

mpiexec ./bridge.elf

