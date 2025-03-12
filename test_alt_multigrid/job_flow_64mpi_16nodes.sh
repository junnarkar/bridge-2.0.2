#!/bin/bash
#PJM -L rscunit=fx
#PJM -L rscgrp=fx-debug
#PJM -L node=16
#PJM -L elapse=0:20:00
#PJM --mpi proc=64
#PJM -j
#PJM -S
##PJM -L rscgrp=fx-debug

export OMP_NUM_THREADS=12

date
echo "----------- this script -----------------"
cat $0
echo "----------- this script, done -----------"
echo
echo "----------- main.yaml -------------------"
cat main.yaml
echo "----------- main.yaml, done -------------"
echo
echo "----------- test_alt_Multigrid_32x64.yaml --------"
cat test_alt_Multigrid_32x64.yaml
echo "----------- test_alt_Multigrid_32x64.yaml --"
echo
BIN=./build/bridge_test_alt
ls -l ${BIN}
mpiexec ${BIN}
date
