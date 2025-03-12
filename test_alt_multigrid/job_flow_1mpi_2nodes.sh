#!/bin/bash
#PJM -L rscunit=fx
#PJM -L rscgrp=fx-extra
#PJM -L node=2
#PJM -L elapse=0:05:00
#PJM --mpi proc=8
#PJM -j
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
echo "----------- test_alt_Coarse.yaml --------"
cat test_alt_Coarse.yaml
echo "----------- test_alt_Coarse.yaml, done --"
echo
BIN=./brdige_test_alt
ls -l ${BIN}
mpiexec -np 8 ${BIN}
date
