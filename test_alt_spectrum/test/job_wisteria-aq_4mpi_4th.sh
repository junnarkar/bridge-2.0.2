#!/bin/bash
#PJM -L rscgrp=share-debug
#PJM -L gpu=4
#PJM -L elapse=0:30:00
#PJM --omp thread=4
####PJM -g wo22i060
#PJM -g ga42
#PJM -j

module load nvidia
module load ompi

mpiexec ${NQSII_MPIOPTS} -np 4 -npernode 4 -bind-to none ./bridge.elf

