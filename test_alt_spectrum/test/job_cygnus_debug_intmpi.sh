#!/bin/bash

module load intmpi/${NQSV_MPI_VER}

cd ${PBS_O_WORKDIR}/

mpirun ${NQSII_MPIOPTS} -np 4 ${PBS_O_WORKDIR}/bridge.elf
