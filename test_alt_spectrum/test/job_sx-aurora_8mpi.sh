#!/bin/sh
#PBS -q debug
#PBS -T necmpi
#PBS -l elapstim_req=00:30:00
#PBS -o stdout_%r.%02j
#PBS -e stderr_%r.%02j
#PBS -venode=1
#PBS --venum-lhost=1

cd ${PBS_O_WORKDIR}

export VE_PROGINF=DETAIL

mpirun -np 8 -venode -nn 1 ./bridge.elf

