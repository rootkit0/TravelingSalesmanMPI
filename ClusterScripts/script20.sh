#!/bin/bash

#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -pe mpich 4
#$ -N MPI_Matrix20

MPICH_MACHINES=$TMPDIR/mpich_machines
cat $PE_HOSTFILE | awk '{print $1":1"}' > $MPICH_MACHINES

mpiexec -np 20 -hostfile $MPICH_MACHINES ./tspsec 20 tsp20.1

rm -rf $MPICH_MACHINES
