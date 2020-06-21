#!/bin/bash

#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -pe mpich 4
#$ -N MPI_Matrix10

MPICH_MACHINES=$TMPDIR/mpich_machines
cat $PE_HOSTFILE | awk '{print $1":1"}' > $MPICH_MACHINES

mpiexec -np 20 -hostfile $MPICH_MACHINES ./tspsec 10 tsp10.1

rm -rf $MPICH_MACHINES
