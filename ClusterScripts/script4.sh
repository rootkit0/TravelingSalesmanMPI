#!/bin/bash

#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -pe mpich 1
#$ -N MPI_Matrix4

MPICH_MACHINES=$TMPDIR/mpich_machines
cat $PE_HOSTFILE | awk '{print $1":1"}' > $MPICH_MACHINES

mpiexec -np 20 -hostfile $MPICH_MACHINES ./tspsec 4 tsp4.1

rm -rf $MPICH_MACHINES
