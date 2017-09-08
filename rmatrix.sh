#!/bin/sh

#PBS -N Rmatrix
#PBS -l nodes=1:ppn=1
#PBS -l walltime=1:00:00

module load devel
cd $PBS_O_WORKDIR

./Rmatrix_calc < parameterrb87wL.inp > OutPut.data
