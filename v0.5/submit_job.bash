#!/bin/bash
#PBS -W group_list=partner712
#PBS -q debugq
#PBS -l walltime=00:09:10
#PBS -l mem=23gb
#PBS -l nodes=1:ppn=12

cd $PBS_O_WORKDIR

module load gsl/1.15
module load fftw-parallel/2.1.5

export OMP_THREADS_NUM=12
time ./halogen input/square_massofparts_withexclusion.input


