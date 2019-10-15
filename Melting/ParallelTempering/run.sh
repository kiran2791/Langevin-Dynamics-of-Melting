#!/bin/bash
#SBATCH -p mpi
#SBATCH -J T12Multi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=14-00
#SBATCH --output=simulation.out
#SBATCH --cpus-per-task=6

echo Runnning on `hostname`
echo Time is `date`
file='LD_woutput.cpp'
#g++ -O3 $file
mpic++ -O3 $file
number=$RANDOM
mpiexec -np 6 ./a.out $number > log.dat
echo Time is `date`
