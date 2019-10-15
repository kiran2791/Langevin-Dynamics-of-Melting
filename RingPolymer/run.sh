#!/bin/bash
#SBATCH -p short
#SBATCH -J 700Lin_Eq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1-00
#SBATCH --output=job.out

file='LD_woutput.cpp'
g++ -O3 $file
number=$RANDOM
./a.out $number > log.dat

