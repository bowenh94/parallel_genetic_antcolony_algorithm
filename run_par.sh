#!/bin/csh
#SBATCH --time=08:59:00
#SBATCH -N 4
#SBATCH -c 2

mpiexec ./par

echo 'done!'
