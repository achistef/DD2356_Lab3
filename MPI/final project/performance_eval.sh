#!/bin/bash -l
# The -l above is required to get the full environment with modules

# The name of the script is myjob
#SBATCH -J performance_eval

# Only 1 hour wall-clock time will be given to this job
#SBATCH -t 1:00:00
#SBATCH -A edu19.DD2356

# Number of nodes
#SBATCH --nodes=8

# Number of MPI processes per node
#SBATCH --ntasks-per-node=32

# Number of cores hosting OpenMP threads
#SBATCH -c 8

#SBATCH -e error_file.e

# TODO add threads
# export OMP_NUM_THREADS=8 


# Strong scalling, 10 GB

aprun -n 32 -N 32 ./bin/project.out -r 5 /cfs/klemming/scratch/s/sergiorg/DD2356/input/wikipedia_10GB.txt > results_32p_10GB.txt

aprun -n 64 -N 32 ./bin/project.out -r 5 /cfs/klemming/scratch/s/sergiorg/DD2356/input/wikipedia_10GB.txt > results_64p_10GB.txt

aprun -n 128 -N 32 ./bin/project.out -r 5 /cfs/klemming/scratch/s/sergiorg/DD2356/input/wikipedia_10GB.txt > results_128p_10GB.txt

aprun -n 256 -N 32 ./bin/project.out -r 5 /cfs/klemming/scratch/s/sergiorg/DD2356/input/wikipedia_10GB.txt > results_256p_10GB.txt


# Strong scalling, 160 GB

aprun -n 32 -N 32 ./bin/project.out -r 5 /cfs/klemming/scratch/s/sergiorg/DD2356/input/wikipedia_160GB.txt > results_32p_160GB.txt

aprun -n 64 -N 32 ./bin/project.out -r 5 /cfs/klemming/scratch/s/sergiorg/DD2356/input/wikipedia_160GB.txt > results_64p_160GB.txt

aprun -n 128 -N 32 ./bin/project.out -r 5 /cfs/klemming/scratch/s/sergiorg/DD2356/input/wikipedia_160GB.txt > results_128p_160GB.txt

aprun -n 256 -N 32 ./bin/project.out -r 5 /cfs/klemming/scratch/s/sergiorg/DD2356/input/wikipedia_160GB.txt > results_256p_160GB.txt