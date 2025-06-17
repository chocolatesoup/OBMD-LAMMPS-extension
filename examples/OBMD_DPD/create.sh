#!/bin/bash

ncpu=24

python3 input.py 

source="/temp/petra/workspace/lmp_4mxsc/develop/RELEASE_v0.1_D6.2_M30/OBMD-LAMMPS-extension-main/code/build_mpi/lmp_mpi"
mpirun --oversubscribe -np ${ncpu} ${source} -in in.simulation 


