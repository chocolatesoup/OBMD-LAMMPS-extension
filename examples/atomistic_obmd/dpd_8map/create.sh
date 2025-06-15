#!/bin/bash

ncpu=24

#/net/gauss/home/petra/anaconda3/bin/python3 input.py 3.0

source="/temp/petra/programi/lmp_4mxsc/develop_pp-lmp_v1.0-april/pp_lmp-main/code/build_mpi/lmp_mpi"
mpirun --oversubscribe -np ${ncpu} ${source} -in in.simulation -l ./../out/data/log.lammps 
