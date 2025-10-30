#!/bin/bash

/usr/bin/python3 input.py

ncpu=24
source="../../../build_mpi/lmp_mpi"
mpirun --oversubscribe -np ${ncpu} ${source} -in in.simulation -l log.lammps


