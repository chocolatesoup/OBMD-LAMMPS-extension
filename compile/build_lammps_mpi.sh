#!/bin/sh

mkdir -p build_mpi

cd build_mpi

cmake ../cmake \
-D PKG_DPD-BASIC=on \
-D PKG_MOLECULE=on \
-D PKG_OBMD=on \
-D PKG_RIGID=on \
-D PKG_EXTRA-DUMP=on \
-D PKG_EXTRA-COMPUTE=on \
-D PKG_EXTRA-MOLECULE=on \
-D BUILD_MPI=yes \
-D CMAKE_BUILD_TYPE=Debug \
-D CMAKE_EXPORT_COMPILE_COMMANDS=on \
-D LAMMPS_MACHINE=mpi

make -j 12
