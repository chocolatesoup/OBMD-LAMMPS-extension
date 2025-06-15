#!/bin/sh

mkdir -p build_kk

cd build_kk

cmake ../cmake \
-D PKG_DPD-BASIC=on \
-D PKG_MOLECULE=on \
-D PKG_OBMD=on \
-D PKG_RIGID=on \
-D PKG_EXTRA-DUMP=on \
-D PKG_EXTRA-MOLECULE=on \
-D BUILD_MPI=yes \
-D LAMMPS_MACHINE=mpi \
-D PKG_KOKKOS=on \
-D Kokkos_ARCH_PASCAL60=no \
-D Kokkos_ARCH_TURING75=yes \
-D Kokkos_ENABLE_CUDA=yes \
-D Kokkos_ENABLE_OPENMP=yes

make -j 12
