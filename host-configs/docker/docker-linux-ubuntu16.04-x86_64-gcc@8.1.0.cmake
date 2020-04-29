##################################
# !!!! This is a generated file, edit at own risk !!!!
##################################

# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
##################################

##################################

# SYS_TYPE: linux-ubuntu16.04-x86_64
# Compiler Spec: gcc@8.1.0
##################################

# CMake executable path: /home/axom/axom_tpls/gcc-8.1.0/cmake-3.10.1/bin/cmake

##############
# Compilers
##############

# C compiler used by spack
set(CMAKE_C_COMPILER "/usr/bin/gcc" CACHE PATH "")

# C++ compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/bin/g++" CACHE PATH "")

# Fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE BOOL "")

set(CMAKE_Fortran_COMPILER "/usr/bin/gfortran" CACHE PATH "")

##############
# TPLs
##############

# Root directory for generated TPLs
set(TPL_ROOT "/home/axom/axom_tpls/gcc-8.1.0" CACHE PATH "")

# conduit from uberenv
set(CONDUIT_DIR "${TPL_ROOT}/conduit-master" CACHE PATH "")

# mfem from uberenv
set(MFEM_DIR "${TPL_ROOT}/mfem-4.0" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.21" CACHE PATH "")

# scr not built by uberenv

# raja from uberenv
set(RAJA_DIR "/usr/local/share/raja" CACHE PATH "")

# umpire from uberenv
set(UMPIRE_DIR "${TPL_ROOT}/umpire-2.1.0" CACHE PATH "")

# python not built by uberenv

set(ENABLE_DOCS OFF CACHE BOOL "")

# shroud not built by uberenv

# uncrustify not built by uberenv

# lcov and genhtml not built by uberenv

# cppcheck not built by uberenv

##############
# MPI
##############

set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_C_COMPILER "/home/axom/axom_tpls/gcc-8.1.0/mpich-3.3.2/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/home/axom/axom_tpls/gcc-8.1.0/mpich-3.3.2/bin/mpic++" CACHE PATH "")

set(MPI_Fortran_COMPILER "/home/axom/axom_tpls/gcc-8.1.0/mpich-3.3.2/bin/mpif90" CACHE PATH "")

set(MPIEXEC "/home/axom/axom_tpls/gcc-8.1.0/mpich-3.3.2/bin/mpiexec" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

##############
# Other machine specifics
##############

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

set(ENABLE_OPENMP ON CACHE BOOL "")


