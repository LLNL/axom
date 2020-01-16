#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# SYS_TYPE: blueos_3_ppc64le_ib_p9
# Compiler Spec: gcc@7.3.1
#------------------------------------------------------------------------------

# CMake executable path: /usr/WS1/axom/thirdparty_libs/builds/2019_11_21_10_32_05/gcc-7.3.1/cmake-3.9.6/bin/cmake

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------

# C compiler used by spack
set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-7.3.1/bin/gcc" CACHE PATH "")

# C++ compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/tce/packages/gcc/gcc-7.3.1/bin/g++" CACHE PATH "")

# Fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE BOOL "")

set(CMAKE_Fortran_COMPILER "/usr/tce/packages/gcc/gcc-7.3.1/bin/gfortran" CACHE PATH "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

# Root directory for generated TPLs
set(TPL_ROOT "/usr/WS1/axom/thirdparty_libs/builds/2019_11_21_10_32_05/gcc-7.3.1" CACHE PATH "")

# conduit from uberenv
set(CONDUIT_DIR "${TPL_ROOT}/conduit-master" CACHE PATH "")

# mfem from uberenv
set(MFEM_DIR "${TPL_ROOT}/mfem-4.0" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.19" CACHE PATH "")

# scr not built by uberenv

# raja from uberenv
set(RAJA_DIR "${TPL_ROOT}/raja-0.9.0/share/raja/cmake" CACHE PATH "")

# umpire from uberenv
set(UMPIRE_DIR "${TPL_ROOT}/umpire-1.0.0/share/umpire/cmake" CACHE PATH "")

# python not built by uberenv

set(ENABLE_DOCS OFF CACHE BOOL "")

# shroud not built by uberenv

# uncrustify not built by uberenv

# lcov and genhtml not built by uberenv

# cppcheck not built by uberenv

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_C_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-gcc-7.3.1/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-gcc-7.3.1/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-gcc-7.3.1/bin/mpif90" CACHE PATH "")

set(MPIEXEC "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-gcc-7.3.1/bin/mpirun" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-np" CACHE PATH "")

set(BLT_MPI_COMMAND_APPEND "mpibind" CACHE PATH "")

#------------------------------------------------------------------------------
# Other machine specifics
#------------------------------------------------------------------------------

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

set(ENABLE_OPENMP ON CACHE BOOL "")


