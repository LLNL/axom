#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# SYS_TYPE: toss_3_x86_64_ib
# Compiler Spec: gcc@8.1.0
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/packages/cmake/cmake-3.14.5/bin/cmake
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------

set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-8.1.0/bin/gcc" CACHE PATH "")

set(CMAKE_CXX_COMPILER "/usr/tce/packages/gcc/gcc-8.1.0/bin/g++" CACHE PATH "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

set(CMAKE_Fortran_COMPILER "/usr/tce/packages/gcc/gcc-8.1.0/bin/gfortran" CACHE PATH "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

# Root directory for generated TPLs
set(TPL_ROOT "/usr/WS1/axom/libs/toss_3_x86_64_ib/2020_04_28_18_37_10/gcc-8.1.0" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.5.1" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.0" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.21" CACHE PATH "")

# SCR not built

set(RAJA_DIR "${TPL_ROOT}/raja-0.11.0" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2.1.0" CACHE PATH "")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3-gcc-8.1.0/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3-gcc-8.1.0/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3-gcc-8.1.0/bin/mpif90" CACHE PATH "")

set(MPIEXEC_EXECUTABLE "/usr/bin/srun" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

# Root directory for generated developer tools
set(DEVTOOLS_ROOT "/usr/WS1/axom/devtools/toss_3_x86_64_ib/2020_04_27_21_23_14/gcc-8.1.0" CACHE PATH "")

set(PYTHON_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.7.4/bin/python" CACHE PATH "")

set(ENABLE_DOCS ON CACHE BOOL "")

set(DOXYGEN_EXECUTABLE "${DEVTOOLS_ROOT}/doxygen-1.8.14/bin/doxygen" CACHE PATH "")

set(SPHINX_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.7.4/bin/sphinx-build" CACHE PATH "")

set(SHROUD_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.7.4/bin/shroud" CACHE PATH "")

set(UNCRUSTIFY_EXECUTABLE "${DEVTOOLS_ROOT}/uncrustify-0.61/bin/uncrustify" CACHE PATH "")

set(CPPCHECK_EXECUTABLE "${DEVTOOLS_ROOT}/cppcheck-1.87/bin/cppcheck" CACHE PATH "")

#------------------------------------------------------------------------------
# Other machine specifics
#------------------------------------------------------------------------------

set(ENABLE_OPENMP ON CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")


