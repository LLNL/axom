##################################
# !!!! This is a generated file, edit at own risk !!!!
##################################

##################################
#
# Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory.
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
#
##################################

##################################

# SYS_TYPE: blueos_3_ppc64le_ib
# Compiler Spec: clang@coral_xlf
##################################

# CMake executable path: /usr/workspace/wsrzc/axom/thirdparty_libs/builds/2018_11_14_15_17_37/clang-coral_xlf/cmake-3.9.6/bin/cmake

##############
# Compilers
##############

# C compiler used by spack
set(CMAKE_C_COMPILER "/usr/tce/packages/clang/clang-coral-2018.05.23/bin/clang" CACHE PATH "")

# C++ compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/tce/packages/clang/clang-coral-2018.05.23/bin/clang++" CACHE PATH "")

# Fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE BOOL "")

set(CMAKE_Fortran_COMPILER "/usr/tce/packages/xl/xl-2018.05.18/bin/xlf2003" CACHE PATH "")

##############
# TPLs
##############

# Root directory for generated TPLs
set(TPL_ROOT "/usr/workspace/wsrzc/axom/thirdparty_libs/builds/2018_11_14_15_17_37/clang-coral_xlf" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.19" CACHE PATH "")

# scr not built by uberenv

# conduit from uberenv
set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.3.1" CACHE PATH "")

# mfem from uberenv
set(MFEM_DIR "${TPL_ROOT}/mfem-3.4.0" CACHE PATH "")

# python not built by uberenv

set(ENABLE_DOCS OFF CACHE BOOL "")

# shroud not built by uberenv

# uncrustify not built by uberenv

# lcov and genhtml not built by uberenv

##############
# MPI
##############

set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_C_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-clang-coral-2018.05.23/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-clang-coral-2018.05.23/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-clang-coral-2018.05.23/bin/mpif90" CACHE PATH "")

set(MPIEXEC "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-clang-coral-2018.05.23/bin/mpirun" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-np" CACHE PATH "")

set(BLT_MPI_COMMAND_APPEND "mpibind" CACHE PATH "")

##############
# Other machine specifics
##############

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

set(CMAKE_Fortran_COMPILER_ID "XL" CACHE PATH "All of BlueOS compilers report clang due to nvcc, override to proper compiler family")

set(BLT_FORTRAN_FLAGS "-WF,-C!" CACHE PATH "Converts C-style comments to Fortran style in preprocessed files")

set(BLT_EXE_LINKER_FLAGS "-Wl,-rpath,/usr/tce/packages/xl/xl-2018.05.18/lib/" CACHE PATH "Adds a missing rpath for libraries associated with the fortran compiler")


