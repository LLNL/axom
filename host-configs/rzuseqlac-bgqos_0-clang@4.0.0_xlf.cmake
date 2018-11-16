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

# SYS_TYPE: bgqos_0
# Compiler Spec: clang@4.0.0_xlf
##################################

# CMake executable path: /collab/usr/global/tools/cmake/bgqos_0/cmake-3.8.2/bin/cmake

##############
# Compilers
##############

# Note: we build TPLs with the serial compiler then use MPI wrappers on bgq
# Serial compilers used by spack:
# C compiler: /collab/usr/gapps/opnsrc/gnu/dev/lnx-2.12-ppc/bgclang/r284961-stable/bin/bgclang
# C++ compiler: /collab/usr/gapps/opnsrc/gnu/dev/lnx-2.12-ppc/bgclang/r284961-stable/bin/bgclang++

set(CMAKE_C_COMPILER "/collab/usr/gapps/opnsrc/gnu/dev/lnx-2.12-ppc/bgclang/r284961-stable/llnl/bin/mpiclang" CACHE PATH "")

set(CMAKE_CXX_COMPILER "/collab/usr/gapps/opnsrc/gnu/dev/lnx-2.12-ppc/bgclang/r284961-stable/llnl/bin/mpiclang++" CACHE PATH "")

# Fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE BOOL "")

set(CMAKE_Fortran_COMPILER "/opt/ibmcmp/xlf/bg/14.1/bin/bgxlf2003" CACHE PATH "")

##############
# TPLs
##############

# Root directory for generated TPLs
set(TPL_ROOT "/usr/workspace/wsrzc/axom/thirdparty_libs/builds/2018_11_14_15_08_15/clang-4.0.0_xlf" CACHE PATH "")

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

set(ENABLE_FIND_MPI OFF CACHE BOOL "Use wrapper directly to stop FindMPI returning the wrong linker flags.")

set(ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC ON CACHE BOOL "Ensures that tests will be wrapped with srun to run on the backend nodes")

set(MPI_Fortran_INCLUDE_PATH "/usr/local/tools/deg/drivers/V1R2M0/ppc64/comm/gcc/include" CACHE PATH "Pass in an explicit path to help find mpif.h")

set(MPIEXEC "/usr/bin/srun" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

##############
# Other machine specifics
##############

set(ENABLE_GTEST_DEATH_TESTS OFF CACHE BOOL "")

set(BLT_FORTRAN_FLAGS "-WF,-C!" CACHE PATH "Converts C-style comments to Fortran style in preprocessed files")

# Manually set up HDF5 library dependencies for BGQ to bypass errors from CMake's FindHDF5
set(HDF5_C_LIBRARY_m "-lm" CACHE PATH "")

set(HDF5_C_LIBRARY_dl "-ldl" CACHE PATH "")

set(CMAKE_SKIP_RPATH ON CACHE BOOL "")


