##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# bgqos_0-clang@4.0.0_xlf
##################################

# cmake from uberenv
# cmake executable path: /collab/usr/global/tools/cmake/bgqos_0/cmake-3.8.2/bin/cmake

#######
# using clang@4.0.0_xlf compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/collab/usr/gapps/opnsrc/gnu/dev/lnx-2.12-ppc/bgclang/r284961-stable/llnl/bin/mpiclang" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/collab/usr/gapps/opnsrc/gnu/dev/lnx-2.12-ppc/bgclang/r284961-stable/llnl/bin/mpiclang++" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE BOOL "")

set(CMAKE_Fortran_COMPILER "/opt/ibmcmp/xlf/bg/14.1/bin/bgxlf2003" CACHE PATH "")

# Root directory for generated TPLs
set(TPL_ROOT "/usr/workspace/wsrzc/axom/thirdparty_libs/builds/2018_08_03_15_05_44/spack/opt/spack/bgqos_0/clang-4.0.0_xlf" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.16-4e4odzhtauwbuqgh2rkwp3wnljdqrhri" CACHE PATH "")

# scr not built by uberenv

# conduit from uberenv
set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.3.1-rk3gn24gxjozpvoynmsrehwzsilbge7k" CACHE PATH "")

# mfem from uberenv
set(MFEM_DIR "${TPL_ROOT}/mfem-3.3.2-ln22e3ngntx35dyjs3k2ndzyve5bmv5y" CACHE PATH "")

# python not built by uberenv

# doxygen not built by uberenv

# sphinx not built by uberenv

# shroud not built by uberenv

# uncrustify not built by uberenv

# lcov and genhtml not built by uberenv

##################################
# end uberenv host-config
##################################

##
## Copyright (c) 2017, Lawrence Livermore National Security, LLC.
##
## Produced at the Lawrence Livermore National Laboratory.
##
## LLNL-CODE-741217
##
## All rights reserved.
##
## This file is part of Axom.
##
## For details about use and distribution, please read axom/LICENSE.
##

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################
# Options added manually to 
# lc bgq clang@4.0.0_xlf host configs
##############################################################################

set(ENABLE_DOCS    OFF CACHE BOOL "")

set(CMAKE_SKIP_RPATH TRUE CACHE BOOL "")

# Converts C-style comments to Fortran style in preprocessed files
set(BLT_FORTRAN_FLAGS "-WF,-C!" CACHE STRING "")

# Manually set up HDF5 library dependencies for BGQ to bypass errors from CMake's FindHDF5
set(HDF5_C_LIBRARY_m "-lm" CACHE STRING "")
set(HDF5_C_LIBRARY_dl "-ldl" CACHE STRING "")

##############################################################################
# MPI - manually added for now
##############################################################################
set(ENABLE_MPI      ON CACHE BOOL "")

# Note: On BGQ, CMake uses the wrong linker flags when using FindMPI.
# Disabling FindMPI allows us to use the wrapper directly via the CMake compiler variables.
set(ENABLE_FIND_MPI OFF CACHE BOOL "")

# Pass in an explicit path to help find mpif.h
set(MPI_Fortran_INCLUDE_PATH "/usr/local/tools/deg/drivers/V1R2M0/ppc64/comm/gcc/include" CACHE PATH "")

set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

set(ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC TRUE CACHE BOOL "Ensures that tests will be wrapped with srun to run on the backend nodes")

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

