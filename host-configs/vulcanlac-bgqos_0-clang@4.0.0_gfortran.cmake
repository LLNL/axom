##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# bgqos_0-clang@4.0.0_gfortran
##################################

# cmake from uberenv
# cmake executable path: /collab/usr/global/tools/cmake/bgqos_0/cmake-3.8.2/bin/cmake

#######
# using clang@4.0.0_gfortran compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/collab/usr/gapps/opnsrc/gnu/dev/lnx-2.12-ppc/bgclang/r284961-stable/llnl/bin/mpiclang" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/collab/usr/gapps/opnsrc/gnu/dev/lnx-2.12-ppc/bgclang/r284961-stable/llnl/bin/mpiclang++" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE BOOL "")

set(CMAKE_Fortran_COMPILER "/usr/local/tools/toolchain-4.8.4/gnu-linux-4.8.4/bin/powerpc64-bgq-linux-gfortran" CACHE PATH "")

# Root directory for generated TPLs
set(TPL_ROOT "/usr/workspace/wsa/axom/thirdparty_libs/builds/2018_05_01_15_59_28/spack/opt/spack/bgqos_0/clang-4.0.0_gfortran" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.16-u72ucrahnfkjno3gj3aicmo3iyu7czfg" CACHE PATH "")

# scr not built by uberenv

# conduit from uberenv
set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.3.1-hrqzeb6efxqetvdawg3fokthpk3d4gfg" CACHE PATH "")

# mfem from uberenv
set(MFEM_DIR "${TPL_ROOT}/mfem-3.3.2-niahyxq4gel6yee3vhsaij3aoefvz3ks" CACHE PATH "")

# python not built by uberenv

# lua not built by uberenv

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
# lc bgq clang@4.0.0_gfortran host configs
##############################################################################

set(ENABLE_DOCS    OFF CACHE BOOL "")

set(CMAKE_SKIP_RPATH TRUE CACHE BOOL "")

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

