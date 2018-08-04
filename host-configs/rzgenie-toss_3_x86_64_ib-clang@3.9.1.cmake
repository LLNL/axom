##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# toss_3_x86_64_ib-clang@3.9.1
##################################

# cmake from uberenv
# cmake executable path: /usr/workspace/wsrzc/axom/thirdparty_libs/builds/2018_08_03_13_15_33/spack/opt/spack/toss_3_x86_64_ib/clang-3.9.1/cmake-3.9.6-kx3nvyyx3n57gqudnudr65uxq5soykls/bin/cmake

#######
# using clang@3.9.1 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/tce/packages/clang/clang-3.9.1/bin/clang" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/tce/packages/clang/clang-3.9.1/bin/clang++" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE BOOL "")

set(CMAKE_Fortran_COMPILER "/usr/tce/packages/gcc/gcc-4.9.3/bin/gfortran" CACHE PATH "")

# Root directory for generated TPLs
set(TPL_ROOT "/usr/workspace/wsrzc/axom/thirdparty_libs/builds/2018_08_03_13_15_33/spack/opt/spack/toss_3_x86_64_ib/clang-3.9.1" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.16-svvk7kwygp744sbxd6riqgw5xve72sb4" CACHE PATH "")

# scr not built by uberenv

# conduit from uberenv
set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.3.1-ns53fi42bezk2i7frzslow7hmiun5oai" CACHE PATH "")

# mfem from uberenv
set(MFEM_DIR "${TPL_ROOT}/mfem-3.3.2-25tyrwypoyajb7st7a63n4snbeh4jkbt" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "${TPL_ROOT}/python-2.7.15-5ovgvdg56teo36el2vcmkp4atspgonwv/bin/python" CACHE PATH "")

# doxygen from uberenv
set(DOXYGEN_EXECUTABLE "${TPL_ROOT}/doxygen-1.8.11-gvzzxscgq446lxa57bjmoewkxu52cv7k/bin/doxygen" CACHE PATH "")

# sphinx 1.4.5 from uberenv
set(SPHINX_EXECUTABLE "${TPL_ROOT}/python-2.7.15-5ovgvdg56teo36el2vcmkp4atspgonwv/bin/sphinx-build" CACHE PATH "")

# shroud 0.10.0 from uberenv
set(SHROUD_EXECUTABLE "${TPL_ROOT}/python-2.7.15-5ovgvdg56teo36el2vcmkp4atspgonwv/bin/shroud" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "${TPL_ROOT}/uncrustify-0.61-yfckh2c72c4wx7xpuokpqpqmozyg2zpv/bin/uncrustify" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "${TPL_ROOT}/lcov-1.11-tnhlnv4vfsj7yyq5krefxrm7n4utyyhh/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "${TPL_ROOT}/lcov-1.11-tnhlnv4vfsj7yyq5krefxrm7n4utyyhh/usr/bin/genhtml" CACHE PATH "")

##################################
# end uberenv host-config
##################################

##
## Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
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
# lc toss3 clang@3.9.1  host configs
##############################################################################

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

##############################################################################
# MPI - manually added for now
##############################################################################

set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_HOME             "/usr/tce/packages/mvapich2/mvapich2-2.2-clang-3.9.1" CACHE PATH "")
set(MPI_C_COMPILER       "${MPI_HOME}/bin/mpicc"   CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpicxx"  CACHE PATH "")
set(MPI_Fortran_COMPILER "${MPI_HOME}/bin/mpifort" CACHE PATH "")

set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

