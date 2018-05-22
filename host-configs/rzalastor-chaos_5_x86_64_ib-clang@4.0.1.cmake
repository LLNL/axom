##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# chaos_5_x86_64_ib-clang@4.0.1
##################################

# cmake from uberenv
# cmake executable path: /usr/workspace/wsrzc/axom/thirdparty_libs/builds/2018_05_01_21_02_22/spack/opt/spack/chaos_5_x86_64_ib/clang-4.0.1/cmake-3.8.2-ogx3v5xe7tjiworftog2b65rue7fj5tu/bin/cmake

#######
# using clang@4.0.1 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/global/tools/clang/chaos_5_x86_64_ib/clang-4.0.1/bin/clang" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/global/tools/clang/chaos_5_x86_64_ib/clang-4.0.1/bin/clang++" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE BOOL "")

set(CMAKE_Fortran_COMPILER "/usr/apps/gnu/4.9.3/bin/gfortran" CACHE PATH "")

# Root directory for generated TPLs
set(TPL_ROOT "/usr/workspace/wsrzc/axom/thirdparty_libs/builds/2018_05_01_21_02_22/spack/opt/spack/chaos_5_x86_64_ib/clang-4.0.1" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.16-f354ty3s3eovs5voooegrffh3yg5zde5" CACHE PATH "")

# scr not built by uberenv

# conduit from uberenv
set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.3.1-v2xb25nhtnm2cbmkupq4z3rv53hdxbth" CACHE PATH "")

# mfem from uberenv
set(MFEM_DIR "${TPL_ROOT}/mfem-3.3.2-5ilanyoky7ep7egpzp5cttp5n3t3mehe" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "${TPL_ROOT}/python-2.7.11-jeoiapalyouzfiwq45njnixltqvkqzi6/bin/python" CACHE PATH "")

# lua from uberenv
set(LUA_DIR "${TPL_ROOT}/lua-5.1.5-uhbj2j5qw5rjczce5rner57grpvowwp5" CACHE PATH "")

# doxygen from uberenv
set(DOXYGEN_EXECUTABLE "${TPL_ROOT}/doxygen-1.8.11-os5b5ulif6d2mo2cfbnjqg6q46vxsw74/bin/doxygen" CACHE PATH "")

# sphinx 1.4.5 from uberenv
set(SPHINX_EXECUTABLE "${TPL_ROOT}/python-2.7.11-jeoiapalyouzfiwq45njnixltqvkqzi6/bin/sphinx-build" CACHE PATH "")

# shroud 0.9.0 from uberenv
set(SHROUD_EXECUTABLE "${TPL_ROOT}/python-2.7.11-jeoiapalyouzfiwq45njnixltqvkqzi6/bin/shroud" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "${TPL_ROOT}/uncrustify-0.61-soshhb7er6sanbcrag2i2a4zts3f6vt5/bin/uncrustify" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "${TPL_ROOT}/lcov-1.11-x6zdnfkjncb5qu35zrpzmcfol5xtik63/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "${TPL_ROOT}/lcov-1.11-x6zdnfkjncb5qu35zrpzmcfol5xtik63/usr/bin/genhtml" CACHE PATH "")

# Disable CXX11 on chaos5 intel/clang builds
set(BLT_CXX_STD "c++98" CACHE PATH "")

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
# lc chaos5 clang@4.0.1  host configs
##############################################################################

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

##############################################################################
# MPI - manually added for now.
##############################################################################
set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_HOME             "/usr/global/tools/clang/chaos_5_x86_64_ib/clang-4.0.1" CACHE PATH "")
set(MPI_C_COMPILER       "${MPI_HOME}/bin/mpiclang" CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpiclang++" CACHE PATH "")

set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

