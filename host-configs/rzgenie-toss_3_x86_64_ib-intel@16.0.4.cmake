##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# toss_3_x86_64_ib-intel@16.0.4
##################################

# cmake from uberenv
# cmake executable path: /usr/workspace/wsrzc/axom/thirdparty_libs/builds/2018_08_03_13_15_33/spack/opt/spack/toss_3_x86_64_ib/intel-16.0.4/cmake-3.9.6-rl2m2j5ai7tcpew6vflqaduxh3aerc3o/bin/cmake

#######
# using intel@16.0.4 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/tce/packages/intel/intel-16.0.4/bin/icc" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/tce/packages/intel/intel-16.0.4/bin/icpc" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE BOOL "")

set(CMAKE_Fortran_COMPILER "/usr/tce/packages/intel/intel-16.0.4/bin/ifort" CACHE PATH "")

# Root directory for generated TPLs
set(TPL_ROOT "/usr/workspace/wsrzc/axom/thirdparty_libs/builds/2018_08_03_13_15_33/spack/opt/spack/toss_3_x86_64_ib/intel-16.0.4" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.16-z6nt2tyk2j2eop2df6ghod6dnfhwei6s" CACHE PATH "")

# scr not built by uberenv

# conduit from uberenv
set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.3.1-lztognjko5w4lu4krzqu5aehgjtjg7bj" CACHE PATH "")

# mfem from uberenv
set(MFEM_DIR "${TPL_ROOT}/mfem-3.3.2-3pdfiuw56frxlg5w44ac2ovcexy7qd4o" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "${TPL_ROOT}/python-2.7.15-pdkeouqvzqrfu333xcknssxpyfvftps7/bin/python" CACHE PATH "")

# doxygen from uberenv
set(DOXYGEN_EXECUTABLE "${TPL_ROOT}/doxygen-1.8.11-3lrulsiwppjyanbjnzilafunn3shvy5s/bin/doxygen" CACHE PATH "")

# sphinx 1.4.5 from uberenv
set(SPHINX_EXECUTABLE "${TPL_ROOT}/python-2.7.15-pdkeouqvzqrfu333xcknssxpyfvftps7/bin/sphinx-build" CACHE PATH "")

# shroud 0.10.0 from uberenv
set(SHROUD_EXECUTABLE "${TPL_ROOT}/python-2.7.15-pdkeouqvzqrfu333xcknssxpyfvftps7/bin/shroud" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "${TPL_ROOT}/uncrustify-0.61-nzydqv27ioishmus56myotierqugxn54/bin/uncrustify" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "${TPL_ROOT}/lcov-1.11-fiwdecmhqvashhs323poq5mdr6uqnrio/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "${TPL_ROOT}/lcov-1.11-fiwdecmhqvashhs323poq5mdr6uqnrio/usr/bin/genhtml" CACHE PATH "")

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
# lc toss3 intel@16.0.4  host configs
##############################################################################

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

##############################################################################
# MPI - manually added for now
##############################################################################
set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_HOME             "/usr/tce/packages/mvapich2/mvapich2-2.2-intel-16.0.4" CACHE PATH "")
set(MPI_C_COMPILER       "${MPI_HOME}/bin/mpicc"   CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpicxx"  CACHE PATH "")
set(MPI_Fortran_COMPILER "${MPI_HOME}/bin/mpifort" CACHE PATH "")

set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

