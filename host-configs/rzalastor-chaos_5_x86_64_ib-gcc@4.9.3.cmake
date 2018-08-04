##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# chaos_5_x86_64_ib-gcc@4.9.3
##################################

# cmake from uberenv
# cmake executable path: /usr/workspace/wsrzc/axom/thirdparty_libs/builds/2018_08_01_22_53_23/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/cmake-3.9.6-qevyrceeg47ictjkx7vic7x7xbmvzkxj/bin/cmake

#######
# using gcc@4.9.3 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/apps/gnu/4.9.3/bin/gcc" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/apps/gnu/4.9.3/bin/g++" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE BOOL "")

set(CMAKE_Fortran_COMPILER "/usr/apps/gnu/4.9.3/bin/gfortran" CACHE PATH "")

# Root directory for generated TPLs
set(TPL_ROOT "/usr/workspace/wsrzc/axom/thirdparty_libs/builds/2018_08_01_22_53_23/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.16-2xgsiskrcxzxgs7avics5ghpclty2cr6" CACHE PATH "")

# scr not built by uberenv

# conduit from uberenv
set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.3.1-a53ng434ag6q7gg6pquhzsmqqx6ghxjg" CACHE PATH "")

# mfem from uberenv
set(MFEM_DIR "${TPL_ROOT}/mfem-3.3.2-fcmuhkaq7537pk3y62q4fbetyvuuokxr" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "${TPL_ROOT}/python-2.7.15-kfgq6bl5v33eugneh65oenwjgv5u3ar6/bin/python" CACHE PATH "")

# doxygen from uberenv
set(DOXYGEN_EXECUTABLE "${TPL_ROOT}/doxygen-1.8.11-yupy2hfxy5rmhxmbtx52hrvdkbjc6alq/bin/doxygen" CACHE PATH "")

# sphinx 1.4.5 from uberenv
set(SPHINX_EXECUTABLE "${TPL_ROOT}/python-2.7.15-kfgq6bl5v33eugneh65oenwjgv5u3ar6/bin/sphinx-build" CACHE PATH "")

# shroud 0.10.0 from uberenv
set(SHROUD_EXECUTABLE "${TPL_ROOT}/python-2.7.15-kfgq6bl5v33eugneh65oenwjgv5u3ar6/bin/shroud" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "${TPL_ROOT}/uncrustify-0.61-px2meiscmkbwcnmmom3qnlzdzmf2yx7x/bin/uncrustify" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "${TPL_ROOT}/lcov-1.11-dc44m4x5flevrfvqwytiscdc3dody4zx/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "${TPL_ROOT}/lcov-1.11-dc44m4x5flevrfvqwytiscdc3dody4zx/usr/bin/genhtml" CACHE PATH "")

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
# lc chaos5 gcc@4.9.3  host configs
##############################################################################

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

##############################################################################
# MPI - manually added for now
##############################################################################
set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_HOME             "/usr/local/tools/mvapich2-gnu-2.0" CACHE PATH "")
set(MPI_C_COMPILER       "${MPI_HOME}/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpicxx" CACHE PATH "")
set(MPI_Fortran_COMPILER "${MPI_HOME}/bin/mpif90" CACHE PATH "")

set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

##############################################################################
# GCOV - manually added for now
##############################################################################
set(GCOV_PATH "/usr/apps/gnu/4.9.3/bin/gcov" CACHE PATH "")

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

