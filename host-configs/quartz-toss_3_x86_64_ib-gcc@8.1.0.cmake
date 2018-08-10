##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# toss_3_x86_64_ib-gcc@8.1.0
##################################

# cmake from uberenv
# cmake executable path: /usr/workspace/wsa/axom/thirdparty_libs/builds/2018_08_09_16_30_09/spack/opt/spack/toss_3_x86_64_ib/gcc-8.1.0/cmake-3.9.6-uus2nf4pd6al7iept4gob7blnaqj2e34/bin/cmake

#######
# using gcc@8.1.0 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-8.1.0/bin/gcc" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/tce/packages/gcc/gcc-8.1.0/bin/g++" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE BOOL "")

set(CMAKE_Fortran_COMPILER "/usr/tce/packages/gcc/gcc-8.1.0/bin/gfortran" CACHE PATH "")

# Root directory for generated TPLs
set(TPL_ROOT "/usr/workspace/wsa/axom/thirdparty_libs/builds/2018_08_09_16_30_09/spack/opt/spack/toss_3_x86_64_ib/gcc-8.1.0" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.16-g6qr3fj6tol47g4nruywqvxrmgdpbl5u" CACHE PATH "")

# scr not built by uberenv

# conduit from uberenv
set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.3.1-shqx436yttimbdiiyy4yjbhjzkpviyj4" CACHE PATH "")

# mfem from uberenv
set(MFEM_DIR "${TPL_ROOT}/mfem-3.3.2-lkdmj5s5sw36irno7bz4rxbqc7x6qzlc" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "${TPL_ROOT}/python-2.7.15-gqyr5hii5aplhcqwcc7bfzvri7msyfjs/bin/python" CACHE PATH "")

# doxygen from uberenv
set(DOXYGEN_EXECUTABLE "${TPL_ROOT}/doxygen-1.8.11-sjb3alpfond3pyfcwpufbvrsve7jb55n/bin/doxygen" CACHE PATH "")

# sphinx 1.4.5 from uberenv
set(SPHINX_EXECUTABLE "${TPL_ROOT}/python-2.7.15-gqyr5hii5aplhcqwcc7bfzvri7msyfjs/bin/sphinx-build" CACHE PATH "")

# shroud 0.10.1 from uberenv
set(SHROUD_EXECUTABLE "${TPL_ROOT}/python-2.7.15-gqyr5hii5aplhcqwcc7bfzvri7msyfjs/bin/shroud" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "${TPL_ROOT}/uncrustify-0.61-e3a4vrkbigclb4w5gmmkhnqa35ahrjzk/bin/uncrustify" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "${TPL_ROOT}/lcov-1.11-47te5uagloz3qjlc6qcjueld3d6jtjw2/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "${TPL_ROOT}/lcov-1.11-47te5uagloz3qjlc6qcjueld3d6jtjw2/usr/bin/genhtml" CACHE PATH "")

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
# lc toss3 gcc@8.1.0  host configs
##############################################################################

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

##############################################################################
# MPI - manually added for now
##############################################################################
set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_HOME             "/usr/tce/packages/mvapich2/mvapich2-2.2-gcc-8.1.0" CACHE PATH "")
set(MPI_C_COMPILER       "${MPI_HOME}/bin/mpicc"   CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpicxx"  CACHE PATH "")
set(MPI_Fortran_COMPILER "${MPI_HOME}/bin/mpifort" CACHE PATH "")

set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

