##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# chaos_5_x86_64_ib-intel@15.0.187
##################################

# cmake from uberenv
# cmake executable path: /usr/workspace/wsrzc/axom/thirdparty_libs/builds/2017_11_13_15_58_48/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.187/cmake-3.8.2-j2jm67xewe7ranc7e55soktrlchvey3k/bin/cmake

#######
# using intel@15.0.187 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/local/tools/ic-15.0.187/bin/icc" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/local/tools/ic-15.0.187/bin/icpc" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE BOOL "")

set(CMAKE_Fortran_COMPILER "/usr/local/tools/ic-15.0.187/bin/ifort" CACHE PATH "")

# Root directory for generated TPLs
set(TPL_ROOT "/usr/workspace/wsrzc/axom/thirdparty_libs/builds/2017_11_13_15_58_48/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.187" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.16-dvlk5aqfby3atjbqefigk6ugby5dk5yy" CACHE PATH "")

# scr not built by uberenv

# conduit from uberenv
set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.2.1-ykk4uxptqllgfoeevc7k3s5gabn6txhm" CACHE PATH "")

# mfem from uberenv
set(MFEM_DIR "${TPL_ROOT}/mfem-3.3.2-64q6i7tvzthligr75o7fuyy2spkmphth" CACHE PATH "")

# boost headers from uberenv
set(BOOST_DIR "${TPL_ROOT}/boost-headers-1.58.0-m4zq524o33mi34qimoxtn5rae4amyk6w" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "${TPL_ROOT}/python-2.7.11-2us32xlalxdytptzzhfiniiimqrf4ezh/bin/python" CACHE PATH "")

# lua from uberenv
set(LUA_DIR "${TPL_ROOT}/lua-5.1.5-7gpjuo5khxsxqwidz74wccitsl3hdixj" CACHE PATH "")

# doxygen from uberenv
set(DOXYGEN_EXECUTABLE "${TPL_ROOT}/doxygen-1.8.11-4nntybfqq6li3nbbrjkhw47li2vlj6a6/bin/doxygen" CACHE PATH "")

# sphinx from uberenv
set(SPHINX_EXECUTABLE "${TPL_ROOT}/python-2.7.11-2us32xlalxdytptzzhfiniiimqrf4ezh/bin/sphinx-build" CACHE PATH "")

# shroud from uberenv
set(SHROUD_EXECUTABLE "${TPL_ROOT}/python-2.7.11-2us32xlalxdytptzzhfiniiimqrf4ezh/bin/shroud" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "${TPL_ROOT}/uncrustify-0.61-onqgsfbevtniuj45s564dpfkmlveql33/bin/uncrustify" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "${TPL_ROOT}/lcov-1.11-hyyvo3em3ckyvwg5bz27koz2ulmsbms3/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "${TPL_ROOT}/lcov-1.11-hyyvo3em3ckyvwg5bz27koz2ulmsbms3/usr/bin/genhtml" CACHE PATH "")

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
# lc chaos5 intel@15.0.187 host configs
##############################################################################

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

##############################################################################
# MPI - manually added for now
##############################################################################
set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_HOME             "/usr/local/tools/mvapich2-intel-2.0" CACHE PATH "")
set(MPI_C_COMPILER       "${MPI_HOME}/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpicxx" CACHE PATH "")
set(MPI_Fortran_COMPILER "${MPI_HOME}/bin/mpif90" CACHE PATH "")

set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

