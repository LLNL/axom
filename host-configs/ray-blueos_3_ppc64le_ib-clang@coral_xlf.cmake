##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# blueos_3_ppc64le_ib-clang@coral_xlf
##################################

# cmake from uberenv
# cmake executable path: /usr/workspace/wsa/axom/thirdparty_libs/builds/2018_05_01_15_58_27/spack/opt/spack/blueos_3_ppc64le_ib/clang-coral_xlf/cmake-3.8.2-f4bp77lrsuauuncgtivaoox2vfombs32/bin/cmake

#######
# using clang@coral_xlf compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/tce/packages/clang/clang-coral-2017.10.13/bin/clang" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/tce/packages/clang/clang-coral-2017.10.13/bin/clang++" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE BOOL "")

set(CMAKE_Fortran_COMPILER "/usr/tce/packages/xl/xl-beta-2017.10.13/bin/xlf2003" CACHE PATH "")

# Root directory for generated TPLs
set(TPL_ROOT "/usr/workspace/wsa/axom/thirdparty_libs/builds/2018_05_01_15_58_27/spack/opt/spack/blueos_3_ppc64le_ib/clang-coral_xlf" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.16-73nxx5224rq5zpn4kjrkq2hfmuzu4hmt" CACHE PATH "")

# scr not built by uberenv

# conduit from uberenv
set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.3.1-sl2vrhuyfhhq5r75yoikb2a255skzxcz" CACHE PATH "")

# mfem from uberenv
set(MFEM_DIR "${TPL_ROOT}/mfem-3.3.2-lmjbjzcrjmwvgpgavjki62z4fa42c4bu" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "${TPL_ROOT}/python-2.7.11-2j4lnnqrepn6p2jci2pqv3qtietk5nh4/bin/python" CACHE PATH "")

# lua from uberenv
set(LUA_DIR "${TPL_ROOT}/lua-5.1.5-ff6lol2zlfsqiebnrseluf2kgsbr36fo" CACHE PATH "")

# doxygen from uberenv
set(DOXYGEN_EXECUTABLE "${TPL_ROOT}/doxygen-1.8.11-4ohg4xoyg4gwlmbdjaorce7rrhkzdkds/bin/doxygen" CACHE PATH "")

# sphinx 1.4.5 from uberenv
set(SPHINX_EXECUTABLE "${TPL_ROOT}/python-2.7.11-2j4lnnqrepn6p2jci2pqv3qtietk5nh4/bin/sphinx-build" CACHE PATH "")

# shroud 0.9.0 from uberenv
set(SHROUD_EXECUTABLE "${TPL_ROOT}/python-2.7.11-2j4lnnqrepn6p2jci2pqv3qtietk5nh4/bin/shroud" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "${TPL_ROOT}/uncrustify-0.61-xnpbzzqgafqu22rrdooflush2kdv5yqv/bin/uncrustify" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "${TPL_ROOT}/lcov-1.11-e4socip3p6sucjy3o77n5jmvzu2mqdfk/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "${TPL_ROOT}/lcov-1.11-e4socip3p6sucjy3o77n5jmvzu2mqdfk/usr/bin/genhtml" CACHE PATH "")

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
# lc blueos clang@coral_xlf host configs
##############################################################################

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

# Converts C-style comments to Fortran style in preprocessed files
set(BLT_FORTRAN_FLAGS "-WF,-C!" CACHE STRING "")

# Adds a missing rpath for libraries associated with the fortran compiler
set(BLT_EXE_LINKER_FLAGS "-Wl,-rpath,/usr/tce/packages/xl/xl-beta-2017.10.13/lib" CACHE STRING "")

##############################################################################
# MPI - manually added for now
##############################################################################

set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_HOME                 "/usr/tce/packages/spectrum-mpi/spectrum-mpi-2017.08.24-clang-coral-2017.10.13/" CACHE PATH "")
set(MPI_C_COMPILER           "${MPI_HOME}/bin/mpicc"   CACHE PATH "")
set(MPI_CXX_COMPILER         "${MPI_HOME}/bin/mpicxx"  CACHE PATH "")
set(MPI_Fortran_COMPILER     "${MPI_HOME}/bin/mpifort" CACHE PATH "")

set(MPIEXEC              "mpirun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-np" CACHE PATH "")


##############################################################################
# !---------------------------------------------------------------------------
##############################################################################


