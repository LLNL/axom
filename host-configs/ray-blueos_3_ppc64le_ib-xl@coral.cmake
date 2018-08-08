##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# blueos_3_ppc64le_ib-xl@coral
##################################

# cmake from uberenv
# cmake executable path: /usr/workspace/wsa/axom/thirdparty_libs/builds/2018_06_06_17_07_58/spack/opt/spack/blueos_3_ppc64le_ib/xl-coral/cmake-3.9.6-tgyt2paysrstkcpr2vioirnqrodq66nz/bin/cmake

#######
# using xl@coral compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/tce/packages/xl/xl-2018.05.18/bin/xlc" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/tce/packages/xl/xl-2018.05.18/bin/xlC" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE BOOL "")

set(CMAKE_Fortran_COMPILER "/usr/tce/packages/xl/xl-2018.05.18/bin/xlf2003" CACHE PATH "")

# Root directory for generated TPLs
set(TPL_ROOT "/usr/workspace/wsa/axom/thirdparty_libs/builds/2018_06_06_17_07_58/spack/opt/spack/blueos_3_ppc64le_ib/xl-coral" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.16-m2nzqvxxqop7g6m72zus2iz5dyyljgds" CACHE PATH "")

# scr not built by uberenv

# conduit from uberenv
set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.3.1-aoehnwvmujmlvwbv7smdvp33ymsaepqo" CACHE PATH "")

# mfem from uberenv
set(MFEM_DIR "${TPL_ROOT}/mfem-3.3.2-jts226skzh4yeew5yls3irfcrmo7sfhl" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "${TPL_ROOT}/python-2.7.11-czmmsiz43za3cndaqutjilrozhjwkuma/bin/python" CACHE PATH "")

# doxygen from uberenv
set(DOXYGEN_EXECUTABLE "${TPL_ROOT}/doxygen-1.8.11-wl2hwvov4xx6exbzvljnftuhnw6tuf3s/bin/doxygen" CACHE PATH "")

# sphinx 1.4.5 from uberenv
set(SPHINX_EXECUTABLE "${TPL_ROOT}/python-2.7.11-czmmsiz43za3cndaqutjilrozhjwkuma/bin/sphinx-build" CACHE PATH "")

# shroud 0.9.0 from uberenv
set(SHROUD_EXECUTABLE "${TPL_ROOT}/python-2.7.11-czmmsiz43za3cndaqutjilrozhjwkuma/bin/shroud" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "${TPL_ROOT}/uncrustify-0.61-fzxffsbhyvh7hqdg725ji6x5j7fporn2/bin/uncrustify" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "${TPL_ROOT}/lcov-1.11-x3r6z227aw5uuwtyapupfxwhefls24dc/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "${TPL_ROOT}/lcov-1.11-x3r6z227aw5uuwtyapupfxwhefls24dc/usr/bin/genhtml" CACHE PATH "")

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
# lc blueos xl@coral host configs
##############################################################################

set(CMAKE_C_COMPILER_ID       "XL" CACHE STRING "")
set(CMAKE_CXX_COMPILER_ID     "XL" CACHE STRING "")
set(CMAKE_Fortran_COMPILER_ID "XL" CACHE STRING "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

# Convert C-style comments to Fortran and link fortran exes to C++ libraries
set(BLT_FORTRAN_FLAGS "-WF,-C! -qxlf2003=polymorphic" CACHE STRING "")

##############################################################################
# MPI - manually added for now
##############################################################################

set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_HOME                 "/usr/tce/packages/spectrum-mpi/spectrum-mpi-2018.04.27-xl-2018.05.18/" CACHE PATH "")
set(MPI_C_COMPILER           "${MPI_HOME}/bin/mpixlc"   CACHE PATH "")
set(MPI_CXX_COMPILER         "${MPI_HOME}/bin/mpixlC"   CACHE PATH "")
set(MPI_Fortran_COMPILER     "${MPI_HOME}/bin/mpixlf"   CACHE PATH "")

set(MPIEXEC                  "${MPI_HOME}/bin/mpirun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG     "-np" CACHE PATH "")
set(BLT_MPI_COMMAND_APPEND   "mpibind" CACHE PATH "")

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################


