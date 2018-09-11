##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# blueos_3_ppc64le_ib-gcc@4.9.3
##################################

# cmake from uberenv
# cmake executable path: /usr/workspace/wsa/axom/thirdparty_libs/builds/2018_08_08_11_27_15/spack/opt/spack/blueos_3_ppc64le_ib/gcc-4.9.3/cmake-3.9.6-i462ncjrxqhqrpkmns3pxxmuo4b7qkpe/bin/cmake

#######
# using gcc@4.9.3 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/tcetmp/packages/gcc/gcc-4.9.3/bin/gcc" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/tcetmp/packages/gcc/gcc-4.9.3/bin/g++" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE BOOL "")

set(CMAKE_Fortran_COMPILER "/usr/tcetmp/packages/gcc/gcc-4.9.3/bin/gfortran" CACHE PATH "")

# Root directory for generated TPLs
set(TPL_ROOT "/usr/workspace/wsa/axom/thirdparty_libs/builds/2018_08_08_11_27_15/spack/opt/spack/blueos_3_ppc64le_ib/gcc-4.9.3" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.16-nj464kx3de5a7cnt3rzn57lltdku4x57" CACHE PATH "")

# scr not built by uberenv

# conduit from uberenv
set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.3.1-xppzt3nyq2ngpjwtq3ojlei3dnx65rrd" CACHE PATH "")

# mfem from uberenv
set(MFEM_DIR "${TPL_ROOT}/mfem-3.3.2-4fk5pquwci5yergo67eysgtstybmyme6" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "${TPL_ROOT}/python-2.7.15-fkftablnp4lm263zhkbxfxkoapzlnozo/bin/python" CACHE PATH "")

# doxygen from uberenv
set(DOXYGEN_EXECUTABLE "${TPL_ROOT}/doxygen-1.8.11-7o4ksb5pi6bk6gpfbyqk6to54vfqgy72/bin/doxygen" CACHE PATH "")

# sphinx 1.4.5 from uberenv
set(SPHINX_EXECUTABLE "${TPL_ROOT}/python-2.7.15-fkftablnp4lm263zhkbxfxkoapzlnozo/bin/sphinx-build" CACHE PATH "")

# shroud 0.10.1 from uberenv
set(SHROUD_EXECUTABLE "${TPL_ROOT}/python-2.7.15-fkftablnp4lm263zhkbxfxkoapzlnozo/bin/shroud" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "${TPL_ROOT}/uncrustify-0.61-4yckldhiqow4hnhdp64qy6fswzak5huu/bin/uncrustify" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "${TPL_ROOT}/lcov-1.11-mjtustgukaekv6jjdukymxyrkvp3w2qt/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "${TPL_ROOT}/lcov-1.11-mjtustgukaekv6jjdukymxyrkvp3w2qt/usr/bin/genhtml" CACHE PATH "")

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
# lc blueos gcc@4.9.3  host configs
##############################################################################

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

##############################################################################
# MPI - manually added for now
##############################################################################

set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_HOME                 "/usr/tce/packages/spectrum-mpi/spectrum-mpi-2018.04.27-gcc-4.9.3" CACHE PATH "")
set(MPI_C_COMPILER           "${MPI_HOME}/bin/mpicc"   CACHE PATH "")
set(MPI_CXX_COMPILER         "${MPI_HOME}/bin/mpicxx"  CACHE PATH "")
set(MPI_Fortran_COMPILER     "${MPI_HOME}/bin/mpif90"  CACHE PATH "")

set(MPIEXEC                "${MPI_HOME}/bin/mpirun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG   "-np" CACHE PATH "")
set(BLT_MPI_COMMAND_APPEND "mpibind" CACHE PATH "")
##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

