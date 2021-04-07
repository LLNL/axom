# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#------------------------------------------------------------------------------
# uberenv host-config
#
# This is a generated file, edit at own risk.
#------------------------------------------------------------------------------
# x86_64-gcc@4.9.3
#------------------------------------------------------------------------------

# cmake from uberenv
# cmake executable path: /home/taylor16/apps/cmake-3.8.2/bin/cmake
#------------------------------------------------------------------------------
# using gcc@4.9.3 compiler spec
#------------------------------------------------------------------------------

# c compiler used by spack
set(CMAKE_C_COMPILER "/home/taylor16/apps/gcc-4.9.4/bin/gcc" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/home/taylor16/apps/gcc-4.9.4/bin/g++" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE PATH "")

set(CMAKE_Fortran_COMPILER  "/home/taylor16/apps/gcc-4.9.4/bin/gfortran" CACHE PATH "")

# Root directory for generated TPLs
set(TPL_ROOT "/home/taylor16/tpl/v2" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "${TPL_ROOT}" CACHE PATH "")

# conduit from uberenv
set(CONDUIT_DIR "${TPL_ROOT}" CACHE PATH "")

# doxygen from uberenv
set(DOXYGEN_EXECUTABLE "${TPL_ROOT}/bin/doxygen" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "${TPL_ROOT}/bin/python" CACHE PATH "")

# lua from uberenv
set(LUA_DIR "${TPL_ROOT}" CACHE PATH "")

# sphinx from uberenv
set(SPHINX_EXECUTABLE "${TPL_ROOT}/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "${TPL_ROOT}/bin/uncrustify" CACHE PATH "")

# boost headers from uberenv
set(BOOST_DIR "${TPL_ROOT}" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "${TPL_ROOT}/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "${TPL_ROOT}/usr/bin/genhtml" CACHE PATH "")

#------------------------------------------------------------------------------
# end uberenv host-config
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# MPI - manually added these for now.
#------------------------------------------------------------------------------
#set(ENABLE_MPI ON CACHE PATH "")
#set(MPI_C_COMPILER "/home/taylor16/local/mpich-3.2/bin/mpicc" CACHE PATH "")
#set(MPI_CXX_COMPILER "/home/taylor16/local/mpich-3.2/bin/mpicxx" CACHE PATH "")
#set(MPI_Fortran_COMPILER "/home/taylor16/local/mpich-3.2/bin/mpif90" CACHE PATH "")

#------------------------------------------------------------------------------
# SHROUD - manually added for now. Use a public build add to TPL later
#------------------------------------------------------------------------------
set(SHROUD_EXECUTABLE "/home/taylor16/tpl/shroud/bin/shroud" CACHE PATH "")

