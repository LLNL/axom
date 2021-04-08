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
# cmake executable path: /home/taylor16/tpl/v7mpi/spack/opt/spack/x86_64/gcc-4.9.3/cmake-3.3.1-jkk5puodwd24sh2egavyrp2ony2z3wtn/bin/cmake

#------------------------------------------------------------------------------
# using gcc@4.9.3 compiler spec
#------------------------------------------------------------------------------

# c compiler used by spack
set(CMAKE_C_COMPILER "/home/taylor16/gapps/gcc-4.9.3/bin/gcc" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/home/taylor16/gapps/gcc-4.9.3/bin/g++" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE PATH "")

set(CMAKE_Fortran_COMPILER  "/home/taylor16/gapps/gcc-4.9.3/bin/gfortran" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "/home/taylor16/tpl/v7mpi/spack/opt/spack/x86_64/gcc-4.9.3/hdf5-1.8.16-vhefrjw2f57yhjzb77fuvhlaowx66dht" CACHE PATH "")

# conduit from uberenv
#set(CONDUIT_DIR "/home/taylor16/tpl/v7mpi/spack/opt/spack/x86_64/gcc-4.9.3/conduit-github-sfgmvgy7nlranrdldpoqrfkrcwv3bvid" CACHE PATH "")
set(CONDUIT_DIR "/home/taylor16/conduit/debug-install" CACHE PATH "")

# doxygen from uberenv
set(DOXYGEN_EXECUTABLE "/home/taylor16/tpl/v7mpi/spack/opt/spack/x86_64/gcc-4.9.3/doxygen-1.8.10-7qxnjx5i5cyy2nbygrdihsdgzrnfvs77/bin/doxygen" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "/home/taylor16/tpl/v7mpi/spack/opt/spack/x86_64/gcc-4.9.3/python-2.7.11-hk7h4ntkiplj3ukfwxcbr7oxqtf3n2hy/bin/python" CACHE PATH "")

# lua from uberenv
set(LUA_DIR "/home/taylor16/tpl/v7mpi/spack/opt/spack/x86_64/gcc-4.9.3/lua-5.1.5-44ijmjynmfovg6eteacz54gu6muabngc" CACHE PATH "")

# sphinx from uberenv
set(SPHINX_EXECUTABLE "/home/taylor16/tpl/v7mpi/spack/opt/spack/x86_64/gcc-4.9.3/python-2.7.11-hk7h4ntkiplj3ukfwxcbr7oxqtf3n2hy/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "/home/taylor16/tpl/v7mpi/spack/opt/spack/x86_64/gcc-4.9.3/uncrustify-0.61-sdasmocwkfmf6xupjpoo73rzdddf2nl6/bin/uncrustify" CACHE PATH "")

# sparsehash headers from uberenv
set(SPARSEHASH_DIR "/home/taylor16/tpl/v7mpi/spack/opt/spack/x86_64/gcc-4.9.3/sparsehash-headers-2.0.2-rmpk4wmjkozch4i54wa7l5zjfyru7yr3" CACHE PATH "")

# boost headers from uberenv
set(ENABLE_BOOST ON CACHE PATH "")
set(BOOST_ROOT "/home/taylor16/tpl/v7mpi/spack/opt/spack/x86_64/gcc-4.9.3/boost-headers-1.58.0-6jjq5njm5itr67imya3p4fcnmiaj5nab" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "/home/taylor16/tpl/v7mpi/spack/opt/spack/x86_64/gcc-4.9.3/lcov-1.11-h3ijl43wxo53yrr7j7gsnmpuc4wc5nrn/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "/home/taylor16/tpl/v7mpi/spack/opt/spack/x86_64/gcc-4.9.3/lcov-1.11-h3ijl43wxo53yrr7j7gsnmpuc4wc5nrn/usr/bin/genhtml" CACHE PATH "")

#------------------------------------------------------------------------------
# end uberenv host-config
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# MPI - manually added these for now.
#------------------------------------------------------------------------------
set(ENABLE_MPI ON CACHE PATH "")
set(MPI_C_COMPILER "/home/taylor16/local/mpich-3.2/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "/home/taylor16/local/mpich-3.2/bin/mpicxx" CACHE PATH "")
set(MPI_Fortran_COMPILER "/home/taylor16/local/mpich-3.2/bin/mpif90" CACHE PATH "")
