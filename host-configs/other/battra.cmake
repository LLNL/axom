##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# x86_64-gcc@4.9.3
##################################

# cmake from uberenv
# cmake executable path: /home/taylor16/local/bin/cmake

#######
# using gcc@4.9.3 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/home/taylor16/apps/gcc-4.9.4/bin/gcc" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/home/taylor16/apps/gcc-4.9.4/bin/g++" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE PATH "")

set(CMAKE_Fortran_COMPILER  "/home/taylor16/apps/gcc-4.9.4/bin/gfortran" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "/home/taylor16/tpl/v2" CACHE PATH "")

# conduit from uberenv
set(CONDUIT_DIR "/home/taylor16/tpl/v2" CACHE PATH "")

# doxygen from uberenv
set(DOXYGEN_EXECUTABLE "/home/taylor16/tpl/v2/bin/doxygen" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "/home/taylor16/tpl/v2/bin/python" CACHE PATH "")

# lua from uberenv
set(LUA_DIR "/home/taylor16/tpl/v2" CACHE PATH "")

# sphinx from uberenv
set(SPHINX_EXECUTABLE "/home/taylor16/tpl/v2/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "/home/taylor16/tpl/v2/bin/uncrustify" CACHE PATH "")

# sparsehash headers from uberenv
set(SPARSEHASH_DIR "/home/taylor16/tpl/v2" CACHE PATH "")

# boost headers from uberenv
set(BOOST_DIR "/home/taylor16/tpl/v2" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "/home/taylor16/tpl/v2/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "/home/taylor16/tpl/v2/usr/bin/genhtml" CACHE PATH "")

##################################
# end uberenv host-config
##################################

#######
# MPI - manually added these for now.
#######
#set(ENABLE_MPI ON CACHE PATH "")
#set(MPI_C_COMPILER "/home/taylor16/local/mpich-3.2/bin/mpicc" CACHE PATH "")
#set(MPI_CXX_COMPILER "/home/taylor16/local/mpich-3.2/bin/mpicxx" CACHE PATH "")
#set(MPI_Fortran_COMPILER "/home/taylor16/local/mpich-3.2/bin/mpif90" CACHE PATH "")

##############################################################################
# SHROUD - manually added for now. Use a public build add to TPL later
##############################################################################
set(SHROUD_EXECUTABLE "/home/taylor16/tpl/shroud/bin/shroud" CACHE PATH "")

