###############################################################################
#
# CMake Cache Seed file for bgqos_0 machines using gcc/g++ 4.7.2
#
###############################################################################

#######
# uberenv host-config for asctoolkit
#######
# cmake from uberenv
# cmake executable path: /usr/workspace/wsa/toolkit/thirdparty_libs/builds/mirrorbgq/cmake/cmake_3.3.1/cmake-3.3.1/bin/cmake 
###############################################################################

# Select the c and c++ compiler though the standard CMake Variables.
###############################################################################
set(CMAKE_C_COMPILER "/usr/local/tools/compilers/ibm/mpicc-4.7.2" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/local/tools/compilers/ibm/mpicxx-4.7.2" CACHE PATH "")
set(CMAKE_Fortran_COMPILER "/usr/local/tools/compilers/ibm/mpif90-4.7.2" CACHE PATH "")


# python from uberenv
set(PYTHON_EXECUTABLE "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/mirrorbgq/python/python-2.7.1/Python-2.7.11/python" CACHE PATH "")

###############################################################################
# Additional Compiler Flags
###############################################################################

# additional flags for all configurations 
#set(EXTRA_CXX_FLAGS "-DEXTRA_FLAGS_CXX_DEFINE" CACHE PATH "")

# additional flags for debug builds
# set(EXTRA_CXX_DEBUG_FLAGS "-DEXTRA_CXX_DEBUG_FLAGS_DEFINE" CACHE PATH "")

# additional flags for debug builds
# set(EXTRA_CXX_RELEASE_FLAGS "-DEXTRA_CXX_RELEASE_FLAGS_DEFINE" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/mirrorbgq/hdf5/hdf5-1.8.16" CACHE PATH "")

# conduit from uberenv
set(CONDUIT_DIR "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/mirrorbgq/conduit/conduit" CACHE PATH "")

##################################
# end uberenv host-config
##################################

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################
# Options added manually to 
# lc bgq gcc@4.7.2  host configs
##############################################################################

##############################################################################
# MPI - manually added for now
##############################################################################
set(ENABLE_MPI ON CACHE PATH "")
set(MPI_C_COMPILER "/usr/local/tools/compilers/ibm/mpicxx-4.7.2" CACHE PATH "")
set(MPI_CXX_COMPILER "/usr/local/tools/compilers/ibm/mpicxx-4.7.2" CACHE PATH "")
set(MPI_Fortran_COMPILER "/usr/local/tools/compilers/ibm/mpif90-4.7.2" CACHE PATH "")

##############################################################################
# GCOV - manually added for now
##############################################################################
#set(GCOV_PATH "/usr/apps/gnu/4.9.3/bin/gcov" CACHE PATH "")

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################
