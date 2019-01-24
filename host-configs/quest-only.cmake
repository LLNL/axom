##################################
#
# Copyright (c) 2019, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory.
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
#
##################################
##############################################################################
# Minimal host-config (cmake cache file) to build quest
#
# This host-config disables the sidre component
# and does not require any pre-built third party libraries
# 
# To build the code with this host-config, 
# run the following from your build directory: 
#
#   cmake -C quest-only-toss3-gcc@8.1.0.cmake           \
#         -DCMAKE_BUILD_TYPE={Debug,Release}            \
#         -DENABLE_EXAMPLES={ON,OFF}                    \
#         -DENABLE_TESTS={ON,OFF}                       \
#         -DCMAKE_INSTALL_PREFIX= /path/to/install/dir  \
#         <axom_root>/src
#
##############################################################################
# toss_3_x86_64_ib-gcc@8.1.0
#
##############################################################################

# cmake executable path: /usr/tce/packages/cmake/cmake-3.9.2/bin/cmake

set(ENABLE_FORTRAN ON CACHE BOOL "")

#set(COMPILER_HOME "/usr/tce/packages/gcc/gcc-8.1.0" )
#set(CMAKE_C_COMPILER "${COMPILER_HOME}/bin/gcc" CACHE PATH "")
#set(CMAKE_CXX_COMPILER "${COMPILER_HOME}/bin/g++" CACHE PATH "")
#set(CMAKE_Fortran_COMPILER "${COMPILER_HOME}/bin/gfortran" CACHE PATH "")

set(COMPILER_HOME "/usr/tce/packages/intel/intel-18.0.2" )
set(CMAKE_C_COMPILER "${COMPILER_HOME}/bin/icc" CACHE PATH "")
set(CMAKE_CXX_COMPILER "${COMPILER_HOME}/bin/icpc" CACHE PATH "")
set(CMAKE_Fortran_COMPILER "${COMPILER_HOME}/bin/ifort" CACHE PATH "")

# Make shared lib
set(BUILD_SHARED_LIBS "ON" CACHE BOOL "")

#set(SHROUD_EXECUTABLE "/usr/workspace/wsrzd/olson45/SAND/pyranda/env_genie/bin/shroud" CACHE PATH "")
set(SHROUD_EXECUTABLE "/g/g14/taylor/shroud/shroud-quest/build/temp.linux-x86_64-2.7/venv/bin/shroud" CACHE PATH "")


# Sidre requires conduit and hdf5, so disable it in this host-config
set(AXOM_ENABLE_SIDRE OFF CACHE BOOL "")

##############################################################################
# MPI
##############################################################################
set(ENABLE_MPI ON CACHE BOOL "")

#set(MPI_HOME             "/usr/tce/packages/mvapich2/mvapich2-2.2-gcc-8.1.0" CACHE PATH "")
set(MPI_HOME             "/usr/tce/packages/mvapich2/mvapich2-2.2-intel-18.0.2" CACHE PATH "")
set(MPI_C_COMPILER       "${MPI_HOME}/bin/mpicc"   CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpicxx"  CACHE PATH "")
set(MPI_Fortran_COMPILER "${MPI_HOME}/bin/mpif90" CACHE PATH "")

set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

