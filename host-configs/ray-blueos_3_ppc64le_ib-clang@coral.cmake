##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# blueos_3_ppc64le_ib-clang@coral
##################################

# cmake from uberenv
# cmake executable path: /usr/workspace/wsa/axom/thirdparty_libs/builds/2017_11_13_15_16_15/spack/opt/spack/blueos_3_ppc64le_ib/clang-coral/cmake-3.8.2-ze553yn73s3b32uhop7jlgvjfxo4tjuc/bin/cmake

#######
# using clang@coral compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/tcetmp/packages/clang/clang-coral-2017.05.19/bin/clang" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/tcetmp/packages/clang/clang-coral-2017.05.19/bin/clang++" CACHE PATH "")

# fortran compiler used by spack
# no fortran compiler

set(ENABLE_FORTRAN OFF CACHE BOOL "")

# Root directory for generated TPLs
set(TPL_ROOT "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_11_13_15_16_15/spack/opt/spack/blueos_3_ppc64le_ib/clang-coral" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.16-jybmvmtpcohhrq4trij72mbwj6blad4l" CACHE PATH "")

# scr not built by uberenv

# conduit from uberenv
set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.2.1-yjoop56c7fektxu4afrq6xtmz27ujtzz" CACHE PATH "")

# mfem from uberenv
set(MFEM_DIR "${TPL_ROOT}/mfem-3.3.2-vrr2ba7mdlbbywemnp33tzgsyqj5ldje" CACHE PATH "")

# boost headers from uberenv
set(BOOST_DIR "${TPL_ROOT}/boost-headers-1.58.0-h5ezasjj2b4jbmznzonrpayzubdzb4fl" CACHE PATH "")

# python not built by uberenv

# lua not built by uberenv

# doxygen not built by uberenv

# sphinx not built by uberenv

# shroud not built by uberenv

# uncrustify not built by uberenv

# lcov and genhtml not built by uberenv

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
# lc blueos clang@coral  host configs
##############################################################################

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

##############################################################################
# MPI - manually added for now
##############################################################################

set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_HOME                 "/usr/tcetmp/packages/spectrum-mpi/spectrum-mpi-2017.04.03-clang-coral-2017.05.19" CACHE PATH "")
set(MPI_C_COMPILER           "${MPI_HOME}/bin/mpicc"   CACHE PATH "")
set(MPI_CXX_COMPILER         "${MPI_HOME}/bin/mpicxx"  CACHE PATH "")
set(MPI_Fortran_COMPILER     "${MPI_HOME}/bin/mpifort" CACHE PATH "")

set(MPI_DRIVER_ROOT          "/usr/tcetmp/packages/xlflang/xlflang-2016-11-30/xlf" CACHE PATH "")

set(MPI_LIBS 
    ${MPI_DRIVER_ROOT}/lib/libxl.a
    ${MPI_DRIVER_ROOT}/lib/libxlcuf.a
    ${MPI_DRIVER_ROOT}/lib/libxlcuf_nonshr.a
    ${MPI_DRIVER_ROOT}/lib/libxlcufdevice.a
    ${MPI_DRIVER_ROOT}/lib/libxldevice.a
    ${MPI_DRIVER_ROOT}/lib/libxlmbif.a
    ${MPI_DRIVER_ROOT}/lib/libxlopt.a)

set(MPI_INCLUDE_PATHS 
    ${MPI_DRIVER_ROOT}/include)

set(MPI_Fortran_LIBRARIES    "${MPI_LIBS}" CACHE PATH "")
set(MPI_Fortran_INCLUDE_PATH "${MPI_INCLUDE_PATHS}" CACHE PATH "")

set(MPIEXEC              "mpirun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-np" CACHE PATH "")


##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

