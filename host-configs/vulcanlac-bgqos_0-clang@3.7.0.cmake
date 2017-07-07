##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# bgqos_0-clang@3.7.0
##################################

# cmake from uberenv
# cmake executable path: /usr/global/tools/CMake/bgqos_0/cmake-3.1.2/bin/cmake

#######
# using clang@3.7.0 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/local/bin/bgclang" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/local/bin/bgclang++" CACHE PATH "")

# fortran compiler used by spack
# no fortran compiler

set(ENABLE_FORTRAN OFF CACHE BOOL "")

# Root directory for generated TPLs
set(TPL_ROOT "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_05_15_17_28_06/spack/opt/spack/bgqos_0/clang-3.7.0" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.16-bosaqxj3xd5fhyovqnda3rgj2kjsj4ah" CACHE PATH "")

# conduit from uberenv
set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.2.1-oiiieme5mlcpao7pqwrk2mdquxnaguqm" CACHE PATH "")

# mfem from uberenv
set(MFEM_DIR "${TPL_ROOT}/mfem-3.3-2ctejypbka5bdvs43cnv3twskptqpfjs" CACHE PATH "")

# boost headers from uberenv
set(BOOST_DIR "${TPL_ROOT}/boost-headers-1.58.0-qddl3bajxtossmhy4mazvjpah4zgx5aj" CACHE PATH "")

# python not build by uberenv

# lua not build by uberenv

# doxygen not built by uberenv

# sphinx not built by uberenv

# uncrustify not built by uberenv

# lcov and genhtml from uberenv
# lcov and genhtml not built by uberenv

##################################
# end uberenv host-config
##################################

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################
# Options added manually to 
# lc bgq clang@3.7.0 host configs
##############################################################################

set(ENABLE_DOCS OFF CACHE PATH "")

##############################################################################
# MPI - manually added for now
##############################################################################
set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_HOME             "/usr/apps/gnu/clang/r266321-20160414/mpi/bgclang-mpi3" CACHE PATH "")
set(MPI_C_COMPILER       "${MPI_HOME}/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpicxx" CACHE PATH "")


set(MPI_DRIVER_ROOT      "/bgsys/drivers/V1R2M4/ppc64" CACHE PATH "")

set(MPI_LIBS 
    ${MPI_DRIVER_ROOT}/comm/lib/libmpich-gcc.a
    ${MPI_DRIVER_ROOT}/comm/lib/libopa-gcc.a
    ${MPI_DRIVER_ROOT}/comm/lib/libmpl-gcc.a
    ${MPI_DRIVER_ROOT}/comm/lib/libpami-gcc.a
    ${MPI_DRIVER_ROOT}/spi/lib/libSPI.a
    ${MPI_DRIVER_ROOT}/spi/lib/libSPI_cnk.a
    rt
    pthread
    stdc++
    pthread)

set(MPI_INCLUDE_PATHS 
    ${MPI_DRIVER_ROOT}/comm/include
    ${MPI_DRIVER_ROOT}/comm/lib/gnu
    ${MPI_DRIVER_ROOT}
    ${MPI_DRIVER_ROOT}/comm/sys/include
    ${MPI_DRIVER_ROOT}/spi/include
    ${MPI_DRIVER_ROOT}/spi/include/kernel/cnk )

set(MPI_C_INCLUDE_PATH    ${MPI_INCLUDE_PATHS} CACHE PATH "")
set(MPI_C_LIBRARIES       ${MPI_LIBS} CACHE PATH "")

set(MPI_CXX_INCLUDE_PATH  ${MPI_INCLUDE_PATHS} CACHE PATH "")
set(MPI_CXX_LIBRARIES     ${MPI_LIBS} CACHE PATH "")


set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")


# GTest death tests use forked threads, which does now work on BG/Q 
set(EXTRA_C_FLAGS   -DGTEST_HAS_DEATH_TEST=0 CACHE PATH "")
set(EXTRA_CXX_FLAGS -DGTEST_HAS_DEATH_TEST=0 CACHE PATH "")

set(BLT_ALWAYS_WRAP_TESTS_WITH_MPIEXEC TRUE CACHE BOOL "Ensures that tests will be wrapped with srun to run on the backend nodes")

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

