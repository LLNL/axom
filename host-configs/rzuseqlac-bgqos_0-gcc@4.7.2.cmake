##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# bgqos_0-gcc@4.7.2
##################################

# cmake from uberenv
# cmake executable path: /usr/local/tools/cmake-3.4.3/bin/cmake

#######
# using gcc@4.7.2 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/local/tools/toolchain-4.7.2/scripts/bggcc" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/local/tools/toolchain-4.7.2/scripts/bgg++" CACHE PATH "")

# fortran compiler used by spack
# no fortran compiler

set(ENABLE_FORTRAN OFF CACHE BOOL "")

# Root directory for generated TPLs
set(TPL_ROOT "/usr/workspace/wsrzc/axom/thirdparty_libs/builds/2017_07_19_14_04_28/spack/opt/spack/bgqos_0/gcc-4.7.2" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.16-gsqekrbdiryaddbxrb7lkv5ivwyyylcc" CACHE PATH "")

# conduit from uberenv
set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.2.1-y7kksgazxhxo4qkgbz54vdjfkpmkqu52" CACHE PATH "")

# mfem from uberenv
set(MFEM_DIR "${TPL_ROOT}/mfem-3.3-wb2a3z3lptwlnbbnsh4ssyz7zwa24ekj" CACHE PATH "")

# boost headers from uberenv
set(BOOST_DIR "${TPL_ROOT}/boost-headers-1.58.0-bgbbwvnccff6bzg32vt2dk5343c346jn" CACHE PATH "")

# python not built by uberenv

# lua not built by uberenv

# doxygen not built by uberenv

# sphinx not built by uberenv

# uncrustify not built by uberenv

# lcov and genhtml not built by uberenv

##################################
# end uberenv host-config
##################################

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################
# Options added manually to 
# lc bgq gcc@4.7.2 host configs
##############################################################################

set(ENABLE_DOCS    OFF CACHE BOOL "")
set(ENABLE_PYTHON  OFF CACHE BOOL "")

set(CMAKE_SKIP_RPATH TRUE CACHE BOOL "")

##############################################################################
# MPI - manually added for now
##############################################################################
set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_HOME             "/usr/local/tools/compilers/ibm" CACHE PATH "")
set(MPI_C_COMPILER       "${MPI_HOME}/mpicc-4.7.2-fastmpi" CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_HOME}/mpicxx-4.7.2-fastmpi" CACHE PATH "")
set(MPI_Fortran_COMPILER "${MPI_HOME}/mpigfortran-4.7.2-fastmpi" CACHE PATH "")

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

set(MPI_Fortran_LIBRARIES ${MPI_LIBS} CACHE PATH "")


set(MPIEXEC               "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG  "-n" CACHE PATH "")

set(ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC TRUE CACHE BOOL "Ensures that tests will be wrapped with srun to run on the backend nodes")

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

