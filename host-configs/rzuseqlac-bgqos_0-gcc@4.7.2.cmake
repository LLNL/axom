##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# bgqos_0-gcc@4.7.2
##################################

# cmake from uberenv
# cmake executable path: /collab/usr/global/tools/cmake/bgqos_0/cmake-3.1.2/bin/cmake

#######
# using gcc@4.7.2 compiler spec
#######

# c compiler used by spack
set("CMAKE_C_COMPILER" "/usr/local/tools/toolchain-4.7.2/scripts/bggcc" CACHE PATH "")

# cpp compiler used by spack
set("CMAKE_CXX_COMPILER" "/usr/local/tools/toolchain-4.7.2/scripts/bgg++" CACHE PATH "")

# fortran compiler used by spack
set("ENABLE_FORTRAN" "OFF" CACHE PATH "")

set("CMAKE_Fortran_COMPILER" "/usr/local/tools/toolchain-4.7.2/scripts/bggfortran" CACHE PATH "")

# hdf5 from uberenv
set("HDF5_DIR" "/usr/workspace/wsrzc/axom/thirdparty_libs/builds/2017_04_17_17_46_38/spack/opt/spack/bgqos_0/gcc-4.7.2/hdf5-1.8.16-gsqekrbdiryaddbxrb7lkv5ivwyyylcc" CACHE PATH "")

# conduit from uberenv
set("CONDUIT_DIR" "/usr/workspace/wsrzc/axom/thirdparty_libs/builds/2017_04_17_17_46_38/spack/opt/spack/bgqos_0/gcc-4.7.2/conduit-0.2.1-y7kksgazxhxo4qkgbz54vdjfkpmkqu52" CACHE PATH "")

# boost headers from uberenv
set("BOOST_DIR" "/usr/workspace/wsrzc/axom/thirdparty_libs/builds/2017_04_17_17_46_38/spack/opt/spack/bgqos_0/gcc-4.7.2/boost-headers-1.58.0-bgbbwvnccff6bzg32vt2dk5343c346jn" CACHE PATH "")

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
# lc bgq gcc@4.7.2 host configs
##############################################################################

set(ENABLE_DOCS OFF CACHE PATH "")

##############################################################################
# MPI - manually added for now
##############################################################################
set(ENABLE_MPI ON CACHE PATH "")
set(MPI_C_COMPILER "/usr/local/tools/compilers/ibm/mpicc-4.7.2-fastmpi" CACHE PATH "")
set(MPI_CXX_COMPILER "/usr/local/tools/compilers/ibm/mpicxx-4.7.2-fastmpi" CACHE PATH "")
set(MPI_Fortran_COMPILER  "/usr/local/tools/compilers/ibm/mpigfortran-4.7.2-fastmpi" CACHE PATH "")

set(MPI_LIBS "/bgsys/drivers/V1R2M4/ppc64/comm/lib/libmpich-gcc.a;/bgsys/drivers/V1R2M4/ppc64/comm/lib/libopa-gcc.a;/bgsys/drivers/V1R2M4/ppc64/comm/lib/libmpl-gcc.a;/bgsys/drivers/V1R2M4/ppc64/comm/lib/libpami-gcc.a;/bgsys/drivers/V1R2M4/ppc64/spi/lib/libSPI.a;/bgsys/drivers/V1R2M4/ppc64/spi/lib/libSPI_cnk.a;rt;pthread;stdc++;pthread")

set(MPI_INCLUDE_PATHS "/bgsys/drivers/V1R2M4/ppc64/comm/include;/bgsys/drivers/V1R2M4/ppc64/comm/lib/gnu;/bgsys/drivers/V1R2M4/ppc64;/bgsys/drivers/V1R2M4/ppc64/comm/sys/include;/bgsys/drivers/V1R2M4/ppc64/spi/include;/bgsys/drivers/V1R2M4/ppc64/spi/include/kernel/cnk" )

set(MPI_C_INCLUDE_PATH ${MPI_INCLUDE_PATHS} CACHE PATH "")
set(MPI_C_LIBRARIES ${MPI_LIBS} CACHE PATH "")

set(MPI_CXX_INCLUDE_PATH  ${MPI_INCLUDE_PATHS} CACHE PATH "")
set(MPI_CXX_LIBRARIES ${MPI_LIBS} CACHE PATH "")
set(MPI_Fortran_LIBRARIES ${MPI_LIBS} CACHE PATH "")


set(MPIEXEC "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")


# GTest death tests use forked threads, which does now work on BG/Q 
set(EXTRA_C_FLAGS   -DGTEST_HAS_DEATH_TEST=0 CACHE PATH "")
set(EXTRA_CXX_FLAGS -DGTEST_HAS_DEATH_TEST=0 CACHE PATH "")

set(BLT_ALWAYS_WRAP_TESTS_WITH_MPIEXEC TRUE CACHE PATH "Ensures that tests will be wrapped with srun to run on the backend nodes")

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

