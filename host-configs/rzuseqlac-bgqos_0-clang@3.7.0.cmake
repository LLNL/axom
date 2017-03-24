##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# bgqos_0-clang@3.7.0
##################################

# cmake from uberenv
# cmake executable path: /collab/usr/global/tools/cmake/bgqos_0/cmake-3.1.2/bin/cmake

#######
# using clang@3.7.0 compiler spec
#######

# c compiler used by spack
set("CMAKE_C_COMPILER" "/usr/local/bin/bgclang" CACHE PATH "")

# cpp compiler used by spack
set("CMAKE_CXX_COMPILER" "/usr/local/bin/bgclang++" CACHE PATH "")

# fortran compiler used by spack
# no fortran compiler

set("ENABLE_FORTRAN" "OFF" CACHE PATH "")

# hdf5 from uberenv
set("HDF5_DIR" "/usr/workspace/wsrzc/axom/thirdparty_libs/builds/2017_03_17_12_33_58/spack/opt/spack/bgqos_0/clang-3.7.0/hdf5-1.8.16-bosaqxj3xd5fhyovqnda3rgj2kjsj4ah" CACHE PATH "")

# conduit from uberenv
set("CONDUIT_DIR" "/usr/workspace/wsrzc/axom/thirdparty_libs/builds/2017_03_17_12_33_58/spack/opt/spack/bgqos_0/clang-3.7.0/conduit-0.2.1-oiiieme5mlcpao7pqwrk2mdquxnaguqm" CACHE PATH "")

# sparsehash headers from uberenv
set("SPARSEHASH_DIR" "/usr/workspace/wsrzc/axom/thirdparty_libs/builds/2017_03_17_12_33_58/spack/opt/spack/bgqos_0/clang-3.7.0/sparsehash-headers-2.0.2-jnxnoo3nsy2l6vutmvlvwri4tzgbqmro" CACHE PATH "")

# boost headers from uberenv
set("BOOST_DIR" "/usr/workspace/wsrzc/axom/thirdparty_libs/builds/2017_03_17_12_33_58/spack/opt/spack/bgqos_0/clang-3.7.0/boost-headers-1.58.0-qddl3bajxtossmhy4mazvjpah4zgx5aj" CACHE PATH "")

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

##############################################################################
# MPI - manually added for now
##############################################################################
set(ENABLE_MPI ON CACHE PATH "")
#set(MPI_C_COMPILER "/usr/apps/gnu/clang/llnl/bin/mpiclang-fastmpi" CACHE PATH "")
#set(MPI_CXX_COMPILER "/usr/apps/gnu/clang/llnl/bin/mpiclang++11-fastmpi" CACHE PATH "")

set(MPI_C_COMPILER "/usr/apps/gnu/clang/r266321-20160414/mpi/bgclang-mpi3/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "/usr/apps/gnu/clang/r266321-20160414/mpi/bgclang-mpi3/bin/mpicxx" CACHE PATH "")

set(MPI_LIBS "/bgsys/drivers/V1R2M4/ppc64/comm/lib/libmpich-gcc.a;/bgsys/drivers/V1R2M4/ppc64/comm/lib/libopa-gcc.a;/bgsys/drivers/V1R2M4/ppc64/comm/lib/libmpl-gcc.a;/bgsys/drivers/V1R2M4/ppc64/comm/lib/libpami-gcc.a;/bgsys/drivers/V1R2M4/ppc64/spi/lib/libSPI.a;/bgsys/drivers/V1R2M4/ppc64/spi/lib/libSPI_cnk.a;rt;pthread;stdc++;pthread")

set(MPI_INCLUDE_PATHS "/bgsys/drivers/V1R2M4/ppc64/comm/include;/bgsys/drivers/V1R2M4/ppc64/comm/lib/gnu;/bgsys/drivers/V1R2M4/ppc64;/bgsys/drivers/V1R2M4/ppc64/comm/sys/include;/bgsys/drivers/V1R2M4/ppc64/spi/include;/bgsys/drivers/V1R2M4/ppc64/spi/include/kernel/cnk" )

set(MPI_C_INCLUDE_PATH ${MPI_INCLUDE_PATHS} CACHE PATH "")
set(MPI_C_LIBRARIES ${MPI_LIBS} CACHE PATH "")

set(MPI_CXX_INCLUDE_PATH  ${MPI_INCLUDE_PATHS} CACHE PATH "")
set(MPI_CXX_LIBRARIES ${MPI_LIBS} CACHE PATH "")


set(MPIEXEC "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")


# GTest death tests use forked threads, which does now work on BG/Q 
set(EXTRA_C_FLAGS   -DGTEST_HAS_DEATH_TEST=0 CACHE PATH "")
set(EXTRA_CXX_FLAGS -DGTEST_HAS_DEATH_TEST=0 CACHE PATH "")


##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

