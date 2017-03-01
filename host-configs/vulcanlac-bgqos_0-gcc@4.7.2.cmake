##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# bgqos_0-gcc@4.7.2
##################################

# cmake from uberenv
# cmake executable path: /usr/global/tools/CMake/bgqos_0/cmake-3.1.2/bin/cmake

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
set("HDF5_DIR" "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_01_13_40_58/spack/opt/spack/bgqos_0/gcc-4.7.2/hdf5-1.8.16-gsqekrbdiryaddbxrb7lkv5ivwyyylcc" CACHE PATH "")

# conduit from uberenv
set("CONDUIT_DIR" "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_01_13_40_58/spack/opt/spack/bgqos_0/gcc-4.7.2/conduit-0.2.1-y7kksgazxhxo4qkgbz54vdjfkpmkqu52" CACHE PATH "")

# sparsehash headers from uberenv
set("SPARSEHASH_DIR" "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_01_13_40_58/spack/opt/spack/bgqos_0/gcc-4.7.2/sparsehash-headers-2.0.2-24icetuqyud64fxt3dfx3vxlucy5apk6" CACHE PATH "")

# boost headers from uberenv
set("BOOST_DIR" "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_01_13_40_58/spack/opt/spack/bgqos_0/gcc-4.7.2/boost-headers-1.58.0-bgbbwvnccff6bzg32vt2dk5343c346jn" CACHE PATH "")

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

##############################################################################
# MPI - manually added for now
##############################################################################
set(ENABLE_MPI ON CACHE PATH "")
set(MPI_C_COMPILER "/usr/local/tools/compilers/ibm/mpicc-4.7.2-fastmpi" CACHE PATH "")
set(MPI_CXX_COMPILER "/usr/local/tools/compilers/ibm/mpicxx-4.7.2-fastmpi" CACHE PATH "")
set(MPI_Fortran_COMPILER  "/usr/local/tools/compilers/ibm/mpif90-4.7.2-fastmpi" CACHE PATH "")

set(MPIEXEC "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

