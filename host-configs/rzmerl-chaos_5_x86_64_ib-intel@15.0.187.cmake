##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# chaos_5_x86_64_ib-intel@15.0.187
##################################

# cmake from uberenv
# cmake executable path: /usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_11_10_18_58_00/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.187/cmake-3.3.1-v6a26kd37eqvxt7dcmydhfc4kynkpnhx/bin/cmake

#######
# using intel@15.0.187 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/local/tools/ic-15.0.187/bin/icc" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/local/tools/ic-15.0.187/bin/icpc" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE PATH "")

set(CMAKE_Fortran_COMPILER  "/usr/local/tools/ic-15.0.187/bin/ifort" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_11_10_18_58_00/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.187/hdf5-1.8.16-ibq54wdkvnyrbkxedqk5nou77k6aiyet" CACHE PATH "")

# conduit from uberenv
set(CONDUIT_DIR "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_11_10_18_58_00/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.187/conduit-0.2.0-yd5te3j4pykt5karfj7koopc2sa5fxf4" CACHE PATH "")

# doxygen from uberenv
set(DOXYGEN_EXECUTABLE "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_11_10_18_58_00/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.187/doxygen-1.8.11-f5yxbpnswwocut2dsstpwfuywh2ashak/bin/doxygen" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_11_10_18_58_00/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.187/python-2.7.11-2us32xlalxdytptzzhfiniiimqrf4ezh/bin/python" CACHE PATH "")

# lua from uberenv
set(LUA_DIR "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_11_10_18_58_00/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.187/lua-5.1.5-7gpjuo5khxsxqwidz74wccitsl3hdixj" CACHE PATH "")

# sphinx from uberenv
set(SPHINX_EXECUTABLE "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_11_10_18_58_00/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.187/python-2.7.11-2us32xlalxdytptzzhfiniiimqrf4ezh/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_11_10_18_58_00/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.187/uncrustify-0.61-onqgsfbevtniuj45s564dpfkmlveql33/bin/uncrustify" CACHE PATH "")

# sparsehash headers from uberenv
set(SPARSEHASH_DIR "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_11_10_18_58_00/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.187/sparsehash-headers-2.0.2-4q2ms4ib3yes7y7hrekj5zigrrzh3qhp" CACHE PATH "")

# boost headers from uberenv
set(ENABLE_BOOST ON CACHE PATH "")
set(BOOST_ROOT "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_11_10_18_58_00/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.187/boost-headers-1.58.0-m4zq524o33mi34qimoxtn5rae4amyk6w" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_11_10_18_58_00/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.187/lcov-1.11-hyyvo3em3ckyvwg5bz27koz2ulmsbms3/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_11_10_18_58_00/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.187/lcov-1.11-hyyvo3em3ckyvwg5bz27koz2ulmsbms3/usr/bin/genhtml" CACHE PATH "")

# Temporarily disable CXX11 on intel builds until we resolve issue ATK-619
set(BLT_CXX_STD "c++98" CACHE PATH "")
##################################
# end uberenv host-config
##################################

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################
# Options added manually to 
# lc chaos5 intel@15.0.187 host configs
##############################################################################

##############################################################################
# MPI - manually added for now
##############################################################################
set(ENABLE_MPI ON CACHE PATH "")
set(MPI_C_COMPILER "/usr/local/tools/mvapich2-intel-2.0/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "/usr/local/tools/mvapich2-intel-2.0/bin/mpicxx" CACHE PATH "")
set(MPI_Fortran_COMPILER  "/usr/local/tools/mvapich2-intel-2.0/bin/mpif90" CACHE PATH "")

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

