##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# chaos_5_x86_64_ib-intel@15.0.187
##################################

# cmake from uberenv
# cmake executable path: /usr/workspace/wsa/axom/thirdparty_libs/builds/2017_04_24_15_57_37/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.187/cmake-3.3.1-2ffcsc32b75exq6yu3baa7n7yk5rnael/bin/cmake

#######
# using intel@15.0.187 compiler spec
#######

# c compiler used by spack
set("CMAKE_C_COMPILER" "/usr/local/tools/ic-15.0.187/bin/icc" CACHE PATH "")

# cpp compiler used by spack
set("CMAKE_CXX_COMPILER" "/usr/local/tools/ic-15.0.187/bin/icpc" CACHE PATH "")

# fortran compiler used by spack
set("ENABLE_FORTRAN" "ON" CACHE PATH "")

set("CMAKE_Fortran_COMPILER" "/usr/local/tools/ic-15.0.187/bin/ifort" CACHE PATH "")

# hdf5 from uberenv
set("HDF5_DIR" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_04_24_15_57_37/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.187/hdf5-1.8.16-dvlk5aqfby3atjbqefigk6ugby5dk5yy" CACHE PATH "")

# conduit from uberenv
set("CONDUIT_DIR" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_04_24_15_57_37/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.187/conduit-0.2.1-ldvfo5r46mpnc2taxmwm6lt2sa53qomi" CACHE PATH "")

# boost headers from uberenv
set("BOOST_DIR" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_04_24_15_57_37/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.187/boost-headers-1.58.0-m4zq524o33mi34qimoxtn5rae4amyk6w" CACHE PATH "")

# python from uberenv
set("PYTHON_EXECUTABLE" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_04_24_15_57_37/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.187/python-2.7.11-2us32xlalxdytptzzhfiniiimqrf4ezh/bin/python" CACHE PATH "")

# lua from uberenv
set("LUA_DIR" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_04_24_15_57_37/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.187/lua-5.1.5-7gpjuo5khxsxqwidz74wccitsl3hdixj" CACHE PATH "")

# doxygen from uberenv
set("DOXYGEN_EXECUTABLE" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_04_24_15_57_37/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.187/doxygen-1.8.11-hnig4yzjrv6izpmdpstyr6bxyupnzhpu/bin/doxygen" CACHE PATH "")

# sphinx from uberenv
set("SPHINX_EXECUTABLE" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_04_24_15_57_37/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.187/python-2.7.11-2us32xlalxdytptzzhfiniiimqrf4ezh/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
set("UNCRUSTIFY_EXECUTABLE" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_04_24_15_57_37/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.187/uncrustify-0.61-onqgsfbevtniuj45s564dpfkmlveql33/bin/uncrustify" CACHE PATH "")

# lcov and genhtml from uberenv
set("LCOV_PATH" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_04_24_15_57_37/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.187/lcov-1.11-hyyvo3em3ckyvwg5bz27koz2ulmsbms3/usr/bin/lcov" CACHE PATH "")

set("GENHTML_PATH" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_04_24_15_57_37/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.187/lcov-1.11-hyyvo3em3ckyvwg5bz27koz2ulmsbms3/usr/bin/genhtml" CACHE PATH "")

# Temporarily disable CXX11 on intel builds until we resolve issue ATK-619
set("BLT_CXX_STD" "c++98" CACHE PATH "")

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

set(MPIEXEC "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

##############################################################################
# SHROUD - manually added for now. Use a public build add to TPL later
##############################################################################
set(SHROUD_EXECUTABLE "/usr/apps/shroud/bin/shroud" CACHE PATH "")

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

