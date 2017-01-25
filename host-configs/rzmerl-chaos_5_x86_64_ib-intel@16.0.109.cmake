##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# chaos_5_x86_64_ib-intel@16.0.109
##################################

# cmake from uberenv
# cmake executable path: /usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_12_28_16_37_03/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.109/cmake-3.3.1-gbdavlh2pksdcwqwprmr6ectl2zu7ztv/bin/cmake

#######
# using intel@16.0.109 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/local/tools/ic-16.0.109/bin/icc" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/local/tools/ic-16.0.109/bin/icpc" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE PATH "")

set(CMAKE_Fortran_COMPILER  "/usr/local/tools/ic-16.0.109/bin/ifort" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_12_28_16_37_03/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.109/hdf5-1.8.16-uvurt4h25puxc4woxrgbt6sr3gaxropi" CACHE PATH "")

# conduit from uberenv
set(CONDUIT_DIR "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_12_28_16_37_03/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.109/conduit-0.2.0-sjnz677tsmmc4ohlovin5gtgtwlif6iq" CACHE PATH "")

# doxygen from uberenv
set(DOXYGEN_EXECUTABLE "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_12_28_16_37_03/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.109/doxygen-1.8.11-qafjsbjfdlu52j3ozltk6ddlti235pe6/bin/doxygen" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_12_28_16_37_03/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.109/python-2.7.11-3ylay5hhgxrjdh53fafxesjsh4tllmwh/bin/python" CACHE PATH "")

# lua from uberenv
set(LUA_DIR "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_12_28_16_37_03/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.109/lua-5.1.5-25nzqxaaljfyf6eixkd2bczlmsmpxplv" CACHE PATH "")

# sphinx from uberenv
set(SPHINX_EXECUTABLE "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_12_28_16_37_03/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.109/python-2.7.11-3ylay5hhgxrjdh53fafxesjsh4tllmwh/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_12_28_16_37_03/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.109/uncrustify-0.61-7v4io4dwzqwtplzuykvomoxpibd3cnbl/bin/uncrustify" CACHE PATH "")

# sparsehash headers from uberenv
set(SPARSEHASH_DIR "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_12_28_16_37_03/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.109/sparsehash-headers-2.0.2-cv3sza7qef4mtf4lthelmarncfaz7ygl" CACHE PATH "")

# boost headers from uberenv
set(ENABLE_BOOST ON CACHE PATH "")
set(BOOST_ROOT "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_12_28_16_37_03/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.109/boost-headers-1.58.0-vaxfagykvpueronnzybsyr63rhcmldzj" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_12_28_16_37_03/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.109/lcov-1.11-nrsx4h6lyz4ajjrqi45dphsnbppj3da3/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_12_28_16_37_03/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.109/lcov-1.11-nrsx4h6lyz4ajjrqi45dphsnbppj3da3/usr/bin/genhtml" CACHE PATH "")

# Temporarily disable CXX11 on intel builds until we resolve issue ATK-619
set(BLT_CXX_STD "c++98" CACHE PATH "")
##################################
# end uberenv host-config
##################################

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################
# Options added manually to 
# lc chaos5 intel@16.0.109 host configs
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
# !---------------------------------------------------------------------------
##############################################################################

