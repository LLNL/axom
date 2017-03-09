##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# toss_3_x86_64_ib-intel@16.0.4
##################################

# cmake from uberenv
# cmake executable path: /usr/workspace/wsa/axom/thirdparty_libs/builds/2017_03_08_13_59_54/spack/opt/spack/toss_3_x86_64_ib/intel-16.0.4/cmake-3.3.1-pameyoqbnyiwmdxbeqwmiup7kswibf35/bin/cmake

#######
# using intel@16.0.4 compiler spec
#######

# c compiler used by spack
set("CMAKE_C_COMPILER" "/usr/tce/packages/intel/intel-16.0.4/bin/icc" CACHE PATH "")

# cpp compiler used by spack
set("CMAKE_CXX_COMPILER" "/usr/tce/packages/intel/intel-16.0.4/bin/icpc" CACHE PATH "")

# fortran compiler used by spack
set("ENABLE_FORTRAN" "ON" CACHE PATH "")

set("CMAKE_Fortran_COMPILER" "/usr/tce/packages/intel/intel-16.0.4/bin/ifort" CACHE PATH "")

# hdf5 from uberenv
set("HDF5_DIR" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_03_08_13_59_54/spack/opt/spack/toss_3_x86_64_ib/intel-16.0.4/hdf5-1.8.16-z6nt2tyk2j2eop2df6ghod6dnfhwei6s" CACHE PATH "")

# conduit from uberenv
set("CONDUIT_DIR" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_03_08_13_59_54/spack/opt/spack/toss_3_x86_64_ib/intel-16.0.4/conduit-0.2.1-o2svxknpl6alrrdvuislszxqhljo4k3b" CACHE PATH "")

# sparsehash headers from uberenv
set("SPARSEHASH_DIR" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_03_08_13_59_54/spack/opt/spack/toss_3_x86_64_ib/intel-16.0.4/sparsehash-headers-2.0.2-bjbnaamdnumhjyvmmpk57yxhm2yvuezy" CACHE PATH "")

# boost headers from uberenv
set("BOOST_DIR" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_03_08_13_59_54/spack/opt/spack/toss_3_x86_64_ib/intel-16.0.4/boost-headers-1.58.0-22qg7in5tkbmwxpv2wszzejsmhqdpwmt" CACHE PATH "")

# python from uberenv
set("PYTHON_EXECUTABLE" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_03_08_13_59_54/spack/opt/spack/toss_3_x86_64_ib/intel-16.0.4/python-2.7.11-63apwbc3kjq5b3xpyjkcxpyf6yv5i354/bin/python" CACHE PATH "")

# lua from uberenv
set("LUA_DIR" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_03_08_13_59_54/spack/opt/spack/toss_3_x86_64_ib/intel-16.0.4/lua-5.1.5-glh6zpy27o4lrbzvre4k6qu3rkhjawfu" CACHE PATH "")

# doxygen from uberenv
set("DOXYGEN_EXECUTABLE" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_03_08_13_59_54/spack/opt/spack/toss_3_x86_64_ib/intel-16.0.4/doxygen-1.8.11-tb3urqnfhpui4fp2br7uzk66ahyvmgjd/bin/doxygen" CACHE PATH "")

# sphinx from uberenv
set("SPHINX_EXECUTABLE" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_03_08_13_59_54/spack/opt/spack/toss_3_x86_64_ib/intel-16.0.4/python-2.7.11-63apwbc3kjq5b3xpyjkcxpyf6yv5i354/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
set("UNCRUSTIFY_EXECUTABLE" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_03_08_13_59_54/spack/opt/spack/toss_3_x86_64_ib/intel-16.0.4/uncrustify-0.61-nzydqv27ioishmus56myotierqugxn54/bin/uncrustify" CACHE PATH "")

# lcov and genhtml from uberenv
set("LCOV_PATH" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_03_08_13_59_54/spack/opt/spack/toss_3_x86_64_ib/intel-16.0.4/lcov-1.11-fiwdecmhqvashhs323poq5mdr6uqnrio/usr/bin/lcov" CACHE PATH "")

set("GENHTML_PATH" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_03_08_13_59_54/spack/opt/spack/toss_3_x86_64_ib/intel-16.0.4/lcov-1.11-fiwdecmhqvashhs323poq5mdr6uqnrio/usr/bin/genhtml" CACHE PATH "")

# Temporarily disable CXX11 on intel builds until we resolve issue ATK-619
set("BLT_CXX_STD" "c++98" CACHE PATH "")

##################################
# end uberenv host-config
##################################

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################
# Options added manually to 
# lc toss3 intel@16.0.4  host configs
##############################################################################

##############################################################################
# MPI - manually added for now
##############################################################################
set("ENABLE_MPI" ON CACHE PATH "")
set("MPI_C_COMPILER" "/usr/tce/packages/mvapich2/mvapich2-2.2-intel-16.0.4/bin/mpicc" CACHE PATH "")
set("MPI_CXX_COMPILER" "/usr/tce/packages/mvapich2/mvapich2-2.2-intel-16.0.4/bin/mpiccxx" CACHE PATH "")
set("MPI_Fortran_COMPILER" "/usr/tce/packages/mvapich2/mvapich2-2.2-intel-16.0.4/bin/mpif90" CACHE PATH "")

set("MPI_C_LIBRARIES"   "/usr/tce/packages/mvapich2/mvapich2-2.2-intel-16.0.4/lib/libmpi.a" CACHE PATH "")
set("MPI_CXX_LIBRARIES" "/usr/tce/packages/mvapich2/mvapich2-2.2-intel-16.0.4//lib/libmpicxx.a" CACHE PATH "")
set("MPI_CXX_INCLUDE_PATH" "/usr/tce/packages/mvapich2/mvapich2-2.2-intel-16.0.4/include" CACHE PATH "")

set("MPIEXEC" "/usr/bin/srun" CACHE PATH "")
set("MPIEXEC_NUMPROC_FLAG" "-n" CACHE PATH "")


##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

