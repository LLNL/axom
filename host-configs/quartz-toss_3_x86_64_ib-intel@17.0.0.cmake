##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# toss_3_x86_64_ib-intel@17.0.0
##################################

# cmake from uberenv
# cmake executable path: /usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_04_22_50_48/spack/opt/spack/toss_3_x86_64_ib/intel-17.0.0/cmake-3.3.1-2uuw7ygf74nmfkdbjbwxwwt3oa7yyzht/bin/cmake

#######
# using intel@17.0.0 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/tce/packages/intel/intel-16.0.4/bin/icc" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/tce/packages/intel/intel-16.0.4/bin/icpc" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE PATH "")

set(CMAKE_Fortran_COMPILER  "/usr/tce/packages/intel/intel-16.0.4/bin/ifort" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_04_22_50_48/spack/opt/spack/toss_3_x86_64_ib/intel-17.0.0/hdf5-1.8.16-judt653wkzafz7ke3pbf2bvqp5tspfmr" CACHE PATH "")

# conduit from uberenv
set(CONDUIT_DIR "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_04_22_50_48/spack/opt/spack/toss_3_x86_64_ib/intel-17.0.0/conduit-0.2.0-cdjr3kcdfldi5xeuy5l3cpezk2k5733c" CACHE PATH "")

# doxygen from uberenv
set(DOXYGEN_EXECUTABLE "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_04_22_50_48/spack/opt/spack/toss_3_x86_64_ib/intel-17.0.0/doxygen-1.8.11-2b5e7xoxrsl5z7mzbefco5outba4zlpb/bin/doxygen" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_04_22_50_48/spack/opt/spack/toss_3_x86_64_ib/intel-17.0.0/python-2.7.11-qylimm5rmrrg33kovnlorkgwqf4ctzcf/bin/python" CACHE PATH "")

# lua from uberenv
set(LUA_DIR "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_04_22_50_48/spack/opt/spack/toss_3_x86_64_ib/intel-17.0.0/lua-5.1.5-bo4m2vndivukyjnocq5kgeirq5amzo6z" CACHE PATH "")

# sphinx from uberenv
set(SPHINX_EXECUTABLE "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_04_22_50_48/spack/opt/spack/toss_3_x86_64_ib/intel-17.0.0/python-2.7.11-qylimm5rmrrg33kovnlorkgwqf4ctzcf/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_04_22_50_48/spack/opt/spack/toss_3_x86_64_ib/intel-17.0.0/uncrustify-0.61-soowqg4ifccejpz3dvruzjyajujsgjhr/bin/uncrustify" CACHE PATH "")

# sparsehash headers from uberenv
set(SPARSEHASH_DIR "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_04_22_50_48/spack/opt/spack/toss_3_x86_64_ib/intel-17.0.0/sparsehash-headers-2.0.2-kdirqgdrnp55wwxc6ytoytz2vpwjufjh" CACHE PATH "")

# boost headers from uberenv
set(ENABLE_BOOST ON CACHE PATH "")
set(BOOST_ROOT "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_04_22_50_48/spack/opt/spack/toss_3_x86_64_ib/intel-17.0.0/boost-headers-1.58.0-6u7zrofaj4sdmjiodxelsplx7o32ewqz" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_04_22_50_48/spack/opt/spack/toss_3_x86_64_ib/intel-17.0.0/lcov-1.11-5at4wjezbyejj53ou4logab457thepy3/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_04_22_50_48/spack/opt/spack/toss_3_x86_64_ib/intel-17.0.0/lcov-1.11-5at4wjezbyejj53ou4logab457thepy3/usr/bin/genhtml" CACHE PATH "")

# Temporarily disable CXX11 on intel builds until we resolve issue ATK-619
set(BLT_CXX_STD "c++98" CACHE PATH "")
##################################
# end uberenv host-config
##################################

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################
# Options added manually to 
# lc toss3 intel@17.0.0  host configs
##############################################################################

##############################################################################
# MPI - manually added for now
##############################################################################
set("ENABLE_MPI" ON CACHE PATH "")
set("MPI_C_COMPILER" "/usr/tce/packages/mvapich2/mvapich2-2.2-intel-17.0.0/bin/mpicc" CACHE PATH "")
set("MPI_CXX_COMPILER" "/usr/tce/packages/mvapich2/mvapich2-2.2-intel-17.0.0/bin/mpiccxx" CACHE PATH "")
set("MPI_Fortran_COMPILER" "/usr/tce/packages/mvapich2/mvapich2-2.2-intel-17.0.0/bin/mpif90" CACHE PATH "")

set("MPIEXEC" "/usr/bin/srun" CACHE PATH "")
set("MPIEXEC_NUMPROC_FLAG" "-n" CACHE PATH "")


##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

