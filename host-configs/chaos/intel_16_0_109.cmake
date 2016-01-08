###############################################################################
#
# CMake Cache Seed file for chaos_5_x86_64_ib machines using icc/icpc
#
###############################################################################

##################################
# uberenv host-config
##################################
# chaos_5_x86_64_ib-intel@16.0.0
##################################

# cmake from uberenv
# cmake exectuable path: /usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.0/cmake-3.2.2-xl2rbk2uazutimxpt4dttzyba52aopnl/bin/cmake

#######
# using intel@16.0.0 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/local/bin/icc-16.0.109" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/local/bin/icpc-16.0.109" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE PATH "")

set(CMAKE_Fortran_COMPILER  "/usr/local/bin/ifort-16.0.109" CACHE PATH "")

# conduit from uberenv
set(CONDUIT_DIR "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.0/conduit-github-pnols5i47e5ho5i3xy5nvrxvb4aduszy" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.0/python-2.7.8-osaphsdmxddmc7iesqrvgvbcl5gkiza3/bin/python" CACHE PATH "")

# sphinx from uberenv
set(SPHINX_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.0/python-2.7.8-osaphsdmxddmc7iesqrvgvbcl5gkiza3/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.0/uncrustify-0.61-bdb34qerzc4hyzwcsskx4sqlgszigu6i/bin/uncrustify" CACHE PATH "")

# sparsehash headers from uberenv
set(SPARSEHASH_DIR "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.0/sparsehash-headers-2.0.2-icz2df4kuxs25a43nuklriz5tmbzbqn3" CACHE PATH "")

# boost headers from uberenv
set(ENABLE_BOOST ON CACHE PATH "")
set(BOOST_ROOT "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.0/boost-headers-1.58.0-pmcsxznf2dvn6qrtbix5jiycpcc4vuqj" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.0/lcov-1.11-3ykfvl67oshgejtlgkgus66wn7umx32k/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.0/lcov-1.11-3ykfvl67oshgejtlgkgus66wn7umx32k/usr/bin/genhtml" CACHE PATH "")

##################################
# end uberenv host-config
##################################


###############################################################################
# Additional Compiler Flags
###############################################################################

# additional flags for all configurations 
#set(EXTRA_CXX_FLAGS "-DEXTRA_FLAGS_CXX_DEFINE" CACHE PATH "")

# additional flags for debug builds
# set(EXTRA_CXX_DEBUG_FLAGS "-DEXTRA_CXX_DEBUG_FLAGS_DEFINE" CACHE PATH "")

# additional flags for debug builds
# set(EXTRA_CXX_RELEASE_FLAGS "-DEXTRA_CXX_RELEASE_FLAGS_DEFINE" CACHE PATH "")

#######
# MPI 
#######
set(ENABLE_MPI ON CACHE PATH "")
set(MPI_C_COMPILER "/usr/local/tools/mvapich2-intel-2.0/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "/usr/local/tools/mvapich2-intel-2.0/bin/mpicxx" CACHE PATH "")
set(MPI_Fortran_COMPILER  "/usr/local/tools/mvapich2-intel-2.0/bin/mpif90" CACHE PATH "")
