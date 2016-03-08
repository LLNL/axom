###############################################################################
#
# CMake Cache Seed file for chaos_5_x86_64_ib machines using gcc 4.9.3
#
###############################################################################

##################################
# uberenv host-config
##################################
# chaos_5_x86_64_ib-gcc@4.9.3
##################################

# cmake from uberenv
# cmake executable path: /usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/cmake-3.2.2-jnpt24vowmusazje5nogmtkk2mmriu7u/bin/cmake

#######
# using gcc@4.9.3 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/apps/gnu/4.9.3/bin/gcc" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/apps/gnu/4.9.3/bin/g++" CACHE PATH "")
set(GCOV_PATH "/usr/apps/gnu/4.9.3/bin/gcov" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE PATH "")
set(CMAKE_Fortran_COMPILER  "/usr/apps/gnu/4.9.3/bin/gfortran" CACHE PATH "")

# enable google benchmark functionality for gnu compilers
set(ENABLE_BENCHMARK ON CACHE PATH "")

# conduit from uberenv
set(CONDUIT_DIR "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/conduit-github-6y2swpnp6vthtm43ddawh4mn4jefan5j" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/python-2.7.8-wxs6x4molggy6niut2czpjfhvhmrscj5/bin/python" CACHE PATH "")

# sphinx from uberenv
set(SPHINX_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/python-2.7.8-wxs6x4molggy6niut2czpjfhvhmrscj5/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/uncrustify-0.61-um6yrn4vqk3iwjv4cj7f5vt3oqe2gihp/bin/uncrustify" CACHE PATH "")

# sparsehash headers from uberenv
set(SPARSEHASH_DIR "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/sparsehash-headers-2.0.2-tdigejugdkoolzqjke75rc7jd236fcm4" CACHE PATH "")

# boost headers from uberenv
set(ENABLE_BOOST ON CACHE PATH "")
set(BOOST_ROOT "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/boost-headers-1.58.0-jo2qfge5otlfcp7ahwgu27opxweoouga" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/lcov-1.11-mrnut4rtvsbdmkfvjf454flkk5wxilra/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/lcov-1.11-mrnut4rtvsbdmkfvjf454flkk5wxilra/usr/bin/genhtml" CACHE PATH "")

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
set(MPI_C_COMPILER "/usr/local/tools/mvapich2-gnu-2.0/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "/usr/local/tools/mvapich2-gnu-2.0/bin/mpicxx" CACHE PATH "")
set(MPI_Fortran_COMPILER  "/usr/local/tools/mvapich2-gnu-2.0/bin/mpif90" CACHE PATH "")
