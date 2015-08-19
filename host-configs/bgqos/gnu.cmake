###############################################################################
#
# CMake Cache Seed file for bgqos_0 machines using gcc/g++ 4.7.2.
#
###############################################################################

###############################################################################
# Select the c and c++ compiler though the standard CMake Variables.
###############################################################################
set(CMAKE_C_COMPILER "/usr/local/tools/toolchain-4.7.2/scripts/bggcc" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/local/tools/toolchain-4.7.2/scripts/bgg++" CACHE PATH "")
set(CMAKE_Fortran_COMPILER "/usr/local/tools/toolchain-4.7.2/scripts/bggfortran" CACHE PATH "")

#######
# uberenv host-config for asctoolkit
#######
# cmake from uberenv
# cmake exectuable path: /usr/gapps/asctoolkit/thirdparty_libs/spack/opt/bgqos_0/gcc@4.4.7/cmake@3.2.2/bin/cmake

# python from uberenv
#set(PYTHON_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/bgqos_0/gcc@4.4.7/python@2.7.8/bin/python" CACHE PATH "")

# sphinx from uberenv
#set(SPHINX_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/bgqos_0/gcc@4.4.7/python@2.7.8/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
#set(UNCRUSTIFY_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/bgqos_0/gcc@4.4.7/uncrustify@0.61/bin/uncrustify" CACHE PATH "")

# boost headers from uberenv
set(ENABLE_BOOST ON CACHE PATH "")
set(BOOST_ROOT "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/chaos_5_x86_64_ib/gcc@4.4.7/boost-headers@1.58.0" CACHE PATH "")

###############################################################################
# Additional Compiler Flags
###############################################################################

# additional flags for all configurations 
#set(EXTRA_CXX_FLAGS "-DEXTRA_FLAGS_CXX_DEFINE" CACHE PATH "")

# additional flags for debug builds
# set(EXTRA_CXX_DEBUG_FLAGS "-DEXTRA_CXX_DEBUG_FLAGS_DEFINE" CACHE PATH "")

# additional flags for debug builds
# set(EXTRA_CXX_RELEASE_FLAGS "-DEXTRA_CXX_RELEASE_FLAGS_DEFINE" CACHE PATH "")
