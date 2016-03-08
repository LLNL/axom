###############################################################################
#
# CMake Cache Seed file for bgqos_0 machines using xlc.
#
###############################################################################

###############################################################################
# Select the c and c++ compiler though the standard CMake Variables.
###############################################################################
set(CMAKE_C_COMPILER "bgxlc_r" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/local/tools/compilers/ibm/bgxlc++_r-12.1.0.12a" CACHE PATH "")
#change back to floor version after LC has pushed out 12.1.0.12a to floor version.
#set(CMAKE_CXX_COMPILER "bgxlc++_r" CACHE PATH "")
set(CMAKE_Fortran_COMPILER "bgxlf2003_r" CACHE PATH "")

#######
# uberenv host-config for asctoolkit
#######
# cmake from uberenv
# cmake executable path: /usr/gapps/asctoolkit/thirdparty_libs/spack/opt/bgqos_0/gcc@4.4.7/cmake@3.2.2/bin/cmake

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

# additional flag to disable warnings about omp pragmas when omp is not enabled
set(DISABLE_OMP_PRAGMA_WARNINGS " -qignprag=omp " CACHE STRING "Flag to disable warning about pragmas when omp is disabled")

# additional flags for all configurations 
#set(EXTRA_CXX_FLAGS "-DEXTRA_FLAGS_CXX_DEFINE" CACHE PATH "")

# additional flags for debug builds
# set(EXTRA_CXX_DEBUG_FLAGS "-DEXTRA_CXX_DEBUG_FLAGS_DEFINE" CACHE PATH "")

# additional flags for debug builds
# set(EXTRA_CXX_RELEASE_FLAGS "-DEXTRA_CXX_RELEASE_FLAGS_DEFINE" CACHE PATH "")
