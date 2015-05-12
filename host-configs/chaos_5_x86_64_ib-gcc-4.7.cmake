###############################################################################
#
# CMake Cache Seed file for chaos_5_x86_64_ib machines using gcc 4.7
#
###############################################################################

###############################################################################
# Select the c and c++ compiler though the standard CMake Variables.
###############################################################################
set(CMAKE_C_COMPILER "/usr/apps/gnu/4.7.1/bin/gcc" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/apps/gnu/4.7.1/bin/g++" CACHE PATH "")
set(SPHINX_EXECUTABLE "/usr/gapps/asctoolkit/tools/spack/opt/chaos_5_x86_64_ib/gcc\@4.4.7/python\@2.7.8-703c7a96/bin/sphinx-build" CACHE PATH "")


###############################################################################
# Set location for boost library
###############################################################################
set(BOOST_ROOT "/usr/local/tools/boost" CACHE PATH "")

###############################################################################
# Additional Compiler Flags
###############################################################################

# additional flags for all configurations 
#set(EXTRA_CXX_FLAGS "-DEXTRA_FLAGS_CXX_DEFINE" CACHE PATH "")

# additional flags for debug builds
# set(EXTRA_CXX_DEBUG_FLAGS "-DEXTRA_CXX_DEBUG_FLAGS_DEFINE" CACHE PATH "")

# additional flags for debug builds
# set(EXTRA_CXX_RELEASE_FLAGS "-DEXTRA_CXX_RELEASE_FLAGS_DEFINE" CACHE PATH "")
