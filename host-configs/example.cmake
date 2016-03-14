###############################################################################
#
# Example CMake Cache Seed file for the asctoolkit
#
###############################################################################
# 
# When passed to cmake with "-C", this file initializes CMake variables 
# used during the configure process. 
#
# Variables are initialized using  CMake `set` commands, such as:
# 
#  set(VAR_NAME VALUE CACHE PATH "") 
#
# The `CACHE PATH ""` argument is important to ensure the values are cached, 
# and actually used. 
#
###############################################################################


###############################################################################
# Select the c and c++ compiler though the standard CMake Variables.
###############################################################################
set(CMAKE_C_COMPILER "/usr/apps/gnu/4.7.1/bin/gcc" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/apps/gnu/4.7.1/bin/g++" CACHE PATH "")

###############################################################################
# Additional Compiler Flags
###############################################################################

# additional flags for all configurations 
set(EXTRA_CXX_FLAGS "-DEXTRA_FLAGS_CXX_DEFINE" CACHE PATH "")

# additional flags for debug builds
set(EXTRA_CXX_DEBUG_FLAGS "-DEXTRA_CXX_DEBUG_FLAGS_DEFINE" CACHE PATH "")

# additional flags for debug builds
set(EXTRA_CXX_RELEASE_FLAGS "-DEXTRA_CXX_RELEASE_FLAGS_DEFINE" CACHE PATH "")

