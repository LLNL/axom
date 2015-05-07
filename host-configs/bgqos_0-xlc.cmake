###############################################################################
#
# CMake Cache Seed file for bgqos_0 machines using xlc.
#
###############################################################################

###############################################################################
# Select the c and c++ compiler though the standard CMake Variables.
###############################################################################
set(CMAKE_C_COMPILER "bgxlc_r" CACHE PATH "")
set(CMAKE_CXX_COMPILER "bgxlc++_r" CACHE PATH "")


###############################################################################
# Set location for boost library
# Need to provide your own boost and set this by hand, no boost 1.57 on BG/Q.
###############################################################################
#set(BOOST_ROOT "/usr/local/tools/boost" CACHE PATH "")

###############################################################################
# Additional Compiler Flags
###############################################################################

# additional flags for all configurations 
#set(EXTRA_CXX_FLAGS "-DEXTRA_FLAGS_CXX_DEFINE" CACHE PATH "")

# additional flags for debug builds
# set(EXTRA_CXX_DEBUG_FLAGS "-DEXTRA_CXX_DEBUG_FLAGS_DEFINE" CACHE PATH "")

# additional flags for debug builds
# set(EXTRA_CXX_RELEASE_FLAGS "-DEXTRA_CXX_RELEASE_FLAGS_DEFINE" CACHE PATH "")
