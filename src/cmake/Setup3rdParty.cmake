####################################
# Datastore 3rd Party Dependencies
####################################

################################
# Conduit
################################
include(cmake/FindConduit.cmake)


################################
# Sparsehash
################################
include(cmake/FindSparsehash.cmake)

################################
# Documentation Packages
################################

find_package(Doxygen)
include(cmake/FindSphinx.cmake)

################################
# linting via Uncrustify
################################
include(cmake/FindUncrustify.cmake)

################################
# Find boost headers
#
# This should only be enabled if we aren't building with C++11 support.  The
# policy is to use boost only to replace C++11 assets for older (non-c++11)
# compilers.
# For developers wanting to evaluate or explore code additions and want to
# use boost, they override the below logic by modifying this file, but should
# do so in a private branch only.
################################
if (NOT ENABLE_CXX11)
  find_package(Boost
                1.55
                REQUIRED)
  MESSAGE(STATUS "Boost include dir: " ${Boost_INCLUDE_DIR})
  MESSAGE(STATUS "Boost version: " ${Boost_VERSION} )
endif()


