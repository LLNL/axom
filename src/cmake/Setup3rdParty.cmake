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
################################
if (ENABLE_BOOST)
  if (DEFINED BOOST_ROOT)
    find_package(Boost
                 1.55
                 REQUIRED)
    MESSAGE(STATUS "Boost include dir: " ${Boost_INCLUDE_DIR})
    MESSAGE(STATUS "Boost version: " ${Boost_VERSION} )
  else()
    MESSAGE(FATAL_ERROR "ENABLE_BOOST is true, but BOOST_ROOT was not set.  Check your host-config file.")
  endif()
endif()
