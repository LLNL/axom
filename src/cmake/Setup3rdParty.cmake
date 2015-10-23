####################################
# Datastore 3rd Party Dependencies
####################################

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
  find_package(Boost
                1.55
                REQUIRED)
  MESSAGE(STATUS "Boost include dir: " ${Boost_INCLUDE_DIR})
  MESSAGE(STATUS "Boost version: " ${Boost_VERSION} )
endif()

find_package(Sparsehash
              2.0.2
              REQUIRED)
MESSAGE(STATUS "Sparsehash include dir: " ${Sparsehash_INCLUDE_DIR})
