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
option(ENABLE_BOOST "ENABLE BOOST" OFF)
if (ENABLE_BOOST)
    find_package(Boost 1.57  REQUIRED)
    MESSAGE(STATUS "Boost include dir: " ${Boost_INCLUDE_DIR})
    MESSAGE(STATUS "Boost version: " ${Boost_VERSION} )
endif()