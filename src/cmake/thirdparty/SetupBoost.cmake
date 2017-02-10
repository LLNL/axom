###############################################################################
#
# Setup BOOST
# This file defines:
#  BOOST_FOUND - If Boost was found
#  (Standard Boost Include dirs)

# set the vars that the standard cmake boost logic needs
set(ENABLE_BOOT ON)
set(BOOST_ROOT ${BOOST_DIR})

# find boost
find_package(Boost
             1.58
             REQUIRED)

MESSAGE(STATUS "Boost include dir: " ${Boost_INCLUDE_DIR})
MESSAGE(STATUS "Boost version: " ${Boost_VERSION})
# this is what we used to signal boost is on in our build system
set(BOOST_FOUND TRUE)