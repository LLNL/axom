#------------------------------------------------------------------------------
# Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Setup BOOST
#
# This file defines:
#  BOOST_FOUND - If Boost was found
#  (Standard Boost Include dirs)
#------------------------------------------------------------------------------

# set the vars that the standard cmake boost logic needs
set(BOOST_ROOT ${BOOST_DIR})

# find boost
find_package(Boost
             1.54
             REQUIRED)

MESSAGE(STATUS "Boost include dir: " ${Boost_INCLUDE_DIR})
MESSAGE(STATUS "Boost version: " ${Boost_VERSION})

# this is what we used to signal boost is on in our build system
set(BOOST_FOUND ${Boost_FOUND})
