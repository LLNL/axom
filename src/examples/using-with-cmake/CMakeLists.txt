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
###############################################################################
#
# Example that shows how to use an installed instance of Axom in another
# CMake-based build system.
#
# To build:
#  mkdir build
#  cd build
#  cmake -DAXOM_DIR={axom install path} ../
#  make
#  ./example
#
#
# Note:
#  If you use features that leverage external third party libraries 
#  (ex: Conduit features), the imported CMake targets will handle linking 
#  these libraries, however you may need to add additional include
#  include directories to use their headers. 
# 
###############################################################################

cmake_minimum_required(VERSION 3.8)

project(using_with_cmake)

include("FindAxom.cmake")

# setup the axom include path
include_directories(${AXOM_INCLUDE_DIRS})

# create our example 
add_executable(example example.cpp)

# link to axom targets
target_link_libraries(example axom_utils)

