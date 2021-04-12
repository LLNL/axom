# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Examples for using Axom within other build systems: Make, CMake, BLT
#------------------------------------------------------------------------------
# Note: Installed to 'examples/axom' to avoid name conflicts with other libs.
# This will crop up if we use spack to install axom as a python extension.
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# install 'using-with-make' example (no configuration, so just copy files)
#------------------------------------------------------------------------------
install(
    DIRECTORY 
        using-with-make
    DESTINATION 
        examples/axom
)

#------------------------------------------------------------------------------
# configure and install 'using-with-cmake' example
#------------------------------------------------------------------------------
configure_file(
    using-with-cmake/host-config.cmake.in
    ${CMAKE_BINARY_DIR}/examples/using-with-cmake/host-config.cmake)

install(
    FILES
        using-with-cmake/CMakeLists.txt
        using-with-cmake/example.cpp
        ${CMAKE_BINARY_DIR}/examples/using-with-cmake/host-config.cmake
    DESTINATION
        examples/axom/using-with-cmake
)

#------------------------------------------------------------------------------
# configure and install 'using-with-blt' example
#------------------------------------------------------------------------------
configure_file(
    using-with-blt/host-config.cmake.in
    ${CMAKE_BINARY_DIR}/examples/using-with-blt/host-config.cmake)

install(
    FILES
        using-with-blt/CMakeLists.txt
        using-with-blt/example.cpp
        ${CMAKE_BINARY_DIR}/examples/using-with-blt/host-config.cmake
    DESTINATION
        examples/axom/using-with-blt
)

