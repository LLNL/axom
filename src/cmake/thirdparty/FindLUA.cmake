# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Setup Lua
#------------------------------------------------------------------------------
# This file defines:
#  LUA_FOUND - If Lua was found
#  LUA_INCLUDE_DIRS - The Lua include directories
#  LUA_LIBRARY - The Lua library
#------------------------------------------------------------------------------

# first Check for LUA_DIR
if (NOT EXISTS "${LUA_DIR}")
    message(FATAL_ERROR "Given LUA_DIR does not exist: ${LUA_DIR}")
endif()

if (NOT IS_DIRECTORY "${LUA_DIR}")
    message(FATAL_ERROR "Given LUA_DIR is not a directory: ${LUA_DIR}")
endif()

# Find includes directory
find_path( LUA_INCLUDE_DIR lua.hpp
           PATHS  ${LUA_DIR}/include/
                  ${LUA_DIR}/include/lua
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

# Find libraries
find_library( LUA_LIBRARY NAMES lua liblua
              PATHS ${LUA_DIR}/lib
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)


include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LUA_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LUA  DEFAULT_MSG
                                  LUA_INCLUDE_DIR
                                  LUA_LIBRARY )

if(NOT LUA_FOUND)
    message(FATAL_ERROR "LUA_DIR is not a path to a valid Lua install")
endif()

message(STATUS "Lua Includes: ${LUA_INCLUDE_DIR}")
message(STATUS "Lua Libraries: ${LUA_LIBRARY}")
