###############################################################################
# Copyright (c) 2014, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-666778
#
# All rights reserved.
#
# This file is part of Conduit.
#
# For details, see https://lc.llnl.gov/conduit/.
#
# Please also read conduit/LICENSE
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the disclaimer below.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the disclaimer (as noted below) in the
#   documentation and/or other materials provided with the distribution.
#
# * Neither the name of the LLNS/LLNL nor the names of its contributors may
#   be used to endorse or promote products derived from this software without
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
# LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
# IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
###############################################################################

################################
# Standard Build Layout
################################

##
## Defines the layout of the build directory. Namely,
## it indicates the location where the various header files should go,
## where to store libraries (static or shared), the location of the
## bin directory for all executables and the location for fortran moudules.
##

## Set the path where all the header will be stored
 set(HEADER_INCLUDES_DIRECTORY
     ${PROJECT_BINARY_DIR}/include
     CACHE PATH
     "Directory where all headers will go in the build tree"
     )
 include_directories(${HEADER_INCLUDES_DIRECTORY})

 ## Set the path where all the libraries will be stored
 set(LIBRARY_OUTPUT_PATH
     ${PROJECT_BINARY_DIR}/libs
     CACHE PATH
     "Directory where compiled libraries will go in the build tree"
     )

 ## Set the path where all the executables will go
 set(EXECUTABLE_OUTPUT_PATH
     ${PROJECT_BINARY_DIR}/bin
     CACHE PATH
     "Directory where executables will go in the build tree"
     )

 ## Set the Fortran module directory
 set(CMAKE_Fortran_MODULE_DIRECTORY
     ${PROJECT_BINARY_DIR}/FortranModules
     CACHE PATH
     "Directory where all Fortran modules will go in the build tree"
     )

## Mark as advanced
mark_as_advanced(
     LIBRARY_OUTPUT_PATH
     EXECUTABLE_OUTPUT_PATH
     CMAKE_Fortran_MODULE_DIRECTORY
     )


################################
# Standard CMake Options
################################

# enable creation of compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

## Choose static or shared libraries.
option(BUILD_SHARED_LIBS "Build shared libraries." OFF)

option(ENABLE_WARNINGS "Enable Compiler warnings." OFF)
if(ENABLE_WARNINGS)

    # set the warning levels we want to abide by
    if(CMAKE_BUILD_TOOL MATCHES "(msdev|devenv|nmake)")
        add_definitions(/W2)
    else()
        if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR
            "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
            # using clang or gcc
            add_definitions(-Wall -Wextra)
        endif()
    endif()

endif()

################################
# RPath Settings
################################

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
   set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

################################
#  macros
################################

##------------------------------------------------------------------------------
## - Given a library target nane and a list of sources, this macros builds
##   a static or shared library according to a user-supplied BUILD_SHARED_LIBS
##   option.
##
## make_library(libtarget srcs)
##------------------------------------------------------------------------------
macro(make_library libtarget srcs)

   if( BUILD_SHARED_LIBS )
      add_library(${libtarget} SHARED ${srcs})
   else()
      add_library(${libtarget} STATIC ${srcs})
   endif()

endmacro(make_library)

