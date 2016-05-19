###############################################################################
# Copyright (c) 2014, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-666778
#
# All rights reserved.
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
# Prevent in-source builds
################################
include(blt/cmake/PreventInSourceBuilds.cmake)


################################
# Setup build options and their default values
################################
include(blt/cmake/BLTOptions.cmake)

################################
# enable ctest support 
################################
if(ENABLE_TESTS)
  enable_testing()
endif()


################################
#  macros
################################
include(blt/cmake/BLTMacros.cmake)

################################
# standard tpl support
################################
include(blt/cmake/thirdparty/SetupThirdParty.cmake)

################################
# Setup docs targets
################################
include(blt/cmake/SetupDocs.cmake)

################################
# Setup source checks
################################
include(blt/cmake/SetupCodeChecks.cmake)

################################
# Standard Build Layout
################################

##
## Defines the layout of the build directory. Namely,
## it indicates the location where the various header files should go,
## where to store libraries (static or shared), the location of the
## bin directory for all executables and the location for fortran modules.
##

## Set the path where all the header will be stored
 set(HEADER_INCLUDES_DIRECTORY
     ${PROJECT_BINARY_DIR}/include/
     CACHE PATH
     "Directory where all headers will go in the build tree"
     )
 include_directories(${HEADER_INCLUDES_DIRECTORY})

 ## Set the path where all the libraries will be stored
 set(LIBRARY_OUTPUT_PATH
     ${PROJECT_BINARY_DIR}/lib
     CACHE PATH
     "Directory where compiled libraries will go in the build tree"
     )

 ## Set the path where all the install executables will go
 set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
     ${PROJECT_BINARY_DIR}/bin
     CACHE PATH
     "Directory where executables will go in the build tree"
     )

## Set the path were all test executables will go
 set(TEST_OUTPUT_DIRECTORY
     ${PROJECT_BINARY_DIR}/tests
     CACHE PATH
     "Directory where test executables will go in the build tree"
     )

## Set the path were all example test executables will go
 set(EXAMPLE_OUTPUT_DIRECTORY
     ${PROJECT_BINARY_DIR}/examples
     CACHE PATH
     "Directory where example executables will go in the build tree"
     )

 ## Set the Fortran module directory
 set(CMAKE_Fortran_MODULE_DIRECTORY
     ${PROJECT_BINARY_DIR}/lib/fortran
     CACHE PATH
     "Directory where all Fortran modules will go in the build tree"
     )
#
# TODO: (Cyrus's Note I don't think this should be here?)
#
## Set the Lua module directory
 set(BLT_Lua_MODULE_DIRECTORY
     "${PROJECT_BINARY_DIR}/lib/lua"
     CACHE PATH
     "Directory where all Lua modules will go in the build tree"
 )

## Mark as advanced
mark_as_advanced(
     LIBRARY_OUTPUT_PATH
     CMAKE_RUNTIME_OUTPUT_DIRECTORY
     CMAKE_Fortran_MODULE_DIRECTORY
     )

################################
# Setup compiler options
# (must be included after HEADER_INCLUDES_DIRECTORY is set)
################################
include(blt/cmake/SetupCompilerOptions.cmake)

################################
# Setup code metrics -
# profiling, code coverage, etc.
# (must be included after SetupCompilerOptions)
################################
include(blt/cmake/SetupCodeMetrics.cmake)

################################
# builtin third party libs used by blt
################################
add_subdirectory(blt/thirdparty_builtin)

################################
# blt smoke tests
################################
add_subdirectory(blt/tests)




