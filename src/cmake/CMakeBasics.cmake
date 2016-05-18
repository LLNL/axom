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
# Setup build options and their default values
################################
include(SetupCmakeOptions)

################################
# Prevent in-source builds
################################
include(PreventInSourceBuilds)

################################
#  macros
################################
include(ATKMacros)

################################
# Setup 3rd Party Libs
################################
include(SetupThirdParty)

################################
# Setup toolkit docs targets
################################
include(SetupDocs)

################################
# Setup toolkit generate targets
################################
include(SetupGenerate)

################################
# Setup toolkit source checks
################################
include(SetupCodeChecks)


# XXX this should move to some better place and be based on compiler
#set (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG} -fcheck=bounds")

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
     ${PROJECT_BINARY_DIR}/test
     CACHE PATH
     "Directory where test executables will go in the build tree"
     )

## Set the path were all example test executables will go
 set(EXAMPLE_OUTPUT_DIRECTORY
     ${PROJECT_BINARY_DIR}/example
     CACHE PATH
     "Directory where example executables will go in the build tree"
     )

 ## Set the Fortran module directory
 set(CMAKE_Fortran_MODULE_DIRECTORY
     ${PROJECT_BINARY_DIR}/lib/fortran
     CACHE PATH
     "Directory where all Fortran modules will go in the build tree"
     )

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
include(SetupCompilerOptions)

################################
# Setup code metrics -
# profiling, code coverage, etc.
# (must be included after SetupCompilerOptions)
################################
include(SetupCodeMetrics)

################################
# Standard CMake Options
################################

include(ExternalProject)
if (ENABLE_TESTS)
  include(CTest)

  ## add google test
  add_subdirectory(${PROJECT_SOURCE_DIR}/thirdparty/gtest-1.7.0)
  blt_register_library(NAME gtest
                       INCLUDES ${gtest_SOURCE_DIR}/include
                       LIBRARIES gtest_main gtest
                       )

  ## Add Fruit   FortRan UnIT test
  if (ENABLE_FORTRAN)
    add_subdirectory(${PROJECT_SOURCE_DIR}/thirdparty/fruit-3.3.9)
  endif (ENABLE_FORTRAN)

  if(ENABLE_BENCHMARKS)
    ## add google benchmark
    add_subdirectory(${PROJECT_SOURCE_DIR}/thirdparty/gbenchmark)
    set(GBENCHMARK_INCLUDES ${benchmark_SOURCE_DIR}/include ${benchmark_SOURCE_DIR}
            CACHE INTERNAL "Google Benchmark include directories" FORCE)
    set(GBENCHMARK_LIBS benchmark
            CACHE INTERNAL "Google Benchmark link libraries" FORCE)

    #message(STATUS "Google benchmark -- \n\t inc -- ${GBENCHMARK_INCLUDES} -- \n\t lib -- ${GBENCHMARK_LIBS} ")
  
    # This sets up a target to run the benchmarks
    add_custom_target(run_benchmarks COMMAND ctest -C Benchmark -VV -R benchmark)

  endif()

  enable_testing()
  
#  add_dependencies(test run_benchmarks)

endif()

################################
# MPI
################################
message(STATUS "MPI Support is ${ENABLE_MPI}")
if (ENABLE_MPI)
  find_package(MPI REQUIRED)
  message(STATUS "MPI C Compile Flags: ${MPI_C_COMPILE_FLAGS}")
  message(STATUS "MPI C Include Path: ${MPI_C_INCLUDE_PATH}")
  message(STATUS "MPI C Link Flags: ${MPI_C_LINK_FLAGS}")
  message(STATUS "MPI C Libraries: ${MPI_C_LIBRARIES}")
endif()

################################
# OpenMP
################################
message(STATUS "OpenMP Support is ${ENABLE_OPENMP}")
if(ENABLE_OPENMP)
    find_package(OpenMP REQUIRED)
    message(STATUS "OpenMP CXX Flags: ${OpenMP_CXX_FLAGS}")
endif()
