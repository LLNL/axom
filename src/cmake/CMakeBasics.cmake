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
include(cmake/ATKOptions.cmake)

################################
# Setup toolkit generate targets
################################
include(cmake/SetupShroud.cmake)


################################
<<<<<<< HEAD
# Setup toolkit generate targets
=======
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

################################
# Setup compiler options
# (must be included after HEADER_INCLUDES_DIRECTORY and MPI variables are set)
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
include(cmake/thirdparty/SetupATKThirdParty.cmake)


#
# We don't try to use this approach for CMake generators that support
# multiple configurations. See: CZ JIRA: ATK-45
#
if(NOT CMAKE_CONFIGURATION_TYPES)
    ######################################################
    # Add define we can use when debug builds are enabled
    ######################################################
    if( (CMAKE_BUILD_TYPE MATCHES Debug)
        OR (CMAKE_BUILD_TYPE MATCHES RelWithDebInfo )
      )
        add_definitions(-DATK_DEBUG)
    endif()

endif()

################################
# Fortran Configuration
################################
if(ENABLE_FORTRAN)
    # default property to free form
    set(CMAKE_Fortran_FORMAT FREE)

    # Create macros for Fortran name mangling
    include(FortranCInterface)
    FortranCInterface_HEADER(${HEADER_INCLUDES_DIRECTORY}/common/FC.h MACRO_NAMESPACE "FC_")

    if (ENABLE_MPI)
        # Determine if we should use fortran mpif.h header or fortran mpi module
        find_path(mpif_path
            NAMES "mpif.h"
            PATHS ${MPI_Fortran_INCLUDE_PATH}
            NO_DEFAULT_PATH
            )
        
        if(mpif_path)
            set(MPI_Fortran_USE_MPIF ON CACHE PATH "")
            message(STATUS "Using MPI Fortran header: mpif.h")
        else()
            set(MPI_Fortran_USE_MPIF OFF CACHE PATH "")
            message(STATUS "Using MPI Fortran module: mpi.mod")
        endif()
    endif()
endif()


##############################################################################
# Setup some additional compiler options that can be useful in various targets
# These are stored in their own variables.
# Usage: To add one of these sets of flags to some source files:
#   get_source_file_property(_origflags <src_file> COMPILE_FLAGS)
#   set_source_files_properties(<list_of_src_files> 
#        PROPERTIES COMPILE_FLAGS "${_origFlags} ${<flags_variable}" )
##############################################################################

set(custom_compiler_flags_list) # Tracks custom compiler flags for logging

# Flag for disabling warnings about omp pragmas in the code
blt_append_custom_compiler_flag(FLAGS_VAR ATK_DISABLE_OMP_PRAGMA_WARNINGS
                  DEFAULT "-Wno-unknown-pragmas"
                  XL      "-qignprag=omp"
                  INTEL   "-diag-disable 3180"
                  )
list(APPEND custom_compiler_flags_list ATK_DISABLE_OMP_PRAGMA_WARNINGS)

# Flag for disabling warnings about unused parameters.
# Useful when we include external code.
blt_append_custom_compiler_flag(FLAGS_VAR ATK_DISABLE_UNUSED_PARAMETER_WARNINGS
                  DEFAULT "-Wno-unused-parameter"
                  XL      "-qinfo=nopar"
                  )
list(APPEND custom_compiler_flags_list ATK_DISABLE_UNUSED_PARAMETER_WARNINGS)

# Flag for disabling warnings about unused variables
# Useful when we include external code.
blt_append_custom_compiler_flag(FLAGS_VAR ATK_DISABLE_UNUSED_VARIABLE_WARNINGS
                  DEFAULT "-Wno-unused-variable"
                  XL      "-qinfo=nouse"
                  )
list(APPEND custom_compiler_flags_list ATK_DISABLE_UNUSED_VARIABLE_WARNINGS)

# Flag for disabling warnings about variables that may be uninitialized.
# Useful when we are using compiler generated interface code (e.g. in shroud)
blt_append_custom_compiler_flag(FLAGS_VAR ATK_DISABLE_UNINITIALIZED_WARNINGS
                  DEFAULT "-Wno-uninitialized"
                  XL      "-qsuppress=1540-1102"
                  )
list(APPEND custom_compiler_flags_list ATK_DISABLE_UNINITIALIZED_WARNINGS)

# Flag for disabling warnings about strict aliasing.
# Useful when we are using compiler generated interface code (e.g. in shroud)
blt_append_custom_compiler_flag(FLAGS_VAR ATK_DISABLE_ALIASING_WARNINGS
                  DEFAULT "-Wno-strict-aliasing"
                  XL      ""
                  )
list(APPEND custom_compiler_flags_list ATK_DISABLE_ALIASING_WARNINGS)

# Flag for enabling the C preprocessor in fortran.
# (Note KW 5/2016) The XL flag only applies to *.f90 files -- I could not find a more general solution.   
#     xlf only allows one file remapping at a time. If you have *.f files, '-qsuffix=cpp=f' should work.
#     Alternatively, you can rename the file's extension to automatically invoke the preprocessor (e.g. *.f90 ->  *.F90)
blt_append_custom_compiler_flag(FLAGS_VAR ATK_PREPROCESS_FORTRAN
                  DEFAULT "-cpp"
                  XL      "-qsuffix=cpp=f90"  # Note: only invokes the preprocessor on files with extension *.f90
                  )
list(APPEND custom_compiler_flags_list ATK_PREPROCESS_FORTRAN)

   
# message(STATUS "Custom compiler flags:")
# foreach(flag ${custom_compiler_flags_list})
#    message(STATUS "\tvalue of ${flag} is '${${flag}}'")
# endforeach()
