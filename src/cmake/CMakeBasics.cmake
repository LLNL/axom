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

 ## Set the path where all the executables will go
 set(EXECUTABLE_OUTPUT_PATH
     ${PROJECT_BINARY_DIR}/bin
     CACHE PATH
     "Directory where executables will go in the build tree"
     )

 ## Set the Fortran module directory
 set(CMAKE_Fortran_MODULE_DIRECTORY
     ${PROJECT_BINARY_DIR}/lib/fortran
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
# Check if we want to build Fortran support.
################################
option(ENABLE_FORTRAN "Enables Fortran compiler support." ON)

if(ENABLE_FORTRAN)
    add_definitions(-DATK_ENABLE_FORTRAN)

    # if enabled but no fortran compiler, halt the configure
    if(CMAKE_Fortran_COMPILER)
        MESSAGE(STATUS  "Fortran support enabled. (ENABLE_FORTRAN == ON, Fortran compiler found.)")
    else()
        MESSAGE(FATAL_ERROR "Fortran support selected, but no Fortran compiler was found.")
    endif()    
else()
    MESSAGE(STATUS  "Fortran support disabled.  (ENABLE_FORTRAN == OFF)")
endif()
 
################################
# Enable code coverage via gcov
# Note: Only supported for gnu.
################################
option(ENABLE_CODECOV "Enable/disable code coverage via gcov." OFF)

# These should be set in some separate macro later that sets coverage flags for each compiler.
#SET(GCC_COVERAGE_COMPILE_FLAGS "-fprofile-arcs -ftest-coverage")
#SET(GCC_COVERAGE_LINK_FLAGS    "-lgcov")
SET(GCC_COVERAGE_COMPILE_FLAGS "--coverage")
SET(GCC_COVERAGE_LINK_FLAGS    "--coverage")

if (CMAKE_BUILD_TYPE MATCHES Debug)
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        MESSAGE(STATUS "Debug build is using gnu, code coverage via gcov enabled.")
        SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
        SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )
    else()
        MESSAGE(STATUS "Debug build is not using gnu, code coverage is disabled.")
        SET(ENABLE_CODECOV OFF)
    endif()
endif()


################################
# Standard CMake Options
################################

include(ExternalProject)

option(BUILD_TESTING "Builds unit tests" ON)
if (BUILD_TESTING)

  ## add catch
  include_directories(${PROJECT_SOURCE_DIR}/TPL/catch/single_include)

  ## add google test
  add_subdirectory(${PROJECT_SOURCE_DIR}/TPL/gtest-1.7.0)
  set(GTEST_INCLUDES ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR}
            CACHE INTERNAL "GoogleTest include directories" FORCE)
  set(GTEST_LIBS gtest_main gtest
            CACHE INTERNAL "GoogleTest link libraries" FORCE)

  ## Add Fruit   FortRan UuIT test
  add_subdirectory(${PROJECT_SOURCE_DIR}/TPL/fruit-3.3.9)

  enable_testing()

endif()

################################
# MPI
################################
option(ENABLE_MPI "ENABLE MPI" OFF)
if (ENABLE_MPI)
  find_package(MPI REQUIRED)
endif()

################################
# OpenMP
################################
option(ENABLE_OMP "ENABLE OpenMP" OFF)
if(ENABLE_OMP)
    find_package(OpenMP REQUIRED)
    if (OPENMP_FOUND)
        message(STATUS "Found OpenMP")
    else()
        message(STATUS "Could not find OpenMP")
    endif()    
endif()

## Enable ENABLE C++ 11 features
option(ENABLE_CXX11 "Enables C++11 features" OFF)
if (ENABLE_CXX11)
  # define a macro so the code can ifdef accordingly.
  add_definitions("-DUSE_CXX11")
endif()

## Choose static or shared libraries.
option(BUILD_SHARED_LIBS "Build shared libraries." OFF)

option(ENABLE_WARNINGS "Enable Compiler warnings." ON)
if(ENABLE_WARNINGS)

    # set the warning levels we want to abide by
    if(CMAKE_BUILD_TOOL MATCHES "(msdev|devenv|nmake)")
        add_definitions(/W2)
    else()
        if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
            # using gcc
            add_definitions(-Wall -Wextra -Werror)
        endif()
    endif()

endif()


#############################################
# Support extra compiler flags and defines
#############################################
#
# We don't try to use this approach for CMake generators that support
# multiple configurations. See: CZ JIRA: ATK-45
#
if(NOT CMAKE_CONFIGURATION_TYPES)

    ######################################################
    # Add define we can use when debug builds are enabled
    ######################################################
    if(CMAKE_BUILD_TYPE MATCHES Debug)
        add_definitions(-DATK_DEBUG)
    endif()

    ##########################################
    # Support Extra Flags for the C compiler.
    ##########################################
    if(EXTRA_C_FLAGS)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${EXTRA_C_FLAGS}")
    endif()

    # Extra Flags for the debug builds with the C compiler.
    if(EXTRA_C_FLAGS_DEBUG AND CMAKE_BUILD_TYPE MATCHES Debug)
        add_compile_options("${EXTRA_C_FLAGS_DEBUG}")
    endif()

    # Extra Flags for the release builds with the C compiler.
    if(EXTRA_C_FLAGS_RELEASE AND CMAKE_BUILD_TYPE MATCHES RELEASE)
        add_compile_options("${EXTRA_C_FLAGS_RELEASE}")
    endif()

    #############################################
    # Support Extra Flags for the C++ compiler.
    #############################################
    if(EXTRA_CXX_FLAGS)
        add_compile_options("${EXTRA_CXX_FLAGS}")
    endif()

    # Extra Flags for the debug builds with the C++ compiler.
    if(EXTRA_CXX_FLAGS_DEBUG AND CMAKE_BUILD_TYPE MATCHES Debug)
        add_compile_options("${EXTRA_CXX_FLAGS_DEBUG}")
    endif()

    # Extra Flags for the release builds with the C++ compiler.
    if(EXTRA_CXX_FLAGS_RELEASE AND CMAKE_BUILD_TYPE MATCHES RELEASE)
        add_compile_options("${EXTRA_CXX_FLAGS_RELEASE}")
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
## add_component( COMPONENT_NAME <name> DEFAULT_STATE [ON/OFF] )
##
## Adds a project component to the build.
##
## Adds a component to the build given the component's name and default state
## (ON/OFF). This macro also adds an "option" so that the user can control,
## which components to build.
##------------------------------------------------------------------------------
macro(add_component)

    set(options)
    set(singleValueArgs COMPONENT_NAME DEFAULT_STATE )
    set(multiValueArgs)

    ## parse the arugments to the macro
    cmake_parse_arguments(arg
         "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

    ## setup a cmake vars to capture sources added via our macros
    set("${arg_COMPONENT_NAME}_ALL_SOURCES" CACHE PATH "" FORCE)

    ## adds an option so that the user can control whether to build this
    ## component.
    option( ENABLE_${arg_COMPONENT_NAME}
            "Enables ${arg_component_name}"
            ${arg_DEFAULT_STATE})

    if ( ENABLE_${arg_COMPONENT_NAME} )
        add_subdirectory( ${arg_COMPONENT_NAME} )
    endif()

endmacro(add_component)

##------------------------------------------------------------------------------
## add_target_definitions(TO <target> TARGET_DEFINITIONS [FOO BAR ...])
##
## Adds pre-processor definitions to the given target.
##
## Adds pre-processor definitions to a particular. This macro provides very
## similar functionality to cmake's native "add_definitions" command, but,
## it provides more fine-grained scoping for the compile definitions on a
## per target basis. Given a list of definitions, e.g., FOO and BAR, this macro
## adds compiler definitions to the compiler command for the given target, i.e.,
## it will pass -DFOO and -DBAR.
##
## The supplied target must be added via add_executable() or add_library() or
## with the corresponding make_executable() and make_library() macros.
##
## Note, the list of target definitions *SHOULD NOT* include the "-D" flag. This
## flag is added internally by cmake.
##------------------------------------------------------------------------------
macro(add_target_definitions)

   set(options)
   set(singleValueArgs TO)
   set(multiValueArgs TARGET_DEFINITIONS)

   ## parse the arguments to the macro
   cmake_parse_arguments(arg
        "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

  get_target_property(defs ${arg_TO} COMPILE_DEFINITIONS)
  if (defs MATCHES "NOTFOUND")
    set(defs "")
  endif ()
  foreach (def ${defs} ${arg_TARGET_DEFINITIONS})
    list(APPEND deflist ${def})
  endforeach ()
  set_target_properties(${arg_TO} PROPERTIES COMPILE_DEFINITIONS "${deflist}")

endmacro(add_target_definitions)

##------------------------------------------------------------------------------
## make_library( LIBRARY_NAME <libname> LIBRARY_SOURCES [source1 [source2 ...]]
##               [WITH_MPI] [WITH_OPENMP])
##
## Adds a library to the project composed by the given source files.
##
## Adds a library target, called <libname>, to be built from the given sources.
## This macro internally checks if the global option "BUILD_SHARED_LIBS" is
## ON, in which case, it will create a shared library. By default, a static
## library is generated.
##
## Optionally, "WITH_MPI" can be supplied as an argument. When this argument is
## supplied, the MPI include directory will be added to the compiler command
## and the -DUSE_MPI, as well as other compiler flags, will be included to the
## compiler definition.
##
## Optionally, "WITH_OPENMP" can be supplied as an argument. When this argument is
## supplied, the openmp compiler flag will be added to the compiler command
## and the -DUSE_MPI, will be included to the compiler definition.
##------------------------------------------------------------------------------
macro(make_library)

   set(options WITH_MPI WITH_OPENMP)
   set(singleValueArgs LIBRARY_NAME)
   set(multiValueArgs LIBRARY_SOURCES)

   ## parse the arguments
   cmake_parse_arguments(arg
        "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN} )

   ## sanity check MPI
   if ( ${arg_WITH_MPI} AND NOT ${ENABLE_MPI} )
      message( FATAL_ERROR "Building an MPI library, but MPI is disabled!" )
   endif()

   ## sanity check OpenMP
   if ( ${arg_WITH_OPENMP} AND NOT ${ENABLE_OPENMP} )
      message( FATAL_ERROR "Building an OpenMP library, but OpenMP is disabled!" )
   endif()

   if ( BUILD_SHARED_LIBS )
      add_library(${arg_LIBRARY_NAME} SHARED ${arg_LIBRARY_SOURCES})
   else()
      add_library(${arg_LIBRARY_NAME} STATIC ${arg_LIBRARY_SOURCES})
   endif()

   if ( ${arg_WITH_MPI} )

      add_target_definitions( TO ${arg_LIBRARY_NAME}
                              TARGET_DEFINITIONS USE_MPI )

      target_include_directories( ${arg_LIBRARY_NAME} PUBLIC
                                  ${MPI_C_INCLUDE_PATH} )

      if(NOT "${MPI_CXX_COMPILE_FLAGS}" STREQUAL "")
            set_target_properties( ${arg_LIBRARY_NAME} PROPERTIES COMPILE_FLAGS
                                   ${MPI_CXX_COMPILE_FLAGS} )
      endif()

      if(NOT "${MPI_CXX_LINK_FLAGS}" STREQUAL "")
            set_target_properties( ${arg_LIBRARY_NAME} PROPERTIES LINK_FLAGS
                                   ${MPI_CXX_LINK_FLAGS} )
      endif()

   endif()

   if ( ${arg_WITH_OPENMP} )

      add_target_definitions( TO ${arg_LIBRARY_NAME} TARGET_DEFINITIONS USE_OPENMP )

      set_target_properties( ${arg_LIBRARY_NAME} PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS} )

   endif()

   if ( ENABLE_CXX11 )
      ## Note, this requires cmake 3.1 and above
      set_property(TARGET ${arg_LIBRARY_NAME} PROPERTY CXX_STANDARD 11)
   endif()

   foreach(src ${srcs})
       if(IS_ABSOLUTE)
           list(APPEND "${PROJECT_NAME}_ALL_SOURCES" "${src}")
       else()
           list(APPEND "${PROJECT_NAME}_ALL_SOURCES"
                       "${CMAKE_CURRENT_SOURCE_DIR}/${src}")
       endif()
   endforeach()

   set( "${PROJECT_NAME}_ALL_SOURCES" "${${PROJECT_NAME}_ALL_SOURCES}"
        CACHE STRING "" FORCE )

endmacro(make_library)

##------------------------------------------------------------------------------
## make_executable(EXECUTABLE_SOURCE <source> DEPENDS_ON [dep1 ...] [WITH_MPI] [WITH_OPENMP])
##
## Adds an executable to the project.
##
## Adds an executable target, called <name>, where <name> corresponds to the
## the filename of the given <source> without the file extension, e.g., the
## executable of a source, called driver.cxx, will be "driver".
##
## In addition, the target will be linked with the given list of library
## dependencies.
##
## Optionally, "WITH_MPI" can be supplied as an argument for MPI executables.
## When the "WITH_MPI" argument is supplied, the executable will be linked with
## the MPI C libraries and the MPI includes as well as -DUSE_MPI and other flags
## will be added to the compiler command.
##
## Optionally, "WITH_OPENMP" can be supplied as an argument. When this argument is
## supplied, the openmp compiler flag will be added to the compiler command
## and the -DUSE_MPI, will be included to the compiler definition.
##------------------------------------------------------------------------------
macro(make_executable)

   set(options WITH_MPI WITH_OPENMP)
   set(singleValueArgs EXECUTABLE_NAME EXECUTABLE_SOURCE)
   set(multiValueArgs DEPENDS_ON)

   ## parse the arguments to the macro
   cmake_parse_arguments(arg
        "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

   ## sanity check MPI
   if ( ${arg_WITH_MPI} AND NOT ${ENABLE_MPI} )
      message( FATAL_ERROR "Building an MPI executable, but MPI is disabled!" )
   endif()

   ## sanity check OpenMP
   if ( ${arg_WITH_OPENMP} AND NOT ${ENABLE_OPENMP} )
      message( FATAL_ERROR "Building an OpenMP executable, but OpenMP is disabled!" )
   endif()

   # Use the supplied name for the executable (if given), otherwise use the source file's name
   if( NOT "${arg_EXECUTABLE_NAME}" STREQUAL "" )
     set(exe_name ${arg_EXECUTABLE_NAME})
   else()
     get_filename_component(exe_name ${arg_EXECUTABLE_SOURCE} NAME_WE)
   endif()
      
   add_executable( ${exe_name} ${arg_EXECUTABLE_SOURCE} )
   target_link_libraries( ${exe_name} "${arg_DEPENDS_ON}" )

   if ( ENABLE_CXX11 )
     ## Note, this requires cmake 3.1 and above
     set_property(TARGET ${exe_name} PROPERTY CXX_STANDARD 11)
   endif()

   if ( ${arg_WITH_MPI} )
      add_target_definitions( TO ${exe_name}
                              TARGET_DEFINITIONS USE_MPI )

      target_include_directories( ${exe_name} PUBLIC
                                  ${MPI_C_INCLUDE_PATH} )

      if(NOT "${MPI_CXX_COMPILE_FLAGS}" STREQUAL "")
            set_target_properties( ${exe_name} PROPERTIES COMPILE_FLAGS
                                   ${MPI_CXX_COMPILE_FLAGS} )
      endif()

      if(NOT "${MPI_CXX_LINK_FLAGS}" STREQUAL "")
            set_target_properties( ${exe_name} PROPERTIES LINK_FLAGS
                                   ${MPI_CXX_LINK_FLAGS} )
      endif()

      target_link_libraries( ${exe_name} ${MPI_C_LIBRARIES})
   endif()

   if ( ${arg_WITH_OPENMP} )

      add_target_definitions( TO ${exe_name} TARGET_DEFINITIONS USE_OPENMP )

      set_target_properties( ${exe_name} PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS} )
      set_target_properties( ${exe_name} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS} )
      
   endif()

   if(IS_ABSOLUTE)
       list(APPEND "${PROJECT_NAME}_ALL_SOURCES" "${arg_EXECUTABLE_SOURCE}")
   else()
       list(APPEND "${PROJECT_NAME}_ALL_SOURCES"
                   "${CMAKE_CURRENT_SOURCE_DIR}/${arg_EXECUTABLE_SOURCE}")
   endif()

   set( "${PROJECT_NAME}_ALL_SOURCES" "${${PROJECT_NAME}_ALL_SOURCES}"
        CACHE STRING "" FORCE )

    if ( ENABLE_CXX11 )
      ## Note, this requires cmake 3.1 and above
      set_property(TARGET ${exe_name} PROPERTY CXX_STANDARD 11)
    endif()

endmacro(make_executable)

##------------------------------------------------------------------------------
## add_gtest( TEST_SOURCE testX.cxx DEPENDS_ON [dep1 [dep2 ...]] )
##
## Adds a google test to the project.
##------------------------------------------------------------------------------
macro(add_gtest)

   set(options)
   set(singleValueArgs TEST_SOURCE)
   set(multiValueArgs DEPENDS_ON)

   ## parse the arguments to the macro
   cmake_parse_arguments(arg
        "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN} )

   get_filename_component(test_name ${arg_TEST_SOURCE} NAME_WE)
   add_executable( ${test_name} ${arg_TEST_SOURCE} )
   target_include_directories(${test_name} PRIVATE "${GTEST_INCLUDES}")
   target_link_libraries( ${test_name} "${GTEST_LIBS}" )
   target_link_libraries( ${test_name} "${arg_DEPENDS_ON}" )

    if ( ENABLE_CXX11 )
      ## Note, this requires cmake 3.1 and above
      set_property(TARGET ${test_name} PROPERTY CXX_STANDARD 11)
    endif()

    add_test( NAME ${test_name}
              COMMAND ${test_name}
              WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
              )

    # add any passed source files to the running list for this project
    if(IS_ABSOLUTE)
        list(APPEND "${PROJECT_NAME}_ALL_SOURCES" "${arg_TEST_SOURCE}")
    else()
          list(APPEND "${PROJECT_NAME}_ALL_SOURCES"
                      "${CMAKE_CURRENT_SOURCE_DIR}/${arg_TEST_SOURCE}")
    endif()


    set("${PROJECT_NAME}_ALL_SOURCES" "${${PROJECT_NAME}_ALL_SOURCES}"
        CACHE STRING "" FORCE )

endmacro(add_gtest)

##------------------------------------------------------------------------------
## add_catch_test( TEST_SOURCE testX.cxx DEPENDS_ON [dep1 [dep2 ...]] )
##
## Adds a catch test to the project.
##------------------------------------------------------------------------------
macro(add_catch_test)

   set(options)
   set(singleValueArgs TEST_SOURCE)
   set(multiValueArgs DEPENDS_ON)

   ## parse the arguments to the macro
   cmake_parse_arguments(arg
        "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN} )

   get_filename_component(test_name ${arg_TEST_SOURCE} NAME_WE)
   add_executable( ${test_name} ${arg_TEST_SOURCE} )
   target_link_libraries( ${test_name} "${arg_DEPENDS_ON}" )

    if ( ENABLE_CXX11 )
      ## Note, this requires cmake 3.1 and above
      set_property(TARGET ${test_name} PROPERTY CXX_STANDARD 11)
    endif()

    add_test( NAME ${test_name}
              COMMAND ${test_name}
              WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
              )

    # add any passed source files to the running list for this project
    if(IS_ABSOLUTE)
        list(APPEND "${PROJECT_NAME}_ALL_SOURCES" "${arg_TEST_SOURCE}")
    else()
        list(APPEND "${PROJECT_NAME}_ALL_SOURCES"
                    "${CMAKE_CURRENT_SOURCE_DIR}/${arg_TEST_SOURCE}")
    endif()

    set("${PROJECT_NAME}_ALL_SOURCES" "${${PROJECT_NAME}_ALL_SOURCES}"
        CACHE STRING "" FORCE )

endmacro(add_catch_test)


##------------------------------------------------------------------------------
## - Builds and adds a fortran based test.
##
## add_fortran_test( TEST_SOURCE testX.f DEPENDS_ON dep1 dep2... )
##------------------------------------------------------------------------------
macro(add_fortran_test)
    # only add the test if fortran is enabled
    if(ENABLE_FORTRAN)
       set(options)
       set(singleValueArgs TEST_SOURCE)
       set(multiValueArgs DEPENDS_ON)

       ## parse the arguments to the macro
       cmake_parse_arguments(arg
            "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN} )

       get_filename_component(test_name ${arg_TEST_SOURCE} NAME_WE)
       add_executable( ${test_name} fortran_driver.cpp ${arg_TEST_SOURCE} )
       target_link_libraries( ${test_name} "${arg_DEPENDS_ON}" )

        set_target_properties(${test_name}  PROPERTIES Fortran_FORMAT "FREE")

        add_test( NAME ${test_name}
                  COMMAND ${test_name}
                  WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
                  )
       #TODO: we aren't tracking / grouping fortran sources. 
   endif()
endmacro(add_fortran_test)


##------------------------------------------------------------------------------
## copy_headers_target( <proj> <hdrs> <dest> )
##
## Adds a custom "copy_headers" target for the given project
##
## Adds a custom target, <copy_headers_proj>, for the given project. The role
## of this target is to copy the given list of headers, <hdrs>, to the
## destination directory <dest>.
##
## This macro is used to copy the header of each component in to the build
## space, under an "includes" directory.
##------------------------------------------------------------------------------
macro(copy_headers_target proj hdrs dest)

add_custom_target(copy_headers_${proj}
     COMMAND ${CMAKE_COMMAND}
             -DHEADER_INCLUDES_DIRECTORY=${dest}
             -DLIBHEADERS="${hdrs}"
             -P ${CMAKE_MODULE_PATH}/copy_headers.cmake

     DEPENDS
        ${hdrs}

     WORKING_DIRECTORY
        ${PROJECT_SOURCE_DIR}

     COMMENT
        "copy headers"
     )

     # add any passed source files to the running list for this project
     foreach(hdr ${hdrs})
         if(IS_ABSOLUTE)
             list(APPEND "${PROJECT_NAME}_ALL_SOURCES" "${hdr}")
         else()
             list(APPEND "${PROJECT_NAME}_ALL_SOURCES"
                         "${CMAKE_CURRENT_SOURCE_DIR}/${hdr}")
         endif()
     endforeach()

     set("${PROJECT_NAME}_ALL_SOURCES" "${${PROJECT_NAME}_ALL_SOURCES}"
         CACHE STRING "" FORCE )

endmacro(copy_headers_target)

