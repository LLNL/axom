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

include(PrivateMacros)

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

    ## parse the arguments to the macro
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
##               DEPENDS_ON [dep1 ...] 
##               USE_OPENMP <TRUE or FALSE (default)> )
##
## Adds a library to the project composed by the given source files.
##
## Adds a library target, called <libname>, to be built from the given sources.
## This macro internally checks if the global option "BUILD_SHARED_LIBS" is
## ON, in which case, it will create a shared library. By default, a static
## library is generated.
##
## In addition, this macro will add the associated dependencies to the given
## library target. Specifically, it will add a dependency to the library's
## "copy_headers_target" if it exists and has been defined before the call to
## "make_library", as well as, the corresponding "copy_headers_target" of each
## of the supplied dependencies.
##
## Optionally, "USE_OPENMP" can be supplied as a boolean argument. When this 
## argument is supplied, the openmp compiler flag will be added to the compiler 
## command and the -DUSE_OPENMP, will be included to the compiler definition.
##------------------------------------------------------------------------------
macro(make_library)

   set(singleValueArgs LIBRARY_NAME USE_OPENMP)
   set(multiValueArgs LIBRARY_SOURCES DEPENDS_ON)

   ## parse the arguments
   cmake_parse_arguments(arg
        "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN} )

   # Check for the variable-based options for OpenMP and sanity check
   if(NOT DEFINED arg_USE_OPENMP)
      set(arg_USE_OPENMP FALSE)
   endif()

   if ( BUILD_SHARED_LIBS )
      add_library(${arg_LIBRARY_NAME} SHARED ${arg_LIBRARY_SOURCES})
   else()
      add_library(${arg_LIBRARY_NAME} STATIC ${arg_LIBRARY_SOURCES})
   endif()

    ## handle MPI 
   setup_mpi_target( BUILD_TARGET ${arg_LIBRARY_NAME} )
   
   ## handle OpenMP
   setup_openmp_target( BUILD_TARGET ${arg_LIBRARY_NAME}
                        USE_OPENMP ${arg_USE_OPENMP} )
   
   ## update project sources                     
   update_project_sources( TARGET_SOURCES ${arg_LIBRARY_SOURCES})
   
   ## setup dependencies
   set(lib_header_target "copy_headers_${arg_LIBRARY_NAME}")
   if (TARGET ${lib_header_target})
      add_dependencies( ${arg_LIBRARY_NAME} ${lib_header_target})
   endif()

   foreach(dependency ${arg_DEPENDS_ON})
     
     if (TARGET ${dependency})
        target_link_libraries(${arg_LIBRARY_NAME} ${dependency})
     endif()
     
     set(header_target "copy_headers_${dependency}")
     if (TARGET ${header_target})
        add_dependencies( ${arg_LIBRARY_NAME} ${header_target} )
     endif()
     
   endforeach()

endmacro(make_library)

##------------------------------------------------------------------------------
## make_executable( EXECUTABLE_SOURCE <source> DEPENDS_ON [dep1 ...]
##                  USE_OPENMP < TRUE or FALSE (default)>
##                  [ADD_CTEST])
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
## Optionally, "USE_OPENMP" can be supplied as a boolean argument. When this 
## argument is supplied, the openmp compiler flag will be added to the compiler 
## command and the -DUSE_OPENMP, will be included to the compiler definition.
##
## Optionally, "ADD_CTEST" can be supplied as an argument. When this argument
## is supplied, the executable is added as a ctest with no command line options.
## If you need command line options for your ctest, add it manually with 
## ADD_TEST().
##------------------------------------------------------------------------------
macro(make_executable)

  set(options ADD_CTEST)
  set(singleValueArgs EXECUTABLE_NAME EXECUTABLE_SOURCE USE_OPENMP)
  set(multiValueArgs DEPENDS_ON)

  ## parse the arguments to the macro
  cmake_parse_arguments(arg
      "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

  # Check for the variable-based options for OpenMP and sanity check
  if ( NOT DEFINED arg_USE_OPENMP )
     set(arg_USE_OPENMP FALSE)
  endif()

   # Use the supplied name for the executable (if given), otherwise use the
   # source file's name
   if( NOT "${arg_EXECUTABLE_NAME}" STREQUAL "" )
     set(exe_name ${arg_EXECUTABLE_NAME})
   else()
     get_filename_component(exe_name ${arg_EXECUTABLE_SOURCE} NAME_WE)
   endif()

   add_executable( ${exe_name} ${arg_EXECUTABLE_SOURCE} )

   target_link_libraries(${exe_name} ${arg_DEPENDS_ON})

   ## Add library and header dependencies
   ##  Want to make this more general as in make_library above
   ## Problem -- we want to add dependencies that are not targets -- e.g. lib rt
   # foreach(dependency ${arg_DEPENDS_ON})
   #  if (TARGET ${dependency})
   #     target_link_libraries(${exe_name} ${dependency})
   #  endif()
   #  set(header_target "copy_headers_${dependency}")
   #  if (TARGET ${header_target})
   #     add_dependencies( ${exe_name} ${header_target} )
   #  endif()
   #endforeach()

    ## Handle MPI
   setup_mpi_target( BUILD_TARGET ${exe_name} )
   
   ## Handle OpenMP
   setup_openmp_target( BUILD_TARGET ${exe_name} USE_OPENMP ${arg_USE_OPENMP} ) 
   
   if ( ${arg_ADD_CTEST} )
     add_test( NAME ${exe_name}
               COMMAND ${exe_name}
               WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
               )
   endif()

   ## update project sources
   update_project_sources( TARGET_SOURCES ${arg_EXECUTABLE_SOURCE} )

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

   get_filename_component(test_name_base ${arg_TEST_SOURCE} NAME_WE)
   set(test_name ${test_name_base}_gtest)
   add_executable( ${test_name} ${arg_TEST_SOURCE} )
   target_include_directories(${test_name} PRIVATE "${GTEST_INCLUDES}")
   target_link_libraries( ${test_name} "${GTEST_LIBS}" )
   target_link_libraries( ${test_name} "${arg_DEPENDS_ON}" )

   add_test( NAME ${test_name}
             COMMAND ${test_name}
             WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
             )

   # add any passed source files to the running list for this project
   update_project_sources( TARGET_SOURCES ${arg_TEST_SOURCE} )

endmacro(add_gtest)


##------------------------------------------------------------------------------
## add_benchmark( TEST_SOURCE testX.cxx TEST_ARGS <commandLineArguments> 
##                DEPENDS_ON [dep1 [dep2 ...]] )
##
## Adds a (google) benchmark test to the project.
## TEST_ARGS is a string containing a space delimited set of command line 
## arguments for the test
##        
##  add_benchmark( 
##          TEST_SOURCE slamBench.cpp 
##          TEST_ARGS "--benchmark_min_time=0.0 --v=3 --benchmark_format=json"
##          DEPENDS_ON slic slam )
##------------------------------------------------------------------------------
macro(add_benchmark)

   if(ENABLE_BENCHMARK)
      set(options)
      set(singleValueArgs TEST_SOURCE TEST_ARGS)
      set(multiValueArgs DEPENDS_ON)

      ## parse the arguments to the macro
      cmake_parse_arguments(arg
           "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN} )

      get_filename_component(test_name_base ${arg_TEST_SOURCE} NAME_WE)
      set(test_name ${test_name_base}_benchmark)
      add_executable( ${test_name} ${arg_TEST_SOURCE} )
      target_include_directories(${test_name} PRIVATE "${GBENCHMARK_INCLUDES}")
      target_link_libraries( ${test_name} "${GBENCHMARK_LIBS}" )
      target_link_libraries( ${test_name} "${arg_DEPENDS_ON}" )

      # Add the command line arguments, if present
      set(test_command "${test_name}")
      if( arg_TEST_ARGS )
         separate_arguments(argList UNIX_COMMAND ${arg_TEST_ARGS})
         list(APPEND test_command ${argList})
      endif()
      
      # The 'CONFIGURATIONS Benchmark' line excludes benchmarks 
      # from the general list of tests
      add_test( NAME ${test_name}
                COMMAND ${test_command}
                CONFIGURATIONS Benchmark   
                WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
                )

      add_dependencies(run_benchmarks ${test_name})

      
      # add any passed source files to the running list for this project
      update_project_sources( TARGET_SOURCES ${arg_TEST_SOURCE} )
      
   endif(ENABLE_BENCHMARK)
endmacro(add_benchmark)

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

       get_filename_component(test_name_base ${arg_TEST_SOURCE} NAME_WE)
       set(test_name ${test_name_base}_ftest)
       add_executable( ${test_name} ${arg_TEST_SOURCE} )

       target_include_directories( ${test_name} PUBLIC 
                                   ${CMAKE_Fortran_MODULE_DIRECTORY} )
       target_link_libraries( ${test_name} "${arg_DEPENDS_ON}" )

       set_target_properties(${test_name} PROPERTIES LINKER_LANGUAGE Fortran)
       set_target_properties(${test_name} PROPERTIES Fortran_FORMAT "FREE")

        add_test( NAME ${test_name}
                  COMMAND ${test_name}
                  WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
                  )
       #TODO: we aren't tracking / grouping fortran sources.
    endif(ENABLE_FORTRAN)
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
             -P ${CMAKE_SOURCE_DIR}/cmake/copy_headers.cmake

     DEPENDS
        ${hdrs}

     WORKING_DIRECTORY
        ${PROJECT_SOURCE_DIR}

     COMMENT
        "copy headers"
     )

     # add any passed source files to the running list for this project
     foreach(hdr ${hdrs})
         if(IS_ABSOLUTE ${hdr})
             list(APPEND "${PROJECT_NAME}_ALL_SOURCES" "${hdr}")
         else()
             list(APPEND "${PROJECT_NAME}_ALL_SOURCES"
                         "${CMAKE_CURRENT_SOURCE_DIR}/${hdr}")
         endif()
     endforeach()

     set("${PROJECT_NAME}_ALL_SOURCES" "${${PROJECT_NAME}_ALL_SOURCES}"
         CACHE STRING "" FORCE )

endmacro(copy_headers_target)

