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
## blt_add_component( COMPONENT_NAME <name> DEFAULT_STATE [ON/OFF] )
##
## Adds a project component to the build.
##
## Adds a component to the build given the component's name and default state
## (ON/OFF). This macro also adds an "option" so that the user can control,
## which components to build.
##------------------------------------------------------------------------------
macro(blt_add_component)

    set(options)
    set(singleValueArgs COMPONENT_NAME DEFAULT_STATE )
    set(multiValueArgs)

    # Parse the arguments to the macro
    cmake_parse_arguments(arg
         "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

    # Setup a cmake vars to capture sources added via our macros
    set("${arg_COMPONENT_NAME}_ALL_SOURCES" CACHE PATH "" FORCE)

    # Adds an option so that the user can control whether to build this
    # component.
    # convert the component name to capitals for the ENABLE option.
    string(TOUPPER ${arg_COMPONENT_NAME} COMPONENT_NAME_CAPITALIZED)

    option( ENABLE_${COMPONENT_NAME_CAPITALIZED}
            "Enables ${arg_component_name}"
            ${arg_DEFAULT_STATE})

    if ( ENABLE_${COMPONENT_NAME_CAPITALIZED} )
        add_subdirectory( ${arg_COMPONENT_NAME} )
    endif()

endmacro(blt_add_component)


##------------------------------------------------------------------------------
## blt_add_target_definitions(TO <target> TARGET_DEFINITIONS [FOO [BAR ...]])
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
## with the corresponding blt_add_executable() and blt_add_library() macros.
##
## Note, the list of target definitions *SHOULD NOT* include the "-D" flag. This
## flag is added internally by cmake.
##------------------------------------------------------------------------------
macro(blt_add_target_definitions)

    set(options)
    set(singleValueArgs TO)
    set(multiValueArgs TARGET_DEFINITIONS)

    # Parse the arguments to the macro
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

endmacro(blt_add_target_definitions)


##------------------------------------------------------------------------------
## blt_register_library( NAME <libname>
##                       INCLUDES [include1 [include2 ...]] 
##                       FORTRAN_MODULES [ path1 [ path2 ..]]
##                       LIBRARIES [lib1 [lib2 ...]] )
##
## Registers a library to the project to ease use in other blt macro calls.
##
## Stores information about a library in a specific way that is easily recalled
## in other macros.  For example, after registering gtest, you can add gtest to
## the DEPENDS_ON in your blt_add_executable call and it will add the INCLUDES
## and LIBRARIES to that executable.
##
## This does not actually build the library.  This is strictly to ease use after
## discovering it on your system or building it yourself inside your project.
##
## Output variables (name = "foo"):
##  BLT_FOO_INCLUDES
##  BLT_FOO_FORTRAN_MODULES
##  BLT_FOO_LIBRARIES
##------------------------------------------------------------------------------
macro(blt_register_library)

    set(singleValueArgs NAME )
    set(multiValueArgs INCLUDES FORTRAN_MODULES LIBRARIES)

    ## parse the arguments
    cmake_parse_arguments(arg
        "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN} )

    string(TOUPPER ${arg_NAME} uppercase_name)

    if( arg_INCLUDES )
        set(BLT_${uppercase_name}_INCLUDES ${arg_INCLUDES})
    endif()

    if( arg_FORTRAN_MODULES )
        set(BLT_${uppercase_name}_FORTRAN_MODULES ${arg_INCLUDES})
    endif()

    if( arg_LIBRARIES )
        set(BLT_${uppercase_name}_LIBRARIES ${arg_LIBRARIES})
    endif()

endmacro(blt_register_library)


##------------------------------------------------------------------------------
## blt_add_library( NAME <libname>
##                  SOURCES [source1 [source2 ...]]
##                  HEADERS [header1 [header2 ...]]
##                  DEPENDS_ON [dep1 ...] 
##                  USE_OPENMP <TRUE or FALSE (default)> )
##
## Adds a library to the project composed by the given source files.
##
## Adds a library target, called <libname>, to be built from the given sources.
## This macro internally checks if the global option "ENABLE_SHARED_LIBS" is
## ON, in which case, it will create a shared library. By default, a static
## library is generated.
##
## If given a HEADERS argument, it creates a "blt_copy_headers_target" for 
## this library and installs them in the include/<component name> folder.
## 
## If given a DEPENDS_ON argument, it will add the necessary includes and 
##  libraries if they are already registered with blt_register_library.  If 
## not it will add them as a cmake target dependency.
##
## In addition, this macro will add the associated dependencies to the given
## library target. Specifically, it will add a dependency to the library's
## "blt_copy_headers_target" and adds a dependency for each DEPENDS_ON target.
##
## Optionally, "USE_OPENMP" can be supplied as a boolean argument. When this 
## argument is supplied, the openmp compiler flag will be added to the compiler 
## command and the -DUSE_OPENMP, will be included to the compiler definition.
##------------------------------------------------------------------------------
macro(blt_add_library)

    set(singleValueArgs NAME USE_OPENMP)
    set(multiValueArgs SOURCES HEADERS DEPENDS_ON)

    ## parse the arguments
    cmake_parse_arguments(arg
        "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN} )

    # Check for the variable-based options for OpenMP and sanity check
    if( NOT DEFINED arg_USE_OPENMP )
        set(arg_USE_OPENMP FALSE)
    endif()

    if ( ENABLE_SHARED_LIBS )
        add_library( ${arg_NAME} SHARED ${arg_SOURCES} ${arg_HEADERS} )
    else()
        add_library( ${arg_NAME} STATIC ${arg_SOURCES} ${arg_HEADERS} )
    endif()

    if ( arg_HEADERS )
        blt_copy_headers_target( ${arg_NAME} "${arg_HEADERS}"
                                ${HEADER_INCLUDES_DIRECTORY}/${PROJECT_NAME})
    endif()

    # Must tell fortran where to look for modules
    # CMAKE_Fortran_MODULE_DIRECTORY is the location of generated modules
    foreach (_file ${arg_SOURCES})
        get_source_file_property(_lang ${_file} LANGUAGE)
        if(_lang STREQUAL Fortran)
            set(_have_fortran TRUE)
        endif()
    endforeach()
    if(_have_fortran)
        include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})
    endif()

    blt_setup_target( NAME ${arg_NAME}
                      DEPENDS_ON ${arg_DEPENDS_ON} )

    # Handle MPI 
    blt_setup_mpi_target( BUILD_TARGET ${arg_NAME} )

    # Handle OpenMP
    blt_setup_openmp_target( BUILD_TARGET ${arg_NAME}
                              USE_OPENMP ${arg_USE_OPENMP} )

    # Update project sources                     
    blt_update_project_sources( TARGET_SOURCES ${arg_SOURCES} ${arg_HEADERS})
   
endmacro(blt_add_library)


##------------------------------------------------------------------------------
## blt_add_executable( NAME <name>
##                     SOURCES [source1 [source2 ...]]
##                     DEPENDS_ON [dep1 [dep2 ...]]
##                     OUTPUT_DIR [dir]
##                     USE_OPENMP < TRUE or FALSE (default)>)
##
## Adds an executable target, called <name>.
##
## If given a DEPENDS_ON argument, it will add the necessary includes and 
## libraries if they are already registered with blt_register_library.  If
## not it will add them as a cmake target dependency.
##
## Optionally, "USE_OPENMP" can be supplied as a boolean argument. When this 
## argument is supplied, the openmp compiler flag will be added to the compiler 
## command and the -DUSE_OPENMP, will be included to the compiler definition.
##
## The OUTPUT_DIR is used to control the build output directory of this 
## executable. This is used to overwrite the default bin directory.
##
## If the first entry in SOURCES is a Fortran source file, the fortran linker 
## is used. (via setting the CMake target property LINKER_LANGUAGE to Fortran )
##
##------------------------------------------------------------------------------
macro(blt_add_executable)

    set(options )
    set(singleValueArgs NAME OUTPUT_DIR USE_OPENMP)
    set(multiValueArgs SOURCES DEPENDS_ON)

    # Parse the arguments to the macro
    cmake_parse_arguments(arg
        "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

    # Check for the variable-based options for OpenMP and sanity check
    if ( NOT DEFINED arg_USE_OPENMP )
        set(arg_USE_OPENMP FALSE)
    endif()

    if( "${arg_NAME}" STREQUAL "" )
        message(FATAL_ERROR "Must specify executable name with argument NAME <name>")
    endif()

    add_executable( ${arg_NAME} ${arg_SOURCES} )

    # CMake wants to load with C++ if any of the libraries are C++.
    # Force to load with Fortran if the first file is Fortran.
    list(GET arg_SOURCES 0 _first)
    get_source_file_property(_lang ${_first} LANGUAGE)
    if(_lang STREQUAL Fortran)
        set_target_properties( ${test_name} PROPERTIES LINKER_LANGUAGE Fortran )
    endif()

    blt_setup_target(NAME ${arg_NAME}
                     DEPENDS_ON ${arg_DEPENDS_ON} )

    # Handle MPI
    blt_setup_mpi_target( BUILD_TARGET ${arg_NAME} )

    # Handle OpenMP
    blt_setup_openmp_target( BUILD_TARGET ${arg_NAME} 
                             USE_OPENMP ${arg_USE_OPENMP} ) 

    # Set output directory
    if ( arg_OUTPUT_DIR )
        set_target_properties(${arg_NAME} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY ${arg_OUTPUT_DIR} )
    endif()

    # Update project sources
    blt_update_project_sources( TARGET_SOURCES ${arg_SOURCES} )

endmacro(blt_add_executable)


##------------------------------------------------------------------------------
## blt_add_test( NAME [name] COMMAND [command] NUM_PROCS [n] )
##
## Adds a cmake test to the project.
##
## NAME is used for the name that CTest reports with.
##
## COMMAND is the command line that will be used to run the test.  This will have
## the RUNTIME_OUTPUT_DIRECTORY prepended to it to fully qualify the path.
##
## NUM_PROCS indicates this is an MPI test and how many processors to use. The
## command line will use MPIEXEC and MPIXEC_NUMPROC_FLAG to create the mpi run line.
## These should be defined in your host-config specific to your platform.
##------------------------------------------------------------------------------
macro(blt_add_test)

    set(options )
    set(singleValueArgs NAME NUM_PROCS)
    set(multiValueArgs COMMAND)

    # Parse the arguments to the macro
    cmake_parse_arguments(arg
        "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN} )

    if ( NOT DEFINED arg_NAME )
        message(FATAL_ERROR "NAME is a required parameter to blt_add_test")
    endif()

    if ( NOT DEFINED arg_COMMAND )
        message(FATAL_ERROR "COMMAND is a required parameter to blt_add_test")
    endif()

    # Generate command
    if ( NOT TARGET ${arg_NAME} )
        # Handle case of running multiple tests against one executable, 
        # the NAME will not be the target
        list(GET arg_COMMAND 0 executable)
        get_target_property(runtime_output_directory ${executable} RUNTIME_OUTPUT_DIRECTORY )
    else()
        get_target_property(runtime_output_directory ${arg_NAME} RUNTIME_OUTPUT_DIRECTORY )
    endif()
    set(test_command ${runtime_output_directory}/${arg_COMMAND} )

    # Handle mpi
    if ( ${arg_NUM_PROCS} )
        set(test_command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${arg_NUM_PROCS} ${test_command} )
    endif()

    add_test( NAME ${arg_NAME}
              COMMAND ${test_command} 
              )

endmacro(blt_add_test)


##------------------------------------------------------------------------------
## blt_add_benchmark( TEST_SOURCE testX.cxx TEST_ARGS <commandLineArguments> 
##                DEPENDS_ON [dep1 [dep2 ...]] )
##
## Adds a (google) benchmark test to the project.
## TEST_ARGS is a string containing a space delimited set of command line 
## arguments for the test
##        
##  blt_add_benchmark( 
##          TEST_SOURCE slamBench.cpp 
##          TEST_ARGS "--benchmark_min_time=0.0 --v=3 --benchmark_format=json"
##          DEPENDS_ON slic slam )
##------------------------------------------------------------------------------
macro(blt_add_benchmark)

   if(ENABLE_BENCHMARKS)
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
      blt_update_project_sources( TARGET_SOURCES ${arg_TEST_SOURCE} )
      
   endif(ENABLE_BENCHMARKS)
endmacro(blt_add_benchmark)

##------------------------------------------------------------------------------
## blt_append_custom_compiler_flag( 
##                    FLAGS_VAR flagsVar     (required)
##                    DEFAULT   defaultFlag  (optional)
##                    GNU       gnuFlag      (optional)
##                    CLANG     clangFlag    (optional)
##                    INTEL     intelFlag    (optional)
##                    XL        xlFlag       (optional)
##                    MSVC      msvcFlag     (optional)
## )
##
## Appends compiler-specific flags to a given variable of flags
##
## If a custom flag is given for the current compiler, we use that,
## Otherwise, we will use the DEFAULT flag (if present)
##------------------------------------------------------------------------------
macro(blt_append_custom_compiler_flag)

   set(options)
   set(singleValueArgs FLAGS_VAR DEFAULT GNU CLANG INTEL XL MSVC)
   set(multiValueArgs)

   # Parse the arguments
   cmake_parse_arguments(arg
        "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN} )


   # Sanity check for required keywords
   if(NOT DEFINED arg_FLAGS_VAR)
      message( FATAL_ERROR "append_custom_compiler_flag macro requires FLAGS_VAR keyword and argument." )
   endif()


   # Set the desired flags based on the compiler family   
   if(DEFINED arg_CLANG AND COMPILER_FAMILY_IS_CLANG )
      set (${arg_FLAGS_VAR} "${${arg_FLAGS_VAR}} ${arg_CLANG} " )
      
   elseif(DEFINED arg_XL AND COMPILER_FAMILY_IS_XL )
      set (${arg_FLAGS_VAR} "${${arg_FLAGS_VAR}} ${arg_XL} " )
      
   elseif(DEFINED arg_INTEL AND COMPILER_FAMILY_IS_INTEL )
      set (${arg_FLAGS_VAR} "${${arg_FLAGS_VAR}} ${arg_INTEL} " )

   elseif(DEFINED arg_GNU AND COMPILER_FAMILY_IS_GNU )
      set (${arg_FLAGS_VAR} "${${arg_FLAGS_VAR}} ${arg_GNU} " )
      
   elseif(DEFINED arg_MSVC AND COMPILER_FAMILY_IS_MSVC )
      set (${arg_FLAGS_VAR} "${${arg_FLAGS_VAR}} ${arg_MSVC} " )
      
   elseif(DEFINED arg_DEFAULT)
      set (${arg_FLAGS_VAR} "${${arg_FLAGS_VAR}} ${arg_DEFAULT} ")
      
   endif()   

   #message(STATUS "After append -- ${arg_FLAGS_VAR} --- ${${arg_FLAGS_VAR}} ")

endmacro(blt_append_custom_compiler_flag)

