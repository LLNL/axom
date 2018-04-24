#------------------------------------------------------------------------------
# Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory.
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
#------------------------------------------------------------------------------

##------------------------------------------------------------------------------
## axom_add_code_checks( PREFIX     <Prefix used for created targets>
##                       EXCLUDES   [path1 [path2 ...]])
##
## Adds code checks to all source files under this directory.
##
## PREFIX is used in the creation of all the underlying targets. For example:
## <PREFIX>_uncrustify_check.
##
## EXCLUDES is used to exclude any files from the code checks. It is done with
## a simple CMake reg exp MATCHES check.
##
##------------------------------------------------------------------------------
macro(axom_add_code_checks)

    set(options)
    set(singleValueArgs PREFIX )
    set(multiValueArgs EXCLUDES )

    # Parse the arguments to the macro
    cmake_parse_arguments(arg
         "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

    set(_all_sources)
    file(GLOB_RECURSE _all_sources
         "*.cpp" "*.hpp" "*.cxx" "*.hxx" "*.cc" "*.c" "*.h" "*.hh"
         "*.F" "*.f" "*.f90" "*.F90")

    # Check for excludes
    if (NOT DEFINED arg_EXCLUDES)
        set(_sources ${_all_sources})
    else()
        set(_sources)
        foreach(_source ${_all_sources})
            set(_to_be_excluded FALSE)
            foreach(_exclude ${arg_EXCLUDES})
                if (${_source} MATCHES ${_exclude})
                    set(_to_be_excluded TRUE)
                    break()
                endif()
            endforeach()

            if (NOT ${_to_be_excluded})
                list(APPEND _sources ${_source})
            endif()
        endforeach()
    endif()

    blt_add_code_checks(PREFIX    ${arg_PREFIX}
                        SOURCES   ${_sources}
                        UNCRUSTIFY_CFG_FILE ${PROJECT_SOURCE_DIR}/src/uncrustify.cfg)

endmacro(axom_add_code_checks)


##------------------------------------------------------------------------------
## axom_add_component( COMPONENT_NAME <name> DEFAULT_STATE [ON/OFF] )
##
## Adds a project component to the build.
##
## Adds a component to the build given the component's name and default state
## (ON/OFF). This macro also adds an "option" so that the user can control,
## which components to build.
##------------------------------------------------------------------------------
macro(axom_add_component)

    set(options)
    set(singleValueArgs COMPONENT_NAME DEFAULT_STATE )
    set(multiValueArgs)

    # Parse the arguments to the macro
    cmake_parse_arguments(arg
         "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

    # Adds an option so that the user can control whether to build this
    # component.
    # convert the component name to capitals for the ENABLE option.
    string(TOUPPER ${arg_COMPONENT_NAME} COMPONENT_NAME_CAPITALIZED)
    string(TOLOWER ${arg_COMPONENT_NAME} COMPONENT_NAME_LOWERED)

    option( ENABLE_${COMPONENT_NAME_CAPITALIZED}
            "Enables ${arg_component_name}"
            ${arg_DEFAULT_STATE})

    if ( ENABLE_${COMPONENT_NAME_CAPITALIZED} )
        add_subdirectory( ${arg_COMPONENT_NAME} )
    endif()

    unset(COMPONENT_NAME_CAPITALIZED)
    unset(COMPONENT_NAME_LOWERED)
endmacro(axom_add_component)


##------------------------------------------------------------------------------
## convert_to_native_escaped_file_path( path output )
##
## This macro converts a cmake path to a platform specific string literal
## usable in C++.  (For example, on windows C:/Path will be come C:\\Path)
##------------------------------------------------------------------------------

macro(convert_to_native_escaped_file_path path output)
    file(TO_NATIVE_PATH ${path} ${output})
    string(REPLACE "\\" "\\\\"  ${output} "${${output}}")
endmacro()



##------------------------------------------------------------------------------
## axom_check_code_compiles
## 
## This macro checks if a snippet of C++ code with a main function compiles.
##
## SOURCE_STRING The source snippet to compile. 
## Must be a valid C++ program with a main() function.
##
## CODE_COMPILES A boolean variable the contains the compilation result
##------------------------------------------------------------------------------
macro(axom_check_code_compiles)

    set(options)
    set(singleValueArgs CODE_COMPILES)
    set(multiValueArgs SOURCE_STRING )

    # Parse the arguments to the macro
    cmake_parse_arguments(arg
         "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})


    # Check the arguments
    if(NOT DEFINED arg_SOURCE_STRING)
        message(FATAL_ERROR "SOURCE_STRING is a required parameter for axom_check_code_compiles macro")
    endif()
    if(NOT DEFINED arg_CODE_COMPILES)
        message(FATAL_ERROR "CODE_COMPILES is a required parameter for axom_check_code_compiles macro")
    endif()    

    # message(STATUS "axom_check_code_compiles source string: ${arg_SOURCE_STRING}")

    # Write string as file, try to compile and then remove file
    string(RANDOM LENGTH 5 _rand)
    set(_fname ${CMAKE_CURRENT_BINARY_DIR}/_axomCheckCompiles${_rand}.cpp)
    file(WRITE ${_fname} "${arg_SOURCE_STRING}")
    try_compile(${arg_CODE_COMPILES}
                ${CMAKE_CURRENT_BINARY_DIR}/CMakeTmp      
                SOURCES ${_fname}
                OUTPUT_VARIABLE _res)
    file(REMOVE ${_fname})

    # message(STATUS "Output of axom_check_code_compiles: \n ${_res}")

    # clear the variables set within the macro
    set(_fname)
    set(_res)

endmacro(axom_check_code_compiles)
