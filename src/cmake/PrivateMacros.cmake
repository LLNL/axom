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

## Internal CMake Macros


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

endmacro(copy_headers_target)



##------------------------------------------------------------------------------
## blt_setup_target( NAME [name] DEPENDS_ON [dep1 ...] )
##------------------------------------------------------------------------------
macro(blt_setup_target)
    set(options)
    set(singleValueArgs NAME)
    set(multiValueArgs DEPENDS_ON)

    # Parse the arguments
    cmake_parse_arguments(arg "${options}" "${singleValueArgs}" 
                        "${multiValueArgs}" ${ARGN} )
                        
    # Ensure that build target is supplied by the caller
    if ( NOT DEFINED arg_NAME )
        message( FATAL_ERROR "Must provide a NAME argument to the macro" )
    endif()

    # Add it's own copy headers target
    if (TARGET "copy_headers_${arg_NAME}")
        add_dependencies( ${arg_NAME} "copy_headers_${arg_NAME}")
    endif()

    # Add dependency's information
    foreach( dependency ${arg_DEPENDS_ON} )
        string(TOUPPER ${dependency} uppercase_dependency )

        if ( DEFINED BLT_${uppercase_dependency}_INCLUDES )
            target_include_directories( ${arg_NAME} PRIVATE
                ${BLT_${uppercase_dependency}_INCLUDES} )
        endif()

        if ( DEFINED BLT_${uppercase_dependency}_LIBRARIES )
            target_link_libraries( ${arg_NAME}
                ${BLT_${uppercase_dependency}_LIBRARIES} )
        else()
            target_link_libraries( ${arg_NAME} ${dependency} )
        endif()

        if (TARGET "copy_headers_${dependency}")
            add_dependencies( ${arg_NAME} "copy_headers_${dependency}" )
        endif()
    endforeach()

endmacro(blt_setup_target)

##------------------------------------------------------------------------------
## setup_mpi_target( BUILD_TARGET <target> )
##------------------------------------------------------------------------------
macro(setup_mpi_target)
  
  set(options)
  set(singleValueArgs BUILD_TARGET)
  set(multiValueArgs)
  
  ## parse the arguments
  cmake_parse_arguments(arg "${options}" "${singleValueArgs}" 
                            "${multiValueArgs}" ${ARGN} )
                            
  ## ensure that build target is supplied by the caller
  if ( NOT DEFINED arg_BUILD_TARGET )
    message( FATAL_ERROR "Must provide a BUILD_TARGET argument to the macro" )
  endif()
  
  if ( ${ENABLE_MPI} )
 
    add_target_definitions( TO ${arg_BUILD_TARGET} TARGET_DEFINITIONS USE_MPI )
    target_include_directories( ${arg_BUILD_TARGET} 
                                PUBLIC ${MPI_C_INCLUDE_PATH} )
    target_include_directories( ${arg_BUILD_TARGET} 
                                PUBLIC ${MPI_Fortran_INCLUDE_PATH} )

    if ( NOT "${MPI_C_COMPILE_FLAGS}" STREQUAL "")
       set_target_properties( ${arg_BUILD_TARGET} 
              PROPERTIES COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS} )
    endif()

    if ( NOT "${MPI_C_LINK_FLAGS}" STREQUAL "")
       set_target_properties( ${arg_BUILD_TARGET} 
              PROPERTIES LINK_FLAGS ${MPI_C_LINK_FLAGS} )
    endif()

    if ( NOT "${MPI_Fortran_LINK_FLAGS}" STREQUAL "" )
      set_target_properties( ${arg_BUILD_TARGET} 
             PROPERTIES LINK_FLAGS ${MPI_Fortran_LINK_FLAGS} )
    endif()
    
    target_link_libraries( ${arg_BUILD_TARGET} ${MPI_C_LIBRARIES} )
    target_link_libraries( ${arg_BUILD_TARGET} ${MPI_Fortran_LIBRARIES} )
  endif()

endmacro(setup_mpi_target)

##------------------------------------------------------------------------------
## setup_openmp_target( TARGET <target> USE_OPENMP <bool> )
##------------------------------------------------------------------------------
macro(setup_openmp_target)

  set(options)
  set(singleValueArgs BUILD_TARGET USE_OPENMP)
  set(multiValueArgs)
  
  ## parse the arguments
  cmake_parse_arguments(arg "${options}" "${singleValueArgs}" 
                            "${multiValueArgs}" ${ARGN} )
  
  if ( NOT DEFINED arg_BUILD_TARGET )
    message ( FATAL_ERROR "Must provide a BUILD_TARGET argument to the macro")
  endif()
  
  if ( NOT DEFINED arg_USE_OPENMP )
    message( FATAL_ERROR "Must provide an OpenMP boolean flag")
  endif()
  
  if ( ${arg_USE_OPENMP} AND NOT ${ENABLE_OPENMP} )
    message( FATAL_ERROR "Building an OpenMP library, but OpenMP is disabled!")
  endif()
  
  if ( ${arg_USE_OPENMP} )

    add_target_definitions( TO ${arg_BUILD_TARGET}
                            TARGET_DEFINITIONS USE_OPENMP )

    set_target_properties( ${arg_BUILD_TARGET}
                           PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS} )
    set_target_properties( ${arg_BUILD_TARGET}
                           PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS} )

  endif()

endmacro(setup_openmp_target)

##------------------------------------------------------------------------------
## update_project_sources( TARGET_SOURCES <souces> )
##------------------------------------------------------------------------------
macro(update_project_sources)

  set(options)
  set(singleValueArgs)
  set(multiValueArgs TARGET_SOURCES)
  
  ## parse the arguments
  cmake_parse_arguments(arg "${options}" "${singleValueArgs}" 
                            "${multiValueArgs}" ${ARGN} )

  ## check arguments
  if ( NOT DEFINED arg_TARGET_SOURCES )
    message( FATAL_ERROR "Must provide target sources" )
  endif()

  ## append the target source to the all project sources
  foreach( src ${arg_TARGET_SOURCES} )
    if(IS_ABSOLUTE ${src})
       list(APPEND "${PROJECT_NAME}_ALL_SOURCES" "${src}")
     else()
       list(APPEND "${PROJECT_NAME}_ALL_SOURCES"
                   "${CMAKE_CURRENT_SOURCE_DIR}/${src}")
     endif()
  endforeach()

  set( "${PROJECT_NAME}_ALL_SOURCES" "${${PROJECT_NAME}_ALL_SOURCES}"
       CACHE STRING "" FORCE )

endmacro(update_project_sources)
