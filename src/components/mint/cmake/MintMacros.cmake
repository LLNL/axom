#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------

##------------------------------------------------------------------------------
## mint_fe_basis( NAME [name] HEADER_FILES [file1 [file2 ...]] )
##
## Creates a custom target that handles copying the headers of a finite element
## basis to an appropriate folder in the build-tree and install-tree.
##
## NOTE: the custom target is appeneded to the "mint_depends" list, so, this
## macro assumes that "mint_depends" is already defined.
##------------------------------------------------------------------------------
macro(mint_fe_basis)

  set(options)
  set(singleValueArgs NAME)
  set(multiValueArgs HEADER_FILES )

  # parse macro arguments
  cmake_parse_arguments( arg
       "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN} )

  ## NOTE: the mine_fe_basis macro, assumes that the
  if ( NOT DEFINED mint_depends )
    message( FATAL_ERROR "macro assumes that mint_depends is defined!" )
  endif()

  # check arguments
  if ( NOT DEFINED arg_NAME  OR "${arg_NAME}" STREQUAL "" )
     message( FATAL_ERROR "Must specify a name for the basis" )
  endif()

  # project paths
  set(blt_dir ${CMAKE_SOURCE_DIR}/blt/cmake)
  set(src_path_root ${PROJECT_SOURCE_DIR}/src/fem/shape_functions)
  set(build_path_root ${HEADER_INCLUDES_DIRECTORY}/mint/fem/shape_functions)
  set(install_path_root ${CMAKE_INSTALL_PREFIX}/include/${build_path_root})

  # paths for this basis
  string(TOLOWER ${arg_NAME} basis_name)
  set(basis_src_path ${src_path_root}/${basis_name})
  set(basis_build_path ${build_path_root}/${basis_name})
  set(basis_install_path ${install_path_root}/${basis_name})

  if ( NOT EXISTS ${basis_src_path} )
    message( FATAL_ERROR "Directory: ${basis_src_path} does not exist!" )
  endif()

  # setup the list of basis files by pre-pending the basis_src_path to them
  set(basis_files)
  foreach (hdr_file ${arg_HEADER_FILES})

    if ( NOT EXISTS ${basis_src_path}/${hdr_file} )
      message( FATAL_ERROR "${hdr_file} does not exist!" )
    endif()

    list(APPEND basis_files "${basis_src_path}/${hdr_file}" )
  endforeach()

  # if the list of basis files is not empty, setup custom target and
  # update "mint_depends" accordingly
  list(LENGTH basis_files _length)
  if ( ${_length} GREATER 0 )

    set(target_name "copy_mint_${basis_name}_basis")

    add_custom_target( ${target_name}
       COMMAND ${CMAKE_COMMAND}
              -DOUTPUT_DIRECTORY=${basis_build_path}
              -DLIBHEADERS="${basis_files}"
              -P ${blt_dir}/copy_headers.cmake
       DEPENDS
            ${basis_files}

       WORKING_DIRECTORY
            ${PROJECT_SOURCE_DIR}/src

       COMMENT
         "copy mint ${basis_name} basis header files"
       )

    ## append custom target to mint_depends
    list(APPEND mint_depends ${target_name})

    ## add rules for the make install target
    install( FILES ${basis_files} DESTINATION ${basis_install_path} )

  endif()

endmacro(mint_fe_basis)