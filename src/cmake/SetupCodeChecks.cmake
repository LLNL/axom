###############################################################################
# Targets related to source code checks (formatting, static analysis, etc)
###############################################################################

add_custom_target(check)

if(UNCRUSTIFY_FOUND)
    # targets for verifying formatting
    add_custom_target(uncrustify_check)
    add_dependencies(check uncrustify_check)
    # targets for forcing formatting
    add_custom_target(uncrustify_inplace)

endif()

##------------------------------------------------------------------------------
## - Macro that helps add all code check targets
##
## add_code_check_targets(cfg)
##
##------------------------------------------------------------------------------
macro(add_code_check_targets cfg_file)

    # Only run uncrustify on C and C++ files
    # Note, we can later extend this by passing in a list of valid types to the macro
    set(_fileTypes ".cpp" ".hpp" ".c" ".h")
    
    # generate the filtered list of source files
    set(_filt_sources)
    foreach(_file ${${PROJECT_NAME}_ALL_SOURCES})
      get_filename_component(_ext ${_file} EXT)
      list(FIND _fileTypes "${_ext}" _index)
      
      if(_index GREATER -1)
         list(APPEND _filt_sources ${_file})
      endif()
    endforeach()

    if(UNCRUSTIFY_FOUND)
        add_uncrustify_check(CFG_FILE ${cfg_file}   SRC_FILES ${_filt_sources})
        add_uncrustify_inplace(CFG_FILE ${cfg_file} SRC_FILES ${_filt_sources})
    endif()
endmacro()
    

##------------------------------------------------------------------------------
## - Macro for invoking uncrustify to check code formatting
##
## add_uncrustify_check( CFG_FILE <uncrusify_configuration_file> 
##                       SRC_FILES <list_of_src_files_to_uncrustify> )
##
##------------------------------------------------------------------------------
macro(add_uncrustify_check)
    
    MESSAGE(STATUS "Creating uncrustify check target: uncrustify_check_${PROJECT_NAME}")

    ## parse the arguments to the macro
    set(options)
    set(singleValueArgs CFG_FILE)
    set(multiValueArgs SRC_FILES)
    cmake_parse_arguments(arg
        "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN} )

    add_custom_target("uncrustify_check_${PROJECT_NAME}"
            ${UNCRUSTIFY_EXECUTABLE}
            -c ${CMAKE_CURRENT_SOURCE_DIR}/${arg_CFG_FILE} --check ${arg_SRC_FILES}
             COMMENT "Running uncrustify source code formatting checks.")
        
    # hook our new target into the check dependency chain
    add_dependencies(uncrustify_check "uncrustify_check_${PROJECT_NAME}")

endmacro(add_uncrustify_check)

##------------------------------------------------------------------------------
## - Macro for invoking uncrustify to apply formatting inplace
##
## add_uncrustify_inplace(CFG_FILE <uncrusify_configuration_file> 
##                        SRC_FILES <list_of_src_files_to_uncrustify> )
##
##------------------------------------------------------------------------------
macro(add_uncrustify_inplace)
    
    MESSAGE(STATUS "Creating uncrustify inplace target: uncrustify_inplace_${PROJECT_NAME}")

    ## parse the arguments to the macro
    set(options)
    set(singleValueArgs CFG_FILE)
    set(multiValueArgs SRC_FILES)
    cmake_parse_arguments(arg
        "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN} )

    add_custom_target("uncrustify_inplace_${PROJECT_NAME}"
            ${UNCRUSTIFY_EXECUTABLE}
            -c ${CMAKE_CURRENT_SOURCE_DIR}/${arg_CFG_FILE} --no-backup ${arg_SRC_FILES}
             COMMENT "Running uncrustify to apply code formatting settings.")
        
    # hook our new target into the uncrustify_inplace dependency chain
    add_dependencies(uncrustify_inplace "uncrustify_inplace_${PROJECT_NAME}")

endmacro(add_uncrustify_inplace)
