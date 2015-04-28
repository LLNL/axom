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
    add_uncrustify_check(${cfg_file})
    add_uncrustify_inplace(${cfg_file})
endmacro()
    

##------------------------------------------------------------------------------
## - Macro for invoking uncrustify to check code formatting
##
## add_uncrustify_check(cfg)
##
##------------------------------------------------------------------------------
macro(add_uncrustify_check cfg_file)
    
    MESSAGE(STATUS "Creating uncrustify check target: uncrustify_check_${PROJECT_NAME}")

    add_custom_target("uncrustify_check_${PROJECT_NAME}"
            ${UNCRUSTIFY_EXECUTABLE}
            -c ${CMAKE_CURRENT_SOURCE_DIR}/${cfg_file} --check ${${PROJECT_NAME}_ALL_SOURCES}
             COMMENT "Running uncrustify source code formatting checks.")
        
    # hook our new target into the check dependency chain
    add_dependencies(uncrustify_check "uncrustify_check_${PROJECT_NAME}")

endmacro(add_uncrustify_check)

##------------------------------------------------------------------------------
## - Macro for invoking uncrustify to apply formatting inplace
##
## add_uncrustify_inplace(cfg)
##
##------------------------------------------------------------------------------
macro(add_uncrustify_inplace cfg_file)
    
    MESSAGE(STATUS "Creating uncrustify inplace target: uncrustify_inplace_${PROJECT_NAME}")

    add_custom_target("uncrustify_inplace_${PROJECT_NAME}"
            ${UNCRUSTIFY_EXECUTABLE}
            -c ${CMAKE_CURRENT_SOURCE_DIR}/${cfg_file} --no-backup ${${PROJECT_NAME}_ALL_SOURCES}
             COMMENT "Running uncrustify to apply code formatting settings.")
        
    # hook our new target into the uncrustify_inplace dependency chain
    add_dependencies(uncrustify_inplace "uncrustify_inplace_${PROJECT_NAME}")

endmacro(add_uncrustify_inplace)




