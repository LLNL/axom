###############################################################################
# Targets related to source code checks (formatting, static analysis, etc)
###############################################################################

add_custom_target(check)

if(UNCRUSTIFY_FOUND)
    add_custom_target(uncrustify_check)
    add_dependencies(check uncrustify_check)
endif()

##------------------------------------------------------------------------------
## - Macro for invoking uncrustify to check code formatting
##
## add_uncrustify_check(project_name cfg source_files)
##
##------------------------------------------------------------------------------
macro(add_uncrustify_check cfg_file)
    
    MESSAGE(STATUS "Creating uncrustify check target ${PROJECT_NAME}_uncrustify_check")

    add_custom_target("${PROJECT_NAME}_uncrustify_check"
            ${UNCRUSTIFY_EXECUTABLE}
            -c ${CMAKE_CURRENT_SOURCE_DIR}/${cfg_file} --check ${${PROJECT_NAME}_ALL_SOURCES}
             COMMENT "Running uncrustify source code formatting checks.")
        
    # hook our new target into the check dependency chain
    add_dependencies(uncrustify_check "${PROJECT_NAME}_uncrustify_check")

endmacro(add_uncrustify_check)





