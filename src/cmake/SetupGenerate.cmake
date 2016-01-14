#
# Setup targets to generate code.
#
# Each package can create their own ${PROJECT}_generate target
#  add_dependencies(generate  ${PROJECT}_generate)

##------------------------------------------------------------------------------
## add_shroud( YAML_INPUT_FILE file
##             TIMESTAMP file
##             DEPENDS_SOURCE file1 ... filen
##             DEPENDS_BINARY file1 ... filen
## )
##
##  YAML_INPUT_FILE - yaml input file to shroud. Required.
##  TIMESTAMP       - output file touched to mark when shroud was last run.
##                    may be used as a dependency on other targets
##  DEPENDS_SOURCE  - splicer files in the source directory
##  DEPENDS_BINARY  - splicer files in the binary directory
##
## Add a shroud target to generate wrappers.
##
##------------------------------------------------------------------------------

macro(add_shroud)

    # Decide where the output files should be written.
    # For now all files are written into the source directory.
    # This allows them to be source controlled and does not require a library user
    # to generate them.  All they have to do is compile them.
    #set(SHROUD_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR})
    set(SHROUD_OUTPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

    set(options)
    set(singleValueArgs YAML_INPUT_FILE TIMESTAMP)
    set(multiValueArgs DEPENDS_SOURCE DEPENDS_BINARY )

    ## parse the arguments to the macro
    cmake_parse_arguments(arg
         "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

    # make sure YAML_INPUT_FILE is defined
    if(NOT arg_YAML_INPUT_FILE)
      message(FATAL_ERROR "add_shroud macro must define YAML_INPUT_FILE")
    endif()

    # default name of TIMESTAMP
    if(NOT arg_TIMESTAMP)
      set(arg_TIMESTAMP ${arg_YAML_INPUT_FILE}.time)
    endif()

    # convert DEPENDS to full paths
    set(shroud_depends)
    foreach (_file ${arg_DEPENDS_SOURCE})
        list(APPEND shroud_depends "${CMAKE_CURRENT_SOURCE_DIR}/${_file}")
    endforeach ()
    foreach (_file ${arg_DEPENDS_BINARY})
        list(APPEND shroud_depends "${CMAKE_CURRENT_BINARY_DIR}/${_file}")
    endforeach ()

    add_custom_command(
        OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/${arg_TIMESTAMP}
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${arg_YAML_INPUT_FILE} ${shroud_depends}
        COMMAND ${EXECUTABLE_OUTPUT_PATH}/shroud
                --logdir ${CMAKE_CURRENT_BINARY_DIR}
                # path controls where to search for splicer files listed in YAML_INPUT_FILE
                --path ${CMAKE_CURRENT_BINARY_DIR}
                --path ${CMAKE_CURRENT_SOURCE_DIR}
                ${CMAKE_CURRENT_SOURCE_DIR}/${arg_YAML_INPUT_FILE}
        COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/${arg_TIMESTAMP}
        COMMENT "Running shroud ${arg_YAML_INPUT_FILE}"
        WORKING_DIRECTORY ${SHROUD_OUTPUT_DIR}
    )

    # Create target to process this Shroud file
    set(_shroud_target shroud_${arg_YAML_INPUT_FILE})
    add_custom_target(${_shroud_target}
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${arg_TIMESTAMP}
    )

#    if(UNCRUSTIFY_FOUND)
#        add_custom_target("uncrustify_inplace_${PROJECT_NAME}"
#	    DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/${arg_TIMESTAMP}
#            ${UNCRUSTIFY_EXECUTABLE}
#            -c ${CMAKE_CURRENT_SOURCE_DIR}/${arg_CFG_FILE} --no-backup ${arg_SRC_FILES}
#             COMMENT "Running uncrustify to apply code formatting settings.")
#    endif(UNCRUSTIFY_FOUND)

    add_dependencies(generate ${_shroud_target})
endmacro(add_shroud)



add_custom_target(generate)
