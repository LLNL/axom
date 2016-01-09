#
# Setup targets to generate code.
#
# Each package can create their own ${PROJECT}_generate target
#  add_dependencies(generate  ${PROJECT}_generate)

##------------------------------------------------------------------------------
## add_shroud( INPUT file
##             TIMESTAMP file
##             DEPENDS_SOURCE file1 ... filen
##             DEPENDS_BINARY file1 ... filen
## )
##
##  INPUT          - yaml input file to shroud. Required.
##  TIMESTAMP      - output file touched to mark when shroud was last run.
##                   may be used as a dependency on other targets
##  DEPENDS_SOURCE - splicer files in the source directory
##  DEPENDS_BINARY - splicer files in the binary directory
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
    set(singleValueArgs INPUT TIMESTAMP)
    set(multiValueArgs DEPENDS_SOURCE DEPENDS_BINARY )

    ## parse the arguments to the macro
    cmake_parse_arguments(arg
         "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

    # make sure INPUT is defined
    if(NOT arg_INPUT)
      message(FATAL_ERROR "add_shroud macro must define INPUT")
    endif()

    # default name of TIMESTAMP
    if(NOT arg_TIMESTAMP)
      set(arg_TIMESTAMP ${arg_INPUT}.time)
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
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${arg_INPUT} ${shroud_depends}
        COMMAND ${EXECUTABLE_OUTPUT_PATH}/shroud
                --logdir ${CMAKE_CURRENT_BINARY_DIR}
                # path controls where to search for splicer files listed in INPUT
                --path ${CMAKE_CURRENT_BINARY_DIR}
                --path ${CMAKE_CURRENT_SOURCE_DIR}
                ${CMAKE_CURRENT_SOURCE_DIR}/${arg_INPUT}
        COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/${arg_TIMESTAMP}
        COMMENT "Running shroud ${arg_INPUT}"
        WORKING_DIRECTORY ${SHROUD_OUTPUT_DIR}
    )
endmacro(add_shroud)



add_custom_target(generate)
