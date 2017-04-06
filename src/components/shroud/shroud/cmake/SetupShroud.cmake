# Setup Shroud
# This file defines:
#  SHROUD_FOUND - If Shroud was found

if(NOT SHROUD_EXECUTABLE)
    MESSAGE(FATAL_ERROR "Could not find Shroud. Shroud requires explicit SHROUD_EXECUTABLE.")
endif()

message(STATUS Found SHROUD: ${SHROUD_EXECUTABLE})

add_custom_target(generate)
set(SHROUD_FOUND TRUE)

#
# Setup targets to generate code.
#
# Each package can create their own ${PROJECT}_generate target
#  add_dependencies(generate  ${PROJECT}_generate)

##------------------------------------------------------------------------------
## add_shroud( YAML_INPUT_FILE file
##             DEPENDS_SOURCE file1 ... filen
##             DEPENDS_BINARY file1 ... filen
##             C_FORTRAN_OUTPUT_DIR dir
##             PYTHON_OUTPUT_DIR dir
## )
##
##  YAML_INPUT_FILE - yaml input file to shroud. Required.
##  DEPENDS_SOURCE  - splicer files in the source directory
##  DEPENDS_BINARY  - splicer files in the binary directory
##  C_FORTRAN_OUTPUT_DIR - directory for C and Fortran wrapper output files.
##  PYTHON_OUTPUT_DIR - directory for Python wrapper output files.
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
    set(singleValueArgs
        YAML_INPUT_FILE
        C_FORTRAN_OUTPUT_DIR
        PYTHON_OUTPUT_DIR
    )
    set(multiValueArgs DEPENDS_SOURCE DEPENDS_BINARY )

    ## parse the arguments to the macro
    cmake_parse_arguments(arg
         "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

    # make sure YAML_INPUT_FILE is defined
    if(NOT arg_YAML_INPUT_FILE)
      message(FATAL_ERROR "add_shroud macro must define YAML_INPUT_FILE")
    endif()

    if(arg_C_FORTRAN_OUTPUT_DIR)
      set(SHROUD_C_FORTRAN_OUTPUT_DIR --outdir-c-fortran ${arg_C_FORTRAN_OUTPUT_DIR})
    endif()

    if(arg_PYTHON_OUTPUT_DIR)
      set(SHROUD_PYTHON_OUTPUT_DIR --outdir-python ${arg_PYTHON_OUTPUT_DIR})
    endif()

    # convert DEPENDS to full paths
    set(shroud_depends)
    foreach (_file ${arg_DEPENDS_SOURCE})
        list(APPEND shroud_depends "${CMAKE_CURRENT_SOURCE_DIR}/${_file}")
    endforeach ()
    foreach (_file ${arg_DEPENDS_BINARY})
        list(APPEND shroud_depends "${CMAKE_CURRENT_BINARY_DIR}/${_file}")
    endforeach ()

    get_filename_component(_basename ${arg_YAML_INPUT_FILE} NAME_WE)
    set(_timestamp  ${CMAKE_CURRENT_BINARY_DIR}/${_basename}.time)
    set(_uncrustify ${CMAKE_CURRENT_BINARY_DIR}/${_basename}.uncrustify)
    set(_cfiles     ${CMAKE_CURRENT_BINARY_DIR}/${_basename}.cfiles)
    set(_ffiles     ${CMAKE_CURRENT_BINARY_DIR}/${_basename}.ffiles)

    set(_cmd
        ${SHROUD_EXECUTABLE}
        --logdir ${CMAKE_CURRENT_BINARY_DIR}
        ${SHROUD_C_FORTRAN_OUTPUT_DIR}
        ${SHROUD_PYTHON_OUTPUT_DIR}
        # path controls where to search for splicer files listed in YAML_INPUT_FILE
        --path ${CMAKE_CURRENT_BINARY_DIR}
        --path ${CMAKE_CURRENT_SOURCE_DIR}
        --cfiles ${_cfiles}
        --ffiles ${_ffiles}
        ${CMAKE_CURRENT_SOURCE_DIR}/${arg_YAML_INPUT_FILE}
    )

    add_custom_command(
        OUTPUT  ${_timestamp}
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${arg_YAML_INPUT_FILE} ${shroud_depends}
        COMMAND ${_cmd}
        COMMAND touch ${_timestamp}
        COMMAND rm -f ${_uncrustify}
        COMMENT "Running shroud ${arg_YAML_INPUT_FILE}"
        WORKING_DIRECTORY ${SHROUD_OUTPUT_DIR}
    )

    # Create target to process this Shroud file
    set(_shroud_target generate_${_basename})
    add_custom_target(${_shroud_target}
        DEPENDS ${_timestamp}
    )

    # Only run uncrustify if shroud has just run
    # XXX uncrustify.cfg is hardwired
    set(_cfg ${PROJECT_SOURCE_DIR}/uncrustify.cfg)
    if(UNCRUSTIFY_FOUND AND (EXISTS ${_cfg}))
        add_custom_command(
            OUTPUT ${_uncrustify}
            DEPENDS  ${_timestamp}
            COMMAND ${UNCRUSTIFY_EXECUTABLE}
                    -c ${_cfg} --no-backup `cat ${_cfiles}`
            COMMAND touch ${_uncrustify}
            COMMENT "Running uncrustify for ${arg_YAML_INPUT_FILE}."
            WORKING_DIRECTORY ${SHROUD_OUTPUT_DIR}
        )
        add_custom_target(${_shroud_target}_uncrustify
            DEPENDS ${_uncrustify}
        )
        add_dependencies(generate ${_shroud_target}_uncrustify)
    endif()

    add_dependencies(generate ${_shroud_target})
endmacro(add_shroud)
