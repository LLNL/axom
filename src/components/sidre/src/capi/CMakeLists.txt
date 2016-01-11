#
# Rules to generate source
#
# A file is generated in the binary directory when a generate command has been run.
# Required for capi.time since it produces several output files.

#
# Generate source
#
add_shroud( YAML_INPUT_FILE api.yaml
            TIMESTAMP capi.time
            DEPENDS_SOURCE csidresplicer.c fsidresplicer.f pysidresplicer.c
            DEPENDS_BINARY genfsidresplicer.f
)

add_custom_command(
    OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/genfsidresplicer.f
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/genfsidresplicer.py
    COMMAND ${PYTHON_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/genfsidresplicer.py fortran > 
            ${CMAKE_CURRENT_BINARY_DIR}/genfsidresplicer.f
)

ADD_CUSTOM_TARGET(sidre_generate
     DEPENDS
         ${CMAKE_CURRENT_BINARY_DIR}/capi.time
)

add_dependencies(generate sidre_generate)

#
# Mark the wappers as done by touching capi.time
# This should avoid generating the wrappers unless api.yaml is changed.
#
execute_process(
    COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/genfsidresplicer.time
    OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/genfsidresplicer.time
)
execute_process(
    COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/capi.time
    OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/capi.time
)
