#
# Rules to generate source
#
# A file is generated in the binary directory when a generate command has been run.
# Required for capi.time since it produces several output files.
# Required for genfsidresplicer.time since CMake wants to delete any generated files during 'make clean'

#set(API_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(API_OUTPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
#
# Generate source
#
add_custom_command(
    OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/capi.time
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/genfsidresplicer.time
            ${CMAKE_CURRENT_SOURCE_DIR}/api.yaml
            ${CMAKE_CURRENT_SOURCE_DIR}/csidresplicer.c
            ${CMAKE_CURRENT_SOURCE_DIR}/fsidresplicer.f
            ${CMAKE_CURRENT_SOURCE_DIR}/genfsidresplicer.f
            ${CMAKE_CURRENT_SOURCE_DIR}/pysidresplicer.c
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/shroud --logdir ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/api.yaml
    COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/capi.time
    COMMENT "Running shroud api.yaml"
    WORKING_DIRECTORY ${API_OUTPUT_DIR}
)

add_custom_command(
    OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/genfsidresplicer.time
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/genfsidresplicer.py
    COMMAND ${PYTHON_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/genfsidresplicer.py > 
            ${CMAKE_CURRENT_SOURCE_DIR}/genfsidresplicer.f
    COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/genfsidresplicer.time
    COMMENT "Running genfsidresplicer.py"
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
