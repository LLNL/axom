set(API_OUTPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
#
# Generate source
#
add_custom_command(
    OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/slic_shroud_api.time
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/slic_api.yaml
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/shroud --logdir ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/slic_api.yaml
    COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/slic_shroud_api.time
    COMMENT "Generating shroud api for SLIC"
    WORKING_DIRECTORY ${API_OUTPUT_DIR}
)

add_custom_target(slic_generate
     DEPENDS
         ${CMAKE_CURRENT_BINARY_DIR}/slic_shroud_api.time
)

add_dependencies(generate slic_generate)