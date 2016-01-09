set(API_OUTPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
#
# Generate source
#
add_shroud( INPUT slic_api.yaml )

add_custom_target(slic_generate
     DEPENDS
         ${CMAKE_CURRENT_BINARY_DIR}/slic_api.yaml.time
)

add_dependencies(generate slic_generate)
