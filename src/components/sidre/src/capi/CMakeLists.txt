#
# sidre/capi
# Specify all headers
#
set(sidre_headers
    sidre.h
    SidreTypes.h
)

set(sidre_generated_headers
    wrapDataStore.h
    wrapDataGroup.h
    wrapDataBuffer.h
    wrapDataView.h
)

#
# Specify all sources
#
set(sidre_generated_source
    wrapDataStore.cpp
    wrapDataGroup.cpp
    wrapDataBuffer.cpp
    wrapDataView.cpp
#    fstr_mod.f
#    wrapfDataStore.f
    )

#set(API_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(API_OUTPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
#
# Generate source
#
add_custom_target(sidre_generate
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/capi.time
    COMMENT "Running shroud api.yaml"
)

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/capi.time
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/api.yaml
    COMMAND sh ${EXECUTABLE_OUTPUT_PATH}/shroud ${CMAKE_CURRENT_SOURCE_DIR}/api.yaml
    COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/capi.time
    WORKING_DIRECTORY ${API_OUTPUT_DIR}
)

#
# Tell parent about our files.
#
set(parent_list)
foreach(file ${sidre_headers})
  list(APPEND parent_list capi/${file} )
endforeach(file)
foreach(file ${sidre_generated_headers})
  list(APPEND parent_list ${API_OUTPUT_DIR}/${file} )
endforeach(file)
set (capi_headers ${parent_list} PARENT_SCOPE)

set(parent_list)
foreach(file ${sidre_source} ${sidre_generated_source})
  list(APPEND parent_list ${API_OUTPUT_DIR}/${file} )
endforeach(file)
set (capi_source ${parent_list} PARENT_SCOPE)
