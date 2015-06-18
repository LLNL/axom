#
# sidre/capi
# Specify all headers
#
set(sidre_headers
    sidre.h
    SidreTypes.h
    SidreWrapperHelpers.hpp
)

set(sidre_generated_headers
    wrapSiDRe.h
    wrapDataStore.h
    wrapDataGroup.h
    wrapDataBuffer.h
    wrapDataView.h
)

#
# Specify all sources
#
set(sidre_generated_source
    wrapSiDRe.cpp
    wrapDataStore.cpp
    wrapDataGroup.cpp
    wrapDataBuffer.cpp
    wrapDataView.cpp
    fstr_mod.f
    wrapfsidre.f
    )

#set(API_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(API_OUTPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
#
# Generate source
#
#add_custom_target(sidre_generate
#    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/capi.time
#)

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/capi.time
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/api.yaml
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/shroud --logdir ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/api.yaml
    COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/capi.time
    COMMENT "Running shroud api.yaml"
    WORKING_DIRECTORY ${API_OUTPUT_DIR}
)

#
# Tell parent about our files.
#

set(generated_list)
set(parent_list)

# static files
foreach(file ${sidre_headers})
  list(APPEND parent_list ${CMAKE_CURRENT_SOURCE_DIR}/${file} )

#  add_custom_command(
#      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${file}
#  )
endforeach(file)

# generated headers
foreach(file ${sidre_generated_headers})
  list(APPEND parent_list ${API_OUTPUT_DIR}/${file} )

  add_custom_command(
      OUTPUT ${API_OUTPUT_DIR}/${file}
      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/capi.time
  )
#  SET_SOURCE_FILES_PROPERTIES(${API_OUTPUT_DIR}/${file}
#       PROPERTIES GENERATED TRUE
#  ) 
endforeach(file)
set (capi_headers ${parent_list} PARENT_SCOPE)
list(APPEND generated_list ${parent_list})


set(parent_list)
# generated source
foreach(file ${sidre_source} ${sidre_generated_source})
  list(APPEND parent_list ${API_OUTPUT_DIR}/${file} )

  add_custom_command(
      OUTPUT ${API_OUTPUT_DIR}/${file}
      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/capi.time
  )
#  SET_SOURCE_FILES_PROPERTIES(${API_OUTPUT_DIR}/${file}
#       PROPERTIES GENERATED TRUE
#  )

endforeach(file)
set (capi_source ${parent_list} PARENT_SCOPE)
list(APPEND generated_list ${parent_list})


ADD_CUSTOM_TARGET(sidre_generate
     DEPENDS ${generated_list}
)

#
# Mark the wappers as done by touching capi.time
# This should avoid generating the wrappers unless api.yaml is changed.
#
execute_process(
    COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/capi.time
#--    COMMAND ${EXECUTABLE_OUTPUT_PATH}/shroud
#--         --logdir ${CMAKE_CURRENT_BINARY_DIR}
#--         --outdir ${API_OUTPUT_DIR}
#--         ${CMAKE_CURRENT_SOURCE_DIR}/api.yaml
#--    RESULT_VARIABLE capi_result
#--    ERROR_VARIABLE capi_stderr
#--    OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/capi.output
    OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/capi.time
)
#--if(capi_result)
#--  message(FATAL_ERROR "Error creating sidre wrappers:\n${capi_stderr}")
#--else()
#--  message(STATUS "Created sidre wrappers")
#--endif()
