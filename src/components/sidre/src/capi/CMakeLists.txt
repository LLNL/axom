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
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/api.yaml
            ${CMAKE_CURRENT_SOURCE_DIR}/csidresplicer.c
            ${CMAKE_CURRENT_SOURCE_DIR}/fsidresplicer.f
            ${CMAKE_CURRENT_SOURCE_DIR}/pysidresplicer.c
	    # generated
            ${CMAKE_CURRENT_BINARY_DIR}/genfsidresplicer.f
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/shroud
            --logdir ${CMAKE_CURRENT_BINARY_DIR}
            --path ${CMAKE_CURRENT_BINARY_DIR}
            --path ${CMAKE_CURRENT_SOURCE_DIR}
            ${CMAKE_CURRENT_SOURCE_DIR}/api.yaml
    COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/capi.time
    COMMENT "Running shroud api.yaml"
    WORKING_DIRECTORY ${API_OUTPUT_DIR}
)

add_custom_command(
    OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/genfsidresplicer.f
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/genfsidresplicer.py
    COMMAND ${PYTHON_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/genfsidresplicer.py fortran > 
            ${CMAKE_CURRENT_BINARY_DIR}/genfsidresplicer.f
)

#
# When using cog, create an intermediate timestamp file since CMake
# will consider it generated and remove it in some situations (make clean)
# Need to keep source around since it is read as input.
#
get_filename_component(PyDir ${PYTHON_EXECUTABLE} DIRECTORY)
set(COG_EXECUTABLE ${PyDir}/cog.py)

add_custom_command(
    OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/SidreAllocatable.cpp.time
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/genfsidresplicer.py
            ${CMAKE_CURRENT_SOURCE_DIR}/SidreAllocatable.cpp
    COMMAND ${COG_EXECUTABLE} -r ${CMAKE_CURRENT_SOURCE_DIR}/SidreAllocatable.cpp
    COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/SidreAllocatable.cpp.time
    WORKING_DIRECTORY ${API_OUTPUT_DIR}
)

add_custom_command(
    OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/SidreAllocatablef.f.time
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/SidreAllocatablef.f
            ${CMAKE_CURRENT_SOURCE_DIR}/genfsidresplicer.py
    COMMAND ${COG_EXECUTABLE} -r ${CMAKE_CURRENT_SOURCE_DIR}/SidreAllocatablef.f
    COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/SidreAllocatablef.f.time
    WORKING_DIRECTORY ${API_OUTPUT_DIR}
)

add_custom_target(sidre_cog
                  DEPENDS
                      ${CMAKE_CURRENT_BINARY_DIR}/SidreAllocatable.cpp.time
                      ${CMAKE_CURRENT_BINARY_DIR}/SidreAllocatablef.f.time
)

ADD_CUSTOM_TARGET(sidre_generate
     DEPENDS
         ${CMAKE_CURRENT_BINARY_DIR}/capi.time
         sidre_cog
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
