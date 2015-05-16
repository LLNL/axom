#
#
#

#
# Write the driver script
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/shroud.in
               ${CMAKE_CURRENT_SOURCE_DIR}/shroud
#               ${EXECUTABLE_OUTPUT_PATH}/shroud
               @ONLY)
