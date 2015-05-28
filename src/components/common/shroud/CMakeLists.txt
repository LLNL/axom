#
#
#

#
# Write the driver script
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/shroud.in
               ${EXECUTABLE_OUTPUT_PATH}/shroud
               @ONLY)
# TODO: shroud also needs to be executable
