#
#
#
if (ENABLE_FORTRAN)
    add_subdirectory(runtime)
endif()

#
# Write the driver script
# configure into CMakeFiles directory, because it must be in a different directory
# since file(COPY) does rename files
set(tmp_file ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/shroud)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/shroud.in
               ${tmp_file}
               @ONLY)
# now copy the temporary into the final destination, setting the permissions
file(COPY ${tmp_file}
    DESTINATION ${EXECUTABLE_OUTPUT_PATH}
    FILE_PERMISSIONS
        OWNER_READ OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE)
file(REMOVE ${tmp_file})

#
# Write the test script
# configure into CMakeFiles directory, because it must be in a different directory
# since file(COPY) does rename files
set(tmp_file ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/do-test)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tests/do-test.in
               ${tmp_file}
               @ONLY)
# now copy the temporary into the final destination, setting the permissions
file(COPY ${tmp_file}
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
    FILE_PERMISSIONS
        OWNER_READ OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE)
file(REMOVE ${tmp_file})



add_test(NAME shroud
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/do-test
)
