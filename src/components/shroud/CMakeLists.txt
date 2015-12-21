#
#
#
project(shroud)

add_subdirectory(docs)

if(ENABLE_FORTRAN)
    # Build runtime library
    add_subdirectory(runtime)
endif()

# Make sure Python and required modules are available
if(PYTHON_EXECUTABLE)
    execute_process(
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/test_python.py
        RESULT_VARIABLE HAVE_PYTHON_MODULES
        OUTPUT_VARIABLE MISSING_MODULES
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(HAVE_PYTHON_MODULES EQUAL 0)
        message(STATUS "Able to run shroud")
        set(ENABLE_SHROUD TRUE CACHE BOOLEAN "" FORCE )
    
    else()
        message(STATUS "Unable to run shroud: ${MISSING_MODULES}")
        set(ENABLE_SHROUD FALSE CACHE BOOLEAN "" FORCE )
    endif()

#####
# get info about Python
#####
    execute_process(
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/python_location.py
        OUTPUT_VARIABLE PYCONFIG
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    message(STATUS "Python Config ${PYCONFIG}")
    message(STATUS "Python Config")

    list(GET PYCONFIG 0  PYTHON_MAJOR_VERSION)
    list(GET PYCONFIG 1  PYTHON_MINOR_VERSION)
    list(GET PYCONFIG 2  PYTHON_MICRO_VERSION)
    list(GET PYCONFIG 3  PYTHON_PREFIX)
    message(STATUS "Python Version ${PYTHON_MAJOR_VERSION}.${PYTHON_MINOR_VERSION}.${PYTHON_MICRO_VERSION}")
    message(STATUS "Python Prefix ${PYTHON_PREFIX}")

    set(PYTHON_INCLUDE ${PYTHON_PREFIX}/include/python${PYTHON_MAJOR_VERSION}.${PYTHON_MINOR_VERSION})


else()
    message(STATUS "Unable to run shroud: Python not found")
    set(ENABLE_SHROUD FALSE CACHE BOOLEAN "" FORCE )
endif()

if(ENABLE_FORTRAN AND ENABLE_SHROUD)
    add_subdirectory(tests/run-tutorial)

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
endif()


################################
# Create CMake importable
# exports for all of our targets
################################
install(EXPORT ${PROJECT_NAME}-targets DESTINATION lib/cmake)
