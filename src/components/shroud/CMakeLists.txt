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

####
    # Run the develop command.
    # The Python files are not installed directly, only an egg.
    # So any changes to the source does not require a new install.
    execute_process(
        COMMAND ${PYTHON_EXECUTABLE} setup.py
          develop
          --install-dir ${CMAKE_Python_MODULE_DIRECTORY}
          --script-dir  ${EXECUTABLE_OUTPUT_PATH}
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
         OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/setup.log
    )

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


################################

#install(CODE "execute_process(COMMAND python setup.py install -f
#--prefix=${CMAKE_INSTALL_PREFIX} WORKING_DIRECTORY
#../source/libmoleculizer-1.1.2/python-src/language_parser)")

#https://mail.python.org/pipermail/distutils-sig/2013-August/022184.html

#https://bloerg.net/2012/11/10/cmake-and-distutils.html
#if (PYTHON)
#    set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
#    set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
#    set(DEPS        "${CMAKE_CURRENT_SOURCE_DIR}/module/__init__.py")
#    set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/build/timestamp")
#
#    configure_file(${SETUP_PY_IN} ${SETUP_PY})
#
#    add_custom_command(OUTPUT ${OUTPUT}
#                       COMMAND ${PYTHON} ${SETUP_PY} build
#                       COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
#                       DEPENDS ${DEPS})
#
#    add_custom_target(target ALL DEPENDS ${OUTPUT})
#
#    install(CODE "execute_process(COMMAND ${PYTHON} ${SETUP_PY} install)")
#endif()
