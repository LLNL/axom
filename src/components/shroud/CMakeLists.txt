#
#
#
project(shroud)

if (BUILD_DOCS)
  add_subdirectory(docs)
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
 else()
	message(FATAL_ERROR "Unable to run shroud: ${MISSING_MODULES}")
 endif()

####
  PYTHON_ADD_DISTUTILS_SETUP(shroud_python
	 src/__init__.py  src/parse_decl.py  src/util.py   src/wrapf.py
	 src/main.py      src/splicer.py     src/wrapc.py  src/wrapp.py
	 src/whelpers.py
  )
else()
  message(FATAL_ERROR "Unable to run shroud: Python not found")
endif()

# Create setup.py in binary directory
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/setup.py
)
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
	  DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
	  FILE_PERMISSIONS
		 OWNER_READ OWNER_EXECUTE
		 GROUP_READ GROUP_EXECUTE
		 WORLD_READ WORLD_EXECUTE)
file(REMOVE ${tmp_file})

#
# Write the test script
# configure into CMakeFiles directory, because it must be in a different directory
# since file(COPY) does rename files
if(BUILD_TESTING)
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


  add_test(
	 NAME shroud
	 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/do-test
  )
endif()

if(ENABLE_FORTRAN)
  if (BUILD_TESTING)
	 add_subdirectory(tests/run-tutorial)
	 add_subdirectory(tests/run-strings)
  endif()
endif()
