################################
# Datastore
################################
project(datastore)

################################
# Add the datastore sources
################################
add_subdirectory(src)

################################
# Add examples
################################
add_subdirectory(examples)

################################
# Add tests
################################
if (BUILD_TESTING)
  add_subdirectory(tests)
endif()


add_uncrustify_check(uncrustify.cfg)

################################
# Create CMake importable
# exports for all of our targets
################################
install(EXPORT ${PROJECT_NAME}-targets DESTINATION lib/cmake)
