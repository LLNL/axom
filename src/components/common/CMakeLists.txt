################################
# Datastore
################################
project(common)


################################
# Add the common sources
################################
add_subdirectory(src)

################################
# Add tests
################################
if (ENABLE_TESTS)
  add_subdirectory(tests)
endif()


add_code_check_targets(uncrustify.cfg)

################################
# Add docs
################################
if (ENABLE_DOCS)
  add_subdirectory(docs)
endif()

################################
# Create CMake importable
# exports for all of our targets
################################
install(EXPORT ${PROJECT_NAME}-targets DESTINATION lib/cmake)
