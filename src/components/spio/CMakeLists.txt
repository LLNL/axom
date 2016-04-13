################################
# IOParallel
################################
project(spio)

################################
# Add the sources
################################
add_subdirectory(src)

################################
# Add examples
################################
if (ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()


################################
# Add tests
################################
if (ENABLE_TESTS)
  add_subdirectory(tests)
endif()


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
