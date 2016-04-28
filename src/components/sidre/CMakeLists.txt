################################
# Datastore
################################
project(sidre)

################################
# Check necessary dependencies
################################
if(NOT CONDUIT_FOUND)
  message(FATAL_ERROR "Sidre requires Conduit. Set CONDUIT_DIR to location of built Conduit.")
endif()

################################
# Add the sidre sources
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
