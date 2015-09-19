################################
# Datastore
################################
project(common)


################################
# Add the common sources
################################
add_subdirectory(src)
#add_subdirectory(chasm)
add_subdirectory(shroud)

################################
# Add tests
################################
if (BUILD_TESTING)
  add_subdirectory(tests)
endif()


add_code_check_targets(uncrustify.cfg)

################################
# Create CMake importable
# exports for all of our targets
################################
#install(EXPORT ${PROJECT_NAME}-targets DESTINATION lib/cmake)
