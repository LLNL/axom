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
if (BUILD_TESTING)
#  No tests for common component yet.
#  We'll have to uncomment the following line 
#  if/when we add tests to common
#
#  add_subdirectory(tests)
endif()


add_code_check_targets(uncrustify.cfg)

################################
# Add docs
################################
add_subdirectory(docs)



################################
# Create CMake importable
# exports for all of our targets
################################
install(EXPORT ${PROJECT_NAME}-targets DESTINATION lib/cmake)
