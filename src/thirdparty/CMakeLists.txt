#
# Our third-party builtin libs are all currently header only,
# so here we are adding tests that make sure they link and can 
# be used
#

################################
# Add tests
################################
if (BUILD_TESTS)
  add_subdirectory(tests)
endif()

