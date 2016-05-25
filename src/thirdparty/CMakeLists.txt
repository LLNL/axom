#
# Our third-party builtin libs are all currently header only,
# so here we are adding tests that make sure they link and can 
# be used
#

set(FMT_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/fmt )
blt_register_library( 
    NAME fmt
    INCLUDES ${FMT_INCLUDE_DIRS}
    )


################################
# Add tests
################################
if (ENABLE_TESTS)
  add_subdirectory(tests)
endif()


