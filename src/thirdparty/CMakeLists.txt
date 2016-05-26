#
# Our third-party builtin libs are all currently header only,
# so here we are adding tests that make sure they link and can 
# be used
#

######################################################################
# Register the stla_io header-only library for reading ascii STL files
######################################################################
set(STLAIO_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/stla_io )
blt_register_library( 
    NAME stla_io
    INCLUDES ${STLAIO_INCLUDE_DIRS}
    )

############################################################
# Register the fmt header-only library for string formatting 
############################################################
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


