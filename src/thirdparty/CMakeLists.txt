#
# Our third-party builtin libs are all currently header only,
# so here we are adding tests that make sure they link and can 
# be used
#

############################################################
# Register the fmt header-only library for string formatting 
############################################################

set(fmtHeaders fmt/fmt.hpp fmt/fmt/format.h fmt/fmt/format.cc)
set(fmtDest ${HEADER_INCLUDES_DIRECTORY}/fmt)
blt_copy_headers_target(NAME        fmt
                        HEADERS     ${fmtHeaders}
                        DESTINATION ${fmtDest})

blt_register_library( 
    NAME fmt
    INCLUDES ${HEADER_INCLUDES_DIRECTORY}
    )


################################
# Add tests
################################
if (ENABLE_TESTS)
  add_subdirectory(tests)
endif()


