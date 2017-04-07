#
# Our third-party builtin libs are all currently header only,
# so here we are adding tests that make sure they link and can 
# be used
#

############################################################
# Register the fmt header-only library for string formatting 
############################################################

set(fmtHeaders 
    fmt/fmt.hpp 
    fmt/fmt/format.h 
    fmt/fmt/format.cc
    fmt/fmt/ostream.h 
    fmt/fmt/ostream.cc    
    fmt/fmt/printf.h    
    )
set(fmtDest ${HEADER_INCLUDES_DIRECTORY}/fmt)
blt_copy_headers_target(NAME        fmt
                        HEADERS     ${fmtHeaders}
                        DESTINATION ${fmtDest})

blt_register_library( 
    NAME fmt
    INCLUDES ${HEADER_INCLUDES_DIRECTORY}
    )

set(FMT_FOUND TRUE CACHE INTERNAL "")
mark_as_advanced(FMT_FOUND)

############################################################
# Register the sparsehash header-only library for hashed map containers
############################################################

set(sparsehash_includes
    ${PROJECT_SOURCE_DIR}/thirdparty/sparsehash)
 
blt_register_library( 
    NAME sparsehash
    INCLUDES ${sparsehash_includes}
    )

set(SPARSEHASH_FOUND TRUE CACHE INTERNAL "")
mark_as_advanced(SPARSEHASH_FOUND)


################################
# Add tests
################################
if (ENABLE_TESTS)
  add_subdirectory(tests)
endif()


