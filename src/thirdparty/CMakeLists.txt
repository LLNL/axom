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

# We have to divide into two "libraries" to get the copied headers in
# the proper directory configuration, some in sparsehash and some in
# sparsehash/internal

set(sparsehash_headers
    sparsehash/dense_hash_map
    sparsehash/dense_hash_set
    sparsehash/sparse_hash_map
    sparsehash/sparse_hash_set
    sparsehash/sparsetable
    sparsehash/template_util.h
    sparsehash/type_traits.h
    )

set(sparsehash_internal_headers
    sparsehash/internal/densehashtable.h
    sparsehash/internal/hashtable-common.h
    sparsehash/internal/libc_allocator_with_realloc.h
    sparsehash/internal/sparsehashtable.h
    sparsehash/internal/sparseconfig.h
    )

set(sparsehash_dest ${HEADER_INCLUDES_DIRECTORY}/sparsehash)
blt_copy_headers_target(NAME        sparsehash
                        HEADERS     ${sparsehash_headers}
                        DESTINATION ${sparsehash_dest})

set(sparsehash_internal_dest ${HEADER_INCLUDES_DIRECTORY}/sparsehash/internal)
blt_copy_headers_target(NAME        sparsehash_internal
                        HEADERS     ${sparsehash_internal_headers}
                        DESTINATION ${sparsehash_internal_dest})

blt_register_library( 
    NAME sparsehash
    INCLUDES ${HEADER_INCLUDES_DIRECTORY}
    )

blt_register_library( 
    NAME sparsehash_internal
    )

set(SPARSEHASH_FOUND TRUE CACHE INTERNAL "")
mark_as_advanced(SPARSEHASH_FOUND)


################################
# Add tests
################################
if (ENABLE_TESTS)
  add_subdirectory(tests)
endif()


