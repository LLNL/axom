#------------------------------------------------------------------------------
# Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
#------------------------------------------------------------------------------

#
# Axom third-party builtin libs are all currently header only.
# We have tests here that make sure they link and can be used.
#

############################################################
# Register the fmt header-only library for string formatting
############################################################

set(fmt_headers
    fmt/format.h
    fmt/ostream.h
    fmt/printf.h
    )

# Do not include the following source files in the headers list even
# though they are explicitly included by the headers. It will cause CMake
# to try to compile them in anything that depends on fmt outside
# of this directory due to CMake's file property scope (directory based)

set(fmt_headers_sources
    fmt/format.cc
    fmt/ostream.cc
    )

blt_add_library(NAME fmt
                HEADERS ${fmt_headers}
                HEADERS_OUTPUT_SUBDIR fmt)

target_include_directories(fmt INTERFACE
            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/thirdparty/fmt>)


# Setup some variables for fmt in Axom's config.hpp
set(FMT_EXCEPTIONS FALSE PARENT_SCOPE)
set(FMT_HEADER_ONLY TRUE PARENT_SCOPE)

if( COMPILER_FAMILY_IS_CLANG )
    set(FMT_USE_EXTERN_TEMPLATES FALSE PARENT_SCOPE)
else()
    set(FMT_USE_EXTERN_TEMPLATES TRUE PARENT_SCOPE)
endif()


install(TARGETS              fmt
        EXPORT               fmt-targets
        INCLUDES DESTINATION include/fmt)
install(DIRECTORY   ${PROJECT_SOURCE_DIR}/src/thirdparty/fmt
        DESTINATION include
        )
install(EXPORT fmt-targets DESTINATION lib/cmake)

# Manually copy the souce files for the headers to the build directory
blt_copy_headers_target(
        NAME        fmt-headers-sources
        HEADERS     ${fmt_headers_sources}
        DESTINATION $<BUILD_INTERFACE:${HEADER_INCLUDES_DIRECTORY}>/fmt)
add_dependencies(blt_copy_headers_fmt blt_copy_headers_fmt-headers-sources)


set(FMT_FOUND TRUE CACHE INTERNAL "")
mark_as_advanced(FMT_FOUND)

############################################################
# Register the built-in sparse hash library
############################################################

if (ENABLE_SPARSEHASH)

  include(CheckIncludeFileCXX)

  # Try to find the necessary include and namespace for hash functions
  # This depends on the compiler and C++ standard library version
  
  set(_found_hash FALSE)

  set(_hdr <functional>)
  set(_ns  std)
  set(_src "#include ${_hdr} \n int main() { int h = ${_ns}::hash<int>()(5)\; return 0\;}")
  axom_check_code_compiles(SOURCE_STRING "${_src}" CODE_COMPILES _found_hash)
 
  if(NOT ${_found_hash})
    set(_hdr <tr1/functional>)
    set(_ns  std::tr1)
    set(_src "#include ${_hdr} \n int main() { int h = ${_ns}::hash<int>()(5)\; return 0\;}")
    axom_check_code_compiles(SOURCE_STRING "${_src}" CODE_COMPILES _found_hash)

  endif()
    
  # Set up sparsehash, if we've found it
  if(${_found_hash})
    message(STATUS "Sparsehash configured with '${_hdr}' header")
    
    set(SPARSEHASH_HASHFUN_HEADER ${_hdr} CACHE INTERNAL "")
    set(SPARSEHASH_HASHFUN_NAMESPACE ${_ns} CACHE INTERNAL "")
    set(SPARSEHASH_FOUND TRUE CACHE INTERNAL "")

    mark_as_advanced(
      SPARSEHASH_FOUND 
      SPARSEHASH_HASHFUN_HEADER 
      SPARSEHASH_HASHFUN_NAMESPACE
    )
    
    # Explicitly list sparsehash files for copy-headers,
    # Separate into two lists to maintain sparsehash's hierarchy
    set(sparsehash_headers
      sparsehash/README
      sparsehash/dense_hash_map
      sparsehash/dense_hash_set
      sparsehash/sparse_hash_map
      sparsehash/sparse_hash_set
      sparsehash/sparsetable
      sparsehash/template_util.h
      sparsehash/type_traits.h     
     )

    set(sparsehash_headers_internal
      sparsehash/internal/densehashtable.h
      sparsehash/internal/hashtable-common.h
      sparsehash/internal/libc_allocator_with_realloc.h
      sparsehash/internal/sparseconfig.h
      sparsehash/internal/sparsehashtable.h
     )
     
    blt_add_library(NAME sparsehash
            HEADERS ${sparsehash_headers}
            HEADERS_OUTPUT_SUBDIR sparsehash)
     
    install(TARGETS              sparsehash
            EXPORT               sparsehash-targets
            INCLUDES DESTINATION include)

    install(DIRECTORY   ${PROJECT_SOURCE_DIR}/src/thirdparty/sparsehash
            DESTINATION include)

    install(EXPORT sparsehash-targets DESTINATION lib/cmake)

    # Manually copy the souce files for the internal headers to the build directory
    blt_copy_headers_target(
            NAME        sparsehash-internal-headers-sources
            HEADERS     ${sparsehash_headers_internal}
            DESTINATION $<BUILD_INTERFACE:${HEADER_INCLUDES_DIRECTORY}>/sparsehash/internal)
    add_dependencies(blt_copy_headers_sparsehash blt_copy_headers_sparsehash-internal-headers-sources)
  else()
    message(STATUS "Could not configure sparsehash")
  endif()
endif()

################################
# Add tests
################################
if (ENABLE_TESTS)
    add_subdirectory(tests)
endif()
