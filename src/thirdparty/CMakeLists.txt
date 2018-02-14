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
# Create SparseHash variables to be used later
# Note: Do not register or add this library because
#   it is a internal header only library
############################################################

if (ENABLE_SPARSEHASH)

  include(CheckIncludeFileCXX)

  check_include_file_cxx(tr1/functional FOUND_TR1_FUNCTIONAL)
  if(FOUND_TR1_FUNCTIONAL)
    set(SPARSEHASH_COMPILER_DEFS
           SPARSEHASH_HASHFUN_HEADER=<tr1/functional>
           SPARSEHASH_HASHFUN_NAMESPACE=std::tr1
        CACHE INTERNAL "")
    set(FOUND_FUNCTIONAL TRUE)
  else()
    check_include_file_cxx(functional FOUND_STD_FUNCTIONAL)
    if (FOUND_STD_FUNCTIONAL)
       set(SPARSEHASH_COMPILER_DEFS
              SPARSEHASH_HASHFUN_HEADER=<functional>
              SPARSEHASH_HASHFUN_NAMESPACE=std
           CACHE INTERNAL "")
       set(FOUND_FUNCTIONAL TRUE)
    endif()
  endif()

  if(FOUND_FUNCTIONAL)
     set(SPARSEHASH_FOUND TRUE CACHE INTERNAL "")
     mark_as_advanced(SPARSEHASH_FOUND)
  endif()

endif()

################################
# Add tests
################################
if (ENABLE_TESTS)
    add_subdirectory(tests)
endif()
