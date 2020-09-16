# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Axom third-party builtin libs are all currently header only.
# We have tests here that make sure they link and can be used.
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Register the sol library for parsing lua
#------------------------------------------------------------------------------

if (LUA_FOUND AND (NOT BLT_CXX_STD STREQUAL "c++11"))
  blt_add_library(NAME    sol
                  HEADERS sol/sol.hpp)

  target_include_directories(sol SYSTEM INTERFACE
              $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/thirdparty>)

  install(TARGETS              sol
          EXPORT               axom-targets
          INCLUDES DESTINATION include)
  install(DIRECTORY   ${PROJECT_SOURCE_DIR}/thirdparty/sol
          DESTINATION include )
  install(EXPORT axom-targets DESTINATION lib/cmake)

  set(SOL_FOUND TRUE CACHE INTERNAL "")
else()
  set(SOL_FOUND FALSE CACHE INTERNAL "")
endif()
mark_as_advanced(SOL_FOUND)

message(STATUS "SOL Enabled: ${SOL_FOUND}")

#------------------------------------------------------------------------------
# Register the built-in CLI11 library for command line option parsing
#------------------------------------------------------------------------------

blt_add_library(NAME cli11
                HEADERS CLI11/CLI11.hpp)

target_include_directories(cli11 SYSTEM INTERFACE
            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/thirdparty>)

install(TARGETS              cli11
        EXPORT               axom-targets
        INCLUDES DESTINATION include)
install(DIRECTORY   ${PROJECT_SOURCE_DIR}/thirdparty/CLI11
        DESTINATION include )
install(EXPORT axom-targets DESTINATION lib/cmake)

set(CLI11_FOUND TRUE CACHE INTERNAL "")
mark_as_advanced(CLI11_FOUND)

#------------------------------------------------------------------------------
# Register the built-in fmt library for string formatting
#------------------------------------------------------------------------------

set(fmt_headers
    fmt/fmt.hpp
    fmt/core.h
    fmt/format.h
    fmt/format-inl.h
    fmt/ostream.h
    fmt/printf.h
    fmt/ranges.h
    fmt/chrono.h
    fmt/color.h
    fmt/locale.h
    fmt/time.h
    )

# Do not include the following source files in the headers list even
# though they are explicitly included by the headers. It will cause CMake
# to try to compile them in anything that depends on fmt outside
# of this directory due to CMake's file property scope (directory based)

set(fmt_headers_sources
    fmt/format.cc
    )

blt_add_library(NAME fmt
                HEADERS ${fmt_headers})

target_include_directories(fmt SYSTEM INTERFACE
            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/thirdparty>)

# Setup some variables for fmt in Axom's config.hpp
set(FMT_EXCEPTIONS FALSE PARENT_SCOPE)
set(FMT_HEADER_ONLY TRUE PARENT_SCOPE)

install(TARGETS              fmt
        EXPORT               axom-targets
        INCLUDES DESTINATION include)
install(DIRECTORY   ${PROJECT_SOURCE_DIR}/thirdparty/fmt
        DESTINATION include
        )
install(EXPORT axom-targets DESTINATION lib/cmake)

set(FMT_FOUND TRUE CACHE INTERNAL "")
mark_as_advanced(FMT_FOUND)

#------------------------------------------------------------------------------
# Register the built-in sparse hash library
#------------------------------------------------------------------------------

if (AXOM_ENABLE_SPARSEHASH)

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
            HEADERS ${sparsehash_headers})
     
    target_include_directories(sparsehash SYSTEM INTERFACE
                $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/thirdparty>)

    # Disable warning introduced in gcc@8.1+ related to how sparsehash casts memory
    # The double guarding for compiler family helps when compiling libraries with GNU
    # then compiling the main code with another compiler (this is happens for example when
    # running clang-query)
    if(C_COMPILER_FAMILY_IS_GNU AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "8.1")
        blt_add_target_compile_flags(TO sparsehash FLAGS
            $<$<AND:$<CXX_COMPILER_ID:GNU>,$<COMPILE_LANGUAGE:CXX>>:-Wno-class-memaccess>)
    endif()

    install(TARGETS              sparsehash
            EXPORT               axom-targets
            INCLUDES DESTINATION include)

    install(DIRECTORY   ${PROJECT_SOURCE_DIR}/thirdparty/sparsehash
            DESTINATION include)

    install(EXPORT axom-targets DESTINATION lib/cmake)

  else()
    message(STATUS "Could not configure sparsehash")
  endif()
endif()

#------------------------------------------------------------------------------
# Add tests
#------------------------------------------------------------------------------
if (AXOM_ENABLE_TESTS)
    add_subdirectory(tests)
endif()
