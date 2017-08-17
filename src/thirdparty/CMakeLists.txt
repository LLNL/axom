#
# Our third-party builtin libs are all currently header only,
# so here we are adding tests that make sure they link and can 
# be used
#

############################################################
# Register the fmt header-only library for string formatting 
############################################################

# Do not include these in the headers list even though they
# are included in the headers, it will cause CMake to try to
# compile them in anything that depends on fmt outside of this
# directory due to CMake's file property scope (directory based)
#set(fmt_sources
#    fmt/format.cc
#    fmt/ostream.cc
#    )

set(fmt_headers
    fmt/format.h
    fmt/ostream.h 
    fmt/printf.h    
    )

blt_add_library(NAME fmt
                HEADERS ${fmt_headers}
                HEADERS_OUTPUT_SUBDIR fmt)

target_include_directories(fmt INTERFACE
            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/thirdparty/fmt>)
target_compile_definitions(fmt INTERFACE 
            FMT_EXCEPTIONS=0
            FMT_HEADER_ONLY=1)
if( COMPILER_FAMILY_IS_CLANG )
    target_compile_definitions(fmt INTERFACE 
                FMT_USE_EXTERN_TEMPLATES=0)
endif()

install(TARGETS              fmt
        EXPORT               fmt-targets
        INCLUDES DESTINATION include/fmt)
install(DIRECTORY   ${PROJECT_SOURCE_DIR}/thirdparty/fmt
        DESTINATION include
        )
install(EXPORT fmt-targets DESTINATION lib/cmake)

set(FMT_FOUND TRUE CACHE INTERNAL "")
mark_as_advanced(FMT_FOUND)

############################################################
# Create SparseHash variables to be used later
# Note: Do not register or add this library because
#   it is a internal header only library
############################################################

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

################################
# Add tests
################################
if (ENABLE_TESTS)
    add_subdirectory(tests)
endif()
