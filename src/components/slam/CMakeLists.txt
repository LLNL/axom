##################################################################################
# SLAM component -- Set-theoretic Lightweight API for Meshes (low level mesh API) 
##################################################################################
project(slam)


### Some additional flags that might be useful
# To avoid warnings from boost: -Wno-missing-braces  
# Other: -fno-strict-aliasing -Wextra

append_custom_compiler_flag(
            FLAGS_VAR SLAM_VECTORIZER_FLAGS
            GNU "-ftree-vectorize -fdump-tree-vect -ftree-vectorizer-verbose=7 -fdump-tree-optimized-details " #-march=native "
            )

### Need to add some extra flags to allow parameters and variables that are only used in ATK_WARNINGs
### TODO: Fix these when the macro becomes available
SET(EXTRA_FLAGS_WHEN_NO_DEBUG "${ATK_DISABLE_UNUSED_PARAMETER_WARNINGS}")

SET(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} ${EXTRA_FLAGS_WHEN_NO_DEBUG}" CACHE INTERNAL "slam") 
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${EXTRA_FLAGS_WHEN_NO_DEBUG}" CACHE INTERNAL "slam") 
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${EXTRA_FLAGS_WHEN_NO_DEBUG}" CACHE INTERNAL "slam")


################################
# Add the slam (meshapi) sources
################################
add_subdirectory(src)

################################
# Add examples
################################
if (BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

################################
# Add tests and benchmarks
################################
if (BUILD_TESTING)
  add_subdirectory(tests)
  if (ENABLE_BENCHMARK)
    add_subdirectory(benchmarks)
  endif()
endif()


add_code_check_targets(uncrustify.cfg)

################################
# Add docs
################################
if (BUILD_DOCS)
  add_subdirectory(docs)
endif()

################################
# Create CMake importable
# exports for all of our targets
################################
install(EXPORT ${PROJECT_NAME}-targets DESTINATION lib/cmake)
