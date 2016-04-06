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

################################
# Add the slam (meshapi) sources
################################
add_subdirectory(src)

################################
# Add examples
################################
if (ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

################################
# Add tests and benchmarks
################################
if (ENABLE_TESTS)
  add_subdirectory(tests)
  if (ENABLE_BENCHMARKS)
    add_subdirectory(benchmarks)
  endif()
endif()


add_code_check_targets(uncrustify.cfg)

################################
# Add docs
################################
if (ENABLE_DOCS)
  add_subdirectory(docs)
endif()

################################
# Create CMake importable
# exports for all of our targets
################################
install(EXPORT ${PROJECT_NAME}-targets DESTINATION lib/cmake)
