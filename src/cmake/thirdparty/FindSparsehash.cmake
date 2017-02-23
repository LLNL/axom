###############################################################################
#
# Setup Sparsehash
# This file defines:
#  SPARSEHASH_FOUND - If Sparsehash was found
#  SPARSEHASH_INCLUDE_DIRS - The Sparsehash include directories

# first Check for SPARSEHASH_DIR

#find includes
find_path( SPARSEHASH_INCLUDE_DIRS sparsehash/type_traits.h
          PATHS  ${SPARSEHASH_DIR}/include
          NO_DEFAULT_PATH
          NO_CMAKE_ENVIRONMENT_PATH
          NO_CMAKE_PATH
          NO_SYSTEM_ENVIRONMENT_PATH
          NO_CMAKE_SYSTEM_PATH)


include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set SPARSEHASH_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(Sparsehash  DEFAULT_MSG
                                  SPARSEHASH_INCLUDE_DIRS)

if(NOT SPARSEHASH_FOUND)
    MESSAGE(WARNING "Could not find Sparsehash. Sparsehash support needs explicit SPARSEHASH_DIR")
endif()
