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
# Specify list of tests
#

set(gtest_utils_tests
   utils_TaskTimer.cpp
   utils_endianness.cpp
   utils_fileUtilities.cpp
   utils_types.cpp
   utils_utilities.cpp
   utils_about.cpp

   axom_array.cpp

   numerics_determinants.cpp
   numerics_eigen_solve.cpp
   numerics_eigen_sort.cpp
   numerics_jacobi_eigensolve.cpp
   numerics_linear_solve.cpp
   numerics_lu.cpp
   numerics_matrix.cpp
   numerics_matvecops.cpp
   numerics_polynomial_solvers.cpp
   )


set(utils_tests_depends
   core
   gtest
   )


#
# Add gtest based tests
#
foreach(test ${gtest_utils_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}_test
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${utils_tests_depends}
                        FOLDER axom/core/tests )
    blt_add_test( NAME ${test_name}
                  COMMAND ${test_name}_test )
endforeach()



#
# Add a test to query properties of the configuration
#
set(utils_config_test_depends ${utils_tests_depends})

if(MFEM_FOUND)
    list(APPEND utils_config_test_depends mfem)
endif()

blt_add_executable( NAME utils_config_test
                    SOURCES utils_config.cpp
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON ${utils_config_test_depends}
                    FOLDER axom/core/tests )

blt_add_test( NAME utils_config
              COMMAND utils_config_test )

# Add flag to ignore unknown openmp pragmas in utils_config
set_property(TARGET utils_config_test
    APPEND_STRING PROPERTY COMPILE_FLAGS "${AXOM_DISABLE_OMP_PRAGMA_WARNINGS}")

# Have to ignore unused parameters for mfem
if(MFEM_FOUND)
    set_property(TARGET utils_config_test
        APPEND_STRING PROPERTY COMPILE_FLAGS "${AXOM_DISABLE_UNUSED_PARAMETER_WARNINGS}")
endif()


if(ENABLE_OPENMP)
    # Add a config test that uses OpenMP
    blt_add_executable( NAME utils_config_openmp_test
                        SOURCES utils_config.cpp
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${utils_config_test_depends} openmp
                        FOLDER axom/core/tests )

    blt_add_test( NAME utils_config_openmp
                  COMMAND utils_config_openmp_test )

    if(MFEM_FOUND)
        set_property(TARGET utils_config_openmp_test
            APPEND_STRING PROPERTY COMPILE_FLAGS "${AXOM_DISABLE_UNUSED_PARAMETER_WARNINGS}")
    endif()

endif()



