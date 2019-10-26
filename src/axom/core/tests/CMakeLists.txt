# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#
# Specify list of tests
#

set(gtest_utils_tests
   core_array.cpp
   core_execution_space.cpp
   core_memory_management.cpp
   core_stack_array.cpp

   utils_TaskTimer.cpp
   utils_endianness.cpp
   utils_fileUtilities.cpp
   utils_utilities.cpp
   utils_about.cpp


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

#
# When MPI is disabled, compile the core_types test with all other tests
#
if ( NOT ENABLE_MPI )
  list( APPEND gtest_utils_tests core_types.cpp )
endif()

set(utils_tests_depends
   core
   gtest
   )

blt_list_append( TO utils_tests_depends ELEMENTS umpire IF ${UMPIRE_FOUND} )
blt_list_append( TO utils_tests_depends ELEMENTS raja IF ${RAJA_FOUD} )
blt_list_append( TO utils_tests_depends ELEMENTS cuda IF ${ENABLE_CUDA} )

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
# Add MPI tests if MPI is enabled
#
if ( ENABLE_MPI )

  set( core_mpi_tests
       core_types.cpp
       )

  foreach( mpitest ${core_mpi_tests} )

    get_filename_component( mpitest_name ${mpitest} NAME_WE )

    blt_add_executable( NAME ${mpitest_name}_test
                        SOURCES ${mpitest}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${utils_tests_depends} mpi
                        FOLDER axom/core/tests )

    blt_add_test( NAME ${mpitest_name}
                  COMMAND ${mpitest_name}_test
                  NUM_MPI_TASKS 1 )

  endforeach()

endif()

#
# Add a test to query properties of the configuration
#
set(utils_config_test_depends core gtest )

blt_list_append( TO utils_config_test_depends ELEMENTS mfem IF ${MFEM_FOUND} )
blt_list_append( TO utils_config_test_depends ELEMENTS openmp IF ${ENABLE_OPENMP} )

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

