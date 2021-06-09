# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Tests for Axom's Core component
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Specify list of tests
#------------------------------------------------------------------------------

set(gtest_utils_tests
    core_array.hpp
    core_execution_for_all.hpp
    core_execution_space.hpp
    core_memory_management.hpp
    core_Path.hpp
    core_stack_array.hpp

    numerics_determinants.hpp
    numerics_eigen_solve.hpp
    numerics_eigen_sort.hpp
    numerics_floating_point_limits.hpp
    numerics_jacobi_eigensolve.hpp
    numerics_linear_solve.hpp
    numerics_lu.hpp
    numerics_matrix.hpp
    numerics_matvecops.hpp
    numerics_polynomial_solvers.hpp

    utils_Timer.hpp
    utils_endianness.hpp
    utils_fileUtilities.hpp
    utils_nvtx_settings.hpp
    utils_utilities.hpp
    utils_about.hpp
    )

set(utils_tests_depends
    axom
    gtest
    )

blt_list_append( TO utils_tests_depends ELEMENTS cuda IF ${ENABLE_CUDA} )

#------------------------------------------------------------------------------
# Serial GTests
#------------------------------------------------------------------------------

# Only enable this test in serial mode if ENABLE_MPI is off
if (NOT ENABLE_MPI)
  list(APPEND gtest_utils_tests core_types.hpp)
endif()

blt_add_executable( NAME       core_serial_tests
                    SOURCES    core_serial_main.cpp
                    HEADERS    ${gtest_utils_tests}
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON ${utils_tests_depends}
                    FOLDER     axom/core/tests )

foreach(test_suite ${gtest_utils_tests})
    get_filename_component( test_suite ${test_suite} NAME_WE )
    axom_add_test( NAME    ${test_suite}_serial
                   COMMAND core_serial_tests --gtest_filter=${test_suite}* )
endforeach()

#------------------------------------------------------------------------------
# MPI GTests
#------------------------------------------------------------------------------
if (ENABLE_MPI)
  set( core_mpi_tests
       core_types.hpp
       )

  blt_add_executable( NAME       core_mpi_tests
                      SOURCES    core_mpi_main.cpp
                      HEADERS    ${core_mpi_tests}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${utils_tests_depends} mpi
                      FOLDER     axom/core/tests )

  foreach( test_suite ${core_mpi_tests} )
    get_filename_component( test_suite ${test_suite} NAME_WE )
    axom_add_test( NAME          ${test_suite}_mpi
                   COMMAND       core_mpi_tests --gtest_filter=${test_suite}*
                   NUM_MPI_TASKS 1 )
  endforeach()
endif()

#------------------------------------------------------------------------------
# Add a test to query properties of the configuration
#------------------------------------------------------------------------------
set(utils_config_test_depends core gtest )

blt_list_append( TO utils_config_test_depends ELEMENTS mfem IF MFEM_FOUND )

blt_add_executable( NAME       utils_config_test
                    SOURCES    utils_config.cpp
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON ${utils_config_test_depends}
                    FOLDER     axom/core/tests )

axom_add_test( NAME    utils_config
               COMMAND utils_config_test )

# Add flag to ignore unknown openmp pragmas in utils_config
set_property(TARGET utils_config_test
    APPEND_STRING PROPERTY COMPILE_FLAGS "${AXOM_DISABLE_OMP_PRAGMA_WARNINGS}")

# Have to ignore unused parameters for mfem
if(MFEM_FOUND)
    set_property(TARGET utils_config_test
        APPEND_STRING PROPERTY COMPILE_FLAGS "${AXOM_DISABLE_UNUSED_PARAMETER_WARNINGS}")
endif()
