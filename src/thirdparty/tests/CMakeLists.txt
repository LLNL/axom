# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Collection of tests for axom's built-in third party libraries
#------------------------------------------------------------------------------

# Create cmake targets for running the tests through axom's code checks
axom_add_code_checks(PREFIX thirdparty_tests)

#------------------------------------------------------------------------------
# Smoke test for Umpire third party library
#------------------------------------------------------------------------------
if ( UMPIRE_FOUND )

  set( umpire_smoke_dependencies umpire gtest )

  blt_list_append( TO umpire_smoke_dependencies ELEMENTS mpi IF ENABLE_MPI)
  blt_list_append( TO umpire_smoke_dependencies ELEMENTS cuda_runtime IF ENABLE_CUDA)

  blt_add_executable( NAME umpire_smoke_test
                      SOURCES umpire_smoke.cpp
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${umpire_smoke_dependencies}
                      FOLDER axom/thirdparty/tests )

  axom_add_test( NAME    umpire_smoke
                 COMMAND umpire_smoke_test )

endif()

#------------------------------------------------------------------------------
# Smoke test for RAJA third party library
#------------------------------------------------------------------------------
if ( RAJA_FOUND )

  set( raja_smoke_dependencies RAJA gtest )

  blt_list_append( TO raja_smoke_dependencies ELEMENTS openmp IF ENABLE_OPENMP)
  blt_list_append( TO raja_smoke_dependencies ELEMENTS umpire IF UMPIRE_FOUND)
  blt_list_append( TO raja_smoke_dependencies ELEMENTS cuda IF ENABLE_CUDA)
  blt_list_append( TO raja_smoke_dependencies ELEMENTS mpi IF ENABLE_MPI)

  blt_add_executable( NAME raja_smoke_test
                      SOURCES raja_smoke.cpp
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${raja_smoke_dependencies}
                      FOLDER axom/thirdparty/tests )

  axom_add_test( NAME    raja_smoke
                 COMMAND raja_smoke_test )

endif()

#------------------------------------------------------------------------------
# Smoke test for hdf5 third party library
#------------------------------------------------------------------------------
if( HDF5_FOUND)

    set( hdf5_smoke_dependencies hdf5 gtest )

    blt_list_append( TO hdf5_smoke_dependencies ELEMENTS mpi IF ENABLE_MPI)

    blt_add_executable(NAME hdf5_smoke_test
                       SOURCES hdf5_smoke.cpp
                       OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                       DEPENDS_ON ${hdf5_smoke_dependencies}
                       FOLDER axom/thirdparty/tests )

    axom_add_test(NAME hdf5_smoke
                 COMMAND hdf5_smoke_test)
endif()

#------------------------------------------------------------------------------
# Smoke test for conduit third party library
#------------------------------------------------------------------------------
if (CONDUIT_FOUND)
    blt_add_executable(NAME       conduit_smoke_test
                       SOURCES    conduit_smoke.cpp
                       OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                       DEPENDS_ON conduit::conduit gtest
                       FOLDER     axom/thirdparty/tests )

    axom_add_test(NAME conduit_smoke
                 COMMAND conduit_smoke_test)

    if (ENABLE_FORTRAN)
        blt_add_executable(NAME       conduit_smoke_F_test
                           SOURCES    f_conduit_smoke.f
                           OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                           DEPENDS_ON conduit::conduit fruit
                           FOLDER     axom/thirdparty/tests )

        axom_add_test(NAME conduit_smoke_F
                     COMMAND conduit_smoke_F_test)
  endif()
endif()

#------------------------------------------------------------------------------
# Smoke test for mfem third party library
#------------------------------------------------------------------------------
if (MFEM_FOUND)
    set( mfem_smoke_dependencies mfem gtest )
    blt_list_append( TO mfem_smoke_dependencies ELEMENTS mpi IF ENABLE_MPI)
    blt_add_executable(NAME mfem_smoke_test
                       SOURCES mfem_smoke.cpp
                       OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                       DEPENDS_ON ${mfem_smoke_dependencies}
                       FOLDER axom/thirdparty/tests )
    set_property(
        TARGET mfem_smoke_test
        APPEND_STRING PROPERTY COMPILE_FLAGS "${AXOM_DISABLE_UNUSED_PARAMETER_WARNINGS}")

    axom_add_test(NAME    mfem_smoke
                  COMMAND mfem_smoke_test)

endif()

#------------------------------------------------------------------------------
# Smoke test for fmt third party library
#------------------------------------------------------------------------------
blt_add_executable(
    NAME fmt_smoke_test
    SOURCES fmt_smoke.cpp
    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
    DEPENDS_ON fmt gtest
    FOLDER axom/thirdparty/tests )
    
axom_add_test(NAME    fmt_smoke
              COMMAND fmt_smoke_test)

#------------------------------------------------------------------------------
# Smoke test for CLI11 third party library
#------------------------------------------------------------------------------
blt_add_executable(
    NAME cli11_smoke_test
    SOURCES cli11_smoke.cpp
    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
    DEPENDS_ON cli11
    FOLDER axom/thirdparty/tests )

axom_add_test(
    NAME cli11_smoke_help
    COMMAND cli11_smoke_test --help)

axom_add_test(
    NAME cli11_smoke_with_options
    COMMAND cli11_smoke_test -b -i 42 --some-float=3.14 --some-string "hello world")

#------------------------------------------------------------------------------
# Smoke test for sparsehash third party library
#------------------------------------------------------------------------------

if (SPARSEHASH_FOUND)
  blt_add_executable(
    NAME sparsehash_smoke_test
    SOURCES sparsehash_smoke.cpp
    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
    DEPENDS_ON sparsehash gtest
    FOLDER axom/thirdparty/tests )

  axom_add_test(NAME    sparsehash_smoke
                COMMAND sparsehash_smoke_test)
endif()

#------------------------------------------------------------------------------
# Smoke test for SOL third party library
#------------------------------------------------------------------------------

if (SOL_FOUND)
    blt_add_executable(
        NAME       sol_smoke_test
        SOURCES    sol_smoke.cpp
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON sol lua gtest
        FOLDER     axom/thirdparty/tests )

    axom_add_test(NAME    sol_smoke
                  COMMAND sol_smoke_test)
endif()

#------------------------------------------------------------------------------
# Smoke test for c2c third party library
#------------------------------------------------------------------------------

if (C2C_FOUND)
    blt_add_executable(
        NAME       c2c_smoke_test
        SOURCES    c2c_smoke.cpp
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON c2c gtest
        FOLDER     axom/thirdparty/tests )

    axom_add_test(NAME    c2c_smoke
                  COMMAND c2c_smoke_test)
endif()

#------------------------------------------------------------------------------
# Add compiler flag tests.
#
# Test whether custom compiler flags actually work as expected on all compilers.
# E.g. if we intend to disable a warning, the warning is actually disabled.
#
# These tests enable all warnings, treat warnings as errors and add the desired
# flag. They should fail if we have the wrong flag set for the given configuration.
# Flag variables are defined in blt's SetupCompilerOptions.cmake
#------------------------------------------------------------------------------

set( BLT_SMOKE_TESTS
        compiler_flag_uninitialized.cpp
        compiler_flag_unused_param.cpp
        compiler_flag_unused_var.cpp
        compiler_flag_omp_pragma.cpp
        compiler_flag_strict_aliasing.cpp
        compiler_flag_unused_local_typedef.cpp
        )
if (ENABLE_FORTRAN)
    list(APPEND BLT_SMOKE_TESTS compiler_flag_fortran_preprocessor.F)
endif()


foreach(test ${BLT_SMOKE_TESTS})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable(
        NAME ${test_name}_test
        SOURCES ${test}
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        USE_OPENMP False
        FOLDER axom/thirdparty/tests)

    axom_add_test(
        NAME ${test_name}
        COMMAND ${test_name}_test)

endforeach()

set( ADDL_SMOKE_FLAGS "${BLT_WARNINGS_AS_ERRORS_FLAG} ${BLT_ENABLE_ALL_WARNINGS_FLAG}")

set_target_properties(compiler_flag_uninitialized_test
    PROPERTIES COMPILE_FLAGS "${ADDL_SMOKE_FLAGS} ${AXOM_DISABLE_UNINITIALIZED_WARNINGS}")
set_target_properties(compiler_flag_unused_param_test
    PROPERTIES COMPILE_FLAGS "${ADDL_SMOKE_FLAGS} ${AXOM_DISABLE_UNUSED_PARAMETER_WARNINGS}")
set_target_properties(compiler_flag_unused_var_test
    PROPERTIES COMPILE_FLAGS "${ADDL_SMOKE_FLAGS} ${AXOM_DISABLE_UNUSED_VARIABLE_WARNINGS}")
set_target_properties(compiler_flag_omp_pragma_test
    PROPERTIES COMPILE_FLAGS "${ADDL_SMOKE_FLAGS} ${AXOM_DISABLE_OMP_PRAGMA_WARNINGS}")
set_target_properties(compiler_flag_unused_local_typedef_test
    PROPERTIES COMPILE_FLAGS "${ADDL_SMOKE_FLAGS} ${AXOM_DISABLE_UNUSED_LOCAL_TYPEDEF}")

# the aliasing test requires optimization, set this up as custom compiler flag
# in case we need more flexibility for this flag (e.g. for MSVC)
blt_append_custom_compiler_flag(FLAGS_VAR TMP_OPT_FLAG
    DEFAULT "-O2"
    MSVC    " "
)

set_target_properties(compiler_flag_strict_aliasing_test
    PROPERTIES COMPILE_FLAGS "${ADDL_SMOKE_FLAGS} ${TMP_OPT_FLAG} ${AXOM_DISABLE_ALIASING_WARNINGS}")


# The following test checks that the C preprocessor is invoked on Fortran files with extension *.F
if(ENABLE_FORTRAN)
    set_target_properties(compiler_flag_fortran_preprocessor_test
                          PROPERTIES LINKER_LANGUAGE Fortran)
endif()
