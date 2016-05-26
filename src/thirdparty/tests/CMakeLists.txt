#
# Specify list of tests
#

################
# Test conduit tests
################
if (CONDUIT_FOUND)
    blt_add_executable(NAME conduit_smoke
                       SOURCES conduit_smoke.cpp
                       OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                       DEPENDS_ON conduit gtest)
    blt_add_test(NAME conduit_smoke
                 COMMAND conduit_smoke)

    if (ENABLE_FORTRAN)
        blt_add_executable(NAME f_conduit_smoke
                           SOURCES f_conduit_smoke.f
                           OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                           DEPENDS_ON conduit conduit_fortran fruit)
        blt_add_test(NAME f_conduit_smoke
                     COMMAND f_conduit_smoke)
  endif()
endif()

################
# Test fmt thirdparty library
################
blt_add_executable(
    NAME fmt_smoke
    SOURCES fmt_smoke.cpp
    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
    DEPENDS_ON fmt gtest)
blt_add_test(NAME fmt_smoke COMMAND fmt_smoke)





################################################################################
# Add compiler flag tests.
#
# Test whether custom compiler flags actually work as expected on all compilers.
# E.g. if we intend to disable a warning, the warning is actually disabled. 
#
# These tests enable all warnings, treat warnings as errors and add the desired 
# flag. They should fail if we have the wrong flag set for the given configuration.  
# Flag variables are defined in blt's SetupCompilerOptions.cmake
################################################################################

set( BLT_SMOKE_TESTS
        compiler_flag_uninitialized.cpp
        compiler_flag_unused_param.cpp
        compiler_flag_unused_var.cpp
        compiler_flag_omp_pragma.cpp
        compiler_flag_strict_aliasing.cpp
        )
if (ENABLE_FORTRAN)                    
    list(APPEND BLT_SMOKE_TESTS compiler_flag_fortran_preprocessor.F)
endif()


foreach(test ${BLT_SMOKE_TESTS})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( 
        NAME ${test_name}
        SOURCES ${test}
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        USE_OPENMP False
        )
    blt_add_test( 
        NAME ${test_name} 
        COMMAND ${test_name} 
        )
endforeach()

set( ADDL_SMOKE_FLAGS "${BLT_WARNINGS_AS_ERRORS_FLAG} ${BLT_ENABLE_ALL_WARNINGS_FLAG}")

set_target_properties(compiler_flag_uninitialized 
    PROPERTIES COMPILE_FLAGS "${ADDL_SMOKE_FLAGS} ${ATK_DISABLE_UNINITIALIZED_WARNINGS}")
set_target_properties(compiler_flag_unused_param 
    PROPERTIES COMPILE_FLAGS "${ADDL_SMOKE_FLAGS} ${ATK_DISABLE_UNUSED_PARAMETER_WARNINGS}")
set_target_properties(compiler_flag_unused_var 
    PROPERTIES COMPILE_FLAGS "${ADDL_SMOKE_FLAGS} ${ATK_DISABLE_UNUSED_VARIABLE_WARNINGS}")
set_target_properties(compiler_flag_omp_pragma 
    PROPERTIES COMPILE_FLAGS "${ADDL_SMOKE_FLAGS} ${ATK_DISABLE_OMP_PRAGMA_WARNINGS}")
    
# the aliasing test requires optimization, set this up as custom compiler flag
# in case we need more flexibility for this flag (e.g. for MSVC)
blt_append_custom_compiler_flag(FLAGS_VAR TMP_OPT_FLAG DEFAULT "-O2")
set_target_properties(compiler_flag_strict_aliasing 
    PROPERTIES COMPILE_FLAGS "${ADDL_SMOKE_FLAGS} ${TMP_OPT_FLAG} ${ATK_DISABLE_ALIASING_WARNINGS}")


# The following test checks that the C preprocessor is invoked on Fortran files with extension *.F 
if(ENABLE_FORTRAN)
    set_target_properties(compiler_flag_fortran_preprocessor 
                          PROPERTIES LINKER_LANGUAGE Fortran)
endif()

