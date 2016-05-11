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



