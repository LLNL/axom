####################################
# Slam (MeshAPI) version of tinyHydro
####################################

set(tinyHydro_lib_headers
    TinyHydroTypes.hpp
    HydroC.hpp
    Part.hpp
    State.hpp
    VectorXY.hpp
    PolygonMeshXY.hpp
    )

set(tinyHydro_lib_sources
    TinyHydroTypes.cpp
    HydroC.cpp
    Part.cpp
    State.cpp
    PolygonMeshXY.cpp
    )



# blt_add_executable depends_on handles our include dir settings
# these teste are using boost so boost needs to be a dep
set(slam_tiny_hydro_lib_depends_on axom_common slam slic gtest)
if (BOOST_FOUND)
    list(APPEND slam_tiny_hydro_lib_depends_on boost)
endif()



#
# Build the tinyHydro library
#
blt_add_library(
    NAME        slamTinyHydroLib
    SOURCES     ${tinyHydro_lib_sources}
                ${tinyHydro_lib_headers}
    DEPENDS_ON  ${slam_tiny_hydro_lib_depends_on}
   )


set(slam_tiny_hydro_tests_depends_on ${slam_tiny_hydro_lib_depends_on})
list(APPEND slam_tiny_hydro_tests_depends_on slamTinyHydroLib)
#

#
# Add gtest based tests for tinyHydro
#
blt_add_executable(
    NAME        slam_tinyHydro_unitTests
    SOURCES     tests/slam_tinyHydro_unitTests.cpp
    OUTPUT_DIR  ${TEST_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${slam_tiny_hydro_tests_depends_on}
    )
blt_add_test(
    NAME        slam_tinyHydro_unitTests
    COMMAND     slam_tinyHydro_unitTests
    )


# Add the standalone executables for the sod example
#
blt_add_executable(
    NAME        slam_tinyHydro_sod
    SOURCES     tests/slam_tinyHydro_sod1DTwoPart.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${slam_tiny_hydro_tests_depends_on}
)
blt_add_test(
    NAME        slam_tinyHydro_sod
    COMMAND     slam_tinyHydro_sod
    )


#
# Add the standalone executables for the sedov example
# Note: Not a test since the runtime is a bit high in debug mode (~40s)
#
blt_add_executable(
    NAME        slam_tinyHydro_sedov
    SOURCES     tests/slam_tinyHydro_sedovTwoPart.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${slam_tiny_hydro_tests_depends_on}
)
