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


#
# Build the tinyHydro library
#
blt_add_library(
    NAME slamTinyHydroLib
    SOURCES "${tinyHydro_lib_sources}"  "${tinyHydro_lib_headers}"
    DEPENDS_ON
        slic
        slam
   ) 
   
   
#
# Add gtest based tests for tinyHydro
#
set(gtest_slam_tinyHydro_tests
    tests/slam_tinyHydro_unitTests.cpp
)

foreach(test ${gtest_slam_tinyHydro_tests})
    get_filename_component(test_name ${test} NAME_WE)
    blt_add_executable(NAME "${test_name}"
                       SOURCES ${test}
                       OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
                       DEPENDS_ON slamTinyHydroLib slic gtest
                       )
    blt_add_test(NAME ${test_name}
                 COMMAND ${test_name}
                 )
endforeach()



#
# Add the standalone executables for the sod example
#
blt_add_executable(
    NAME slam_tinyHydro_sod
    SOURCES tests/slam_tinyHydro_sod1DTwoPart.cpp
    OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON slamTinyHydroLib
)
blt_add_test(NAME slam_tinyHydro_sod
             COMMAND slam_tinyHydro_sod
             )


#
# Add the standalone executables for the sedov example
# Note: Not a test since the runtime is a bit high in debug mode
#
blt_add_executable(
    NAME slam_tinyHydro_sedov
    SOURCES tests/slam_tinyHydro_sedovTwoPart.cpp
    OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON slamTinyHydroLib
)
blt_add_test(NAME slam_tinyHydro_sedov
             COMMAND slam_tinyHydro_sedov
             )


# Add some extra flags to disable warnings about variables used only for SLIC_ASSERTS
# Temporarily used in inaccessible code block to get around warning.
# TODO: Replace with new macros from George..
#SET(EXTRA_WARNING_DISABLE_FLAGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable ")
#get_source_file_property(_origflags ${tinyHydro_dir}/tests/slam_tinyHydro_unitTests.cpp COMPILE_FLAGS)
#set_source_files_properties(
#         ${gtest_slam_tinyHydro_tests}
#         PROPERTIES COMPILE_FLAGS "${_origflags} ${EXTRA_WARNING_DISABLE_FLAGS}"     
#         )


    