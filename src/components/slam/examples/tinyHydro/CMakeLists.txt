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
make_library( 
    LIBRARY_NAME slamTinyHydroLib 
    LIBRARY_SOURCES "${tinyHydro_lib_sources}"  "${tinyHydro_lib_headers}"
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
    add_gtest( TEST_SOURCE ${test} DEPENDS_ON slamTinyHydroLib slic)
endforeach()



#
# Add the standalone executables for the sod example
#
make_executable( 
    EXECUTABLE_NAME "slam_tinyHydro_sod" 
    EXECUTABLE_SOURCE tests/slam_tinyHydro_sod1DTwoPart.cpp
    DEPENDS_ON slamTinyHydroLib
    IS_EXAMPLE
    ADD_CTEST
)

#
# Add the standalone executables for the sedov example
# Note: Not a test since the runtime is a bit high in debug mode
#
make_executable(
    EXECUTABLE_NAME "slam_tinyHydro_sedov" 
    EXECUTABLE_SOURCE tests/slam_tinyHydro_sedovTwoPart.cpp
    DEPENDS_ON slamTinyHydroLib
    IS_EXAMPLE
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


    