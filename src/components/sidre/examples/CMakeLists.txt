
#
# List of example sources
#
set(example_sources
    ShockTube.cpp
    DumpStateFile.cpp
    StressGroups.cpp
    )

#
# Build the example sources
#
foreach(example_source ${example_sources})
    make_executable( EXECUTABLE_SOURCE ${example_source} DEPENDS_ON sidre
conduit rt )
endforeach()
