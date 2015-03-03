
#
# List of example sources
#
set(example_sources
    ShockTube.cc
    )

#
# Build the example sources
#
foreach(example_source ${example_sources})
    make_executable( EXECUTABLE_SOURCE ${example_source} DEPENDS_ON datastore )
endforeach()