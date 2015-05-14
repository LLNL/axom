
#
# List of example sources
#
set(example_sources
	meshapiShockTube.cc
    UnstructMeshField.cc
    )

#
# Build the example sources
#
foreach(example_source ${example_sources})
    make_executable(EXECUTABLE_SOURCE ${example_source} DEPENDS_ON meshapi )
endforeach()
