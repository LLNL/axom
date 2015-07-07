
#
# List of single source file examples
#
set(example_sources
    sidre_shocktube_example.cpp
    )

#
# Build the example sources
#
foreach(example_source ${example_sources})
    make_executable( EXECUTABLE_SOURCE ${example_source} DEPENDS_ON sidre conduit slic )
endforeach()

add_subdirectory(lulesh2)
