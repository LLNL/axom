
#
# List of single source file examples
#
set(example_sources
    StressGroups.cpp
    sidre_shocktube_example.cpp
    F_sidre_shocktube.f
    )


#
# Build the example sources
#
foreach(example_source ${example_sources})
    make_executable( EXECUTABLE_SOURCE ${example_source} DEPENDS_ON sidre
conduit slic rt )
endforeach()

add_subdirectory(lulesh2)
