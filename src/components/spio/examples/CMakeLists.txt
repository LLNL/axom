
#
# List of single source file examples
#
set(example_sources
    IORead.cpp
    IOWrite.cpp
)



#
# Build the example sources
#
foreach(example_source ${example_sources})
    make_executable( EXECUTABLE_SOURCE ${example_source}
        DEPENDS_ON spio sidre conduit slic ${EXTRA_LIBS} WITH_MPI)
endforeach()

