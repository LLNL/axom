
#
# List of single source file examples
#
set(example_sources
    StressGroups.cpp
    sidre_shocktube_example.cpp
#    F_sidre_shocktube.f
    )


#
# Build the example sources
#
foreach(example_source ${example_sources})
    make_executable( EXECUTABLE_SOURCE ${example_source}
        DEPENDS_ON sidre conduit slic rt shroudrt)
endforeach()

set_target_properties(F_sidre_shocktube PROPERTIES LINKER_LANGUAGE Fortran)
#set_target_properties(F_sidre_shocktube  PROPERTIES Fortran_FORMAT "FREE")

add_subdirectory(lulesh2)
