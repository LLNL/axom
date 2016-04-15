
#
# List of single source file examples
#
set(example_sources
    sidre_stressgroups_example.cpp
    sidre_shocktube_example.cpp
)


set(F_example_sources
    F_sidre_shocktube.f
)


if(APPLE)
    set(EXTRA_LIBS "")
elseif(UNIX)
    set(EXTRA_LIBS rt)
endif()

#
# Build the example sources
#
foreach(example_source ${example_sources})
    get_filename_component(exe_name ${example_source} NAME_WE)
    blt_add_executable(
        NAME ${exe_name}
        SOURCES ${example_source}
        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON sidre conduit conduit_io hdf5 slic ${EXTRA_LIBS})
    blt_add_test(
        NAME ${exe_name}
        COMMAND ${exe_name})
endforeach()

if(ENABLE_FORTRAN)
    foreach(example_source ${F_example_sources})
        get_filename_component(exe_name ${example_source} NAME_WE)
        blt_add_executable(
            NAME ${exe_name}
            SOURCES ${example_source}
            OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORy}
            DEPENDS_ON sidre conduit conduit_io hdf5 slic ${EXTRA_LIBS})
        blt_add_test(
            NAME ${exe_name}
            COMMAND ${exe_name})
    endforeach()
    set_target_properties(${exe_name} PROPERTIES LINKER_LANGUAGE Fortran)
#   set_target_properties(${exe_name} PROPERTIES Fortran_FORMAT "FREE")
endif()

add_subdirectory(lulesh2)
