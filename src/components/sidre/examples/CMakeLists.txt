
#
# List of single source file examples
#
set(example_sources
    StressGroups.cpp
    shocktube_example.cpp
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
    make_executable(
        EXECUTABLE_NAME   "sidre_${exe_name}"
        EXECUTABLE_SOURCE ${example_source}
        DEPENDS_ON sidre conduit slic ${EXTRA_LIBS}
        ADD_CTEST)
endforeach()

if(ENABLE_FORTRAN)
    foreach(example_source ${F_example_sources})
        get_filename_component(exe_name ${example_source} NAME_WE)
        make_executable(
            EXECUTABLE_NAME   ${exe_name}
            EXECUTABLE_SOURCE ${example_source}
            DEPENDS_ON sidre conduit slic shroudrt ${EXTRA_LIBS}
            ADD_CTEST)
    endforeach()
    set_target_properties(${exe_name} PROPERTIES LINKER_LANGUAGE Fortran)
#   set_target_properties(${exe_name} PROPERTIES Fortran_FORMAT "FREE")
endif()

add_subdirectory(lulesh2)
