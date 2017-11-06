###############################################################################
# Copyright (c) 2017, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
#
##############################################################################

#
# List of single source file examples
#
set(example_sources
    sidre_stressgroups.cpp
    sidre_shocktube.cpp
)


set(F_example_sources
    sidre_shocktube_F.f
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
        NAME ${exe_name}_ex
        SOURCES ${example_source}
        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON sidre conduit conduit_relay hdf5 slic ${EXTRA_LIBS})
    blt_add_test(
        NAME ${exe_name}
        COMMAND ${exe_name}_ex )
endforeach()

if(ENABLE_FORTRAN)
    foreach(example_source ${F_example_sources})
        get_filename_component(exe_name ${example_source} NAME_WE)
        blt_add_executable(
            NAME ${exe_name}_ex
            SOURCES ${example_source}
            OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
            DEPENDS_ON sidre conduit conduit_relay hdf5 slic ${EXTRA_LIBS})
        blt_add_test(
            NAME ${exe_name}
            COMMAND ${exe_name}_ex )
    endforeach()
    set_target_properties(${exe_name}_ex PROPERTIES LINKER_LANGUAGE Fortran)
#   set_target_properties(${exe_name} PROPERTIES Fortran_FORMAT "FREE")
endif()

add_subdirectory(lulesh2)
