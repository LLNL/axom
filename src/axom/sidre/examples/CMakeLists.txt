# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Sidre and Spio examples
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# List of single source file examples
#------------------------------------------------------------------------------
set(example_sources
    sidre_stressgroups.cpp
    sidre_shocktube.cpp
)

if(ENABLE_MPI)
    set (example_parallel_sources
         sidre_array.cpp
         sidre_createdatastore.cpp
         sidre_external_array.cpp
         sidre_generateindex.cpp)
else()
    list(APPEND example_sources sidre_createdatastore.cpp)
endif()

set(F_example_sources
    sidre_shocktube_F.f
)

if(APPLE)
    set(EXTRA_LIBS "")
elseif(UNIX)
    set(EXTRA_LIBS rt)
endif()

if(MFEM_FOUND AND AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION)
    list(APPEND example_sources
         sidre_mfem_datacollection_vis.cpp
         sidre_mfem_datacollection_restart.cpp
         sidre_mfem_datacollection_materials.cpp)
endif()

#------------------------------------------------------------------------------
# Add targets and tests for serial examples
#------------------------------------------------------------------------------
foreach(example_source ${example_sources})
    get_filename_component(exe_name ${example_source} NAME_WE)
    blt_add_executable(
        NAME ${exe_name}_ex
        SOURCES ${example_source}
        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  axom ${EXTRA_LIBS}
        FOLDER axom/sidre/examples )

    if(AXOM_ENABLE_TESTS)
        axom_add_test(
            NAME ${exe_name}
            COMMAND ${exe_name}_ex )
    endif()
endforeach()

#------------------------------------------------------------------------------
# Add targets and tests for parallel examples
#------------------------------------------------------------------------------
if(ENABLE_MPI)
    foreach(example_source ${example_parallel_sources})
        get_filename_component(exe_name ${example_source} NAME_WE)
        blt_add_executable(
            NAME ${exe_name}_ex
            SOURCES ${example_source}
            OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
            DEPENDS_ON  axom ${EXTRA_LIBS}
            FOLDER axom/sidre/examples)

        if(AXOM_ENABLE_TESTS)
            axom_add_test( NAME ${exe_name}
                           COMMAND ${exe_name}_ex
                           NUM_MPI_TASKS 1
                           )
        endif()
    endforeach()

    blt_add_test( NAME sidre_createdatastore_parallel
                  COMMAND sidre_createdatastore_ex
                  NUM_MPI_TASKS 4
                  )

    blt_add_test( NAME sidre_generateindex_parallel
                  COMMAND sidre_generateindex_ex
                  NUM_MPI_TASKS 4
                  )

endif()

#------------------------------------------------------------------------------
# Add targets and tests for Fortran examples
#------------------------------------------------------------------------------
if(ENABLE_FORTRAN)
    foreach(example_source ${F_example_sources})
        get_filename_component(exe_name ${example_source} NAME_WE)
        blt_add_executable(
            NAME ${exe_name}_ex
            SOURCES ${example_source}
            OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
            DEPENDS_ON sidre ${EXTRA_LIBS}
            FOLDER axom/sidre/examples )

        if(AXOM_ENABLE_TESTS)
            axom_add_test(
                NAME ${exe_name}
                COMMAND ${exe_name}_ex )
        endif()
    endforeach()
endif()

#------------------------------------------------------------------------------
# Add Spio examples
#------------------------------------------------------------------------------
if(ENABLE_MPI)
    add_subdirectory(spio)
endif()

#------------------------------------------------------------------------------
# Sidre version of lulesh 2.0
# Note: Disabled on windows builds due to Lulesh's usage of unix APIs.
#------------------------------------------------------------------------------
if(NOT WIN32)
  add_subdirectory(lulesh2)
endif()
