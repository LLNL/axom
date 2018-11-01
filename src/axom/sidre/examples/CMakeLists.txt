#------------------------------------------------------------------------------
# Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
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
#------------------------------------------------------------------------------

#
# List of single source file examples
#
set(example_sources
    sidre_stressgroups.cpp
    sidre_shocktube.cpp
)

if(ENABLE_MPI)
    set (example_parallel_sources
         sidre_createdatastore.cpp)
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


set(sidre_example_depends
    core
    sidre 
    conduit 
    conduit_blueprint 
    conduit_relay 
    slic
    )

if(HDF5_FOUND)
    list(APPEND sidre_example_depends hdf5)
endif()

if(MPI_FOUND)
    list(APPEND sidre_example_depends mpi)
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
        DEPENDS_ON  ${sidre_example_depends} ${EXTRA_LIBS}
        FOLDER axom/sidre/examples )
    blt_add_test(
        NAME ${exe_name}
        COMMAND ${exe_name}_ex )
endforeach()

if(ENABLE_MPI)
    foreach(example_source ${example_parallel_sources})
        get_filename_component(exe_name ${example_source} NAME_WE)
        blt_add_executable(
            NAME ${exe_name}_ex
            SOURCES ${example_source}
            OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
            DEPENDS_ON  ${sidre_example_depends} ${EXTRA_LIBS}
            FOLDER axom/sidre/examples)

        blt_add_test( NAME ${exe_name} 
                      COMMAND ${exe_name}_ex
                      NUM_MPI_TASKS 1
                      )
    endforeach()
endif()

if(ENABLE_FORTRAN)
    foreach(example_source ${F_example_sources})
        get_filename_component(exe_name ${example_source} NAME_WE)
        blt_add_executable(
            NAME ${exe_name}_ex
            SOURCES ${example_source}
            OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
            DEPENDS_ON ${sidre_example_depends} ${EXTRA_LIBS}
            FOLDER axom/sidre/examples )
        blt_add_test(
            NAME ${exe_name}
            COMMAND ${exe_name}_ex )
    endforeach()
    set_target_properties(${exe_name}_ex PROPERTIES LINKER_LANGUAGE Fortran)
endif()

# Note: Lulesh disabled on windows builds due to usage of unix APIs.
if(NOT WIN32)
  add_subdirectory(lulesh2)
endif()

if(ENABLE_MPI)
    add_subdirectory(spio)
endif()
