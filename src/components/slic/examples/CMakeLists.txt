#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------

set(example_sources
    basic/cpplogger_example.cc
    basic/logging_example.cc
    multicode/driver.cpp
    )

set(mpi_example_sources
    basic/parallel_logging_example.cc
    basic/separate_file_per_rank.cc
    )

if (ENABLE_LUMBERJACK)
    list(APPEND mpi_example_sources basic/lumberjack_logging_example.cc )
endif()

foreach( example_source ${example_sources} )
    get_filename_component( example_name ${example_source} NAME_WE )
    blt_add_executable( NAME "slic_${example_name}" 
                        SOURCES ${example_source} 
                        DEPENDS_ON slic
                        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY} )
    blt_add_test( NAME "slic_${example_name}"
                  COMMAND "slic_${example_name}" )
endforeach()


if ( ENABLE_MPI )
    foreach( example_source ${mpi_example_sources} )
        get_filename_component( example_name ${example_source} NAME_WE )
        blt_add_executable( NAME "slic_mpi_${example_name}" 
                            SOURCES ${example_source} 
                            OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
                            DEPENDS_ON slic )
        blt_add_test( NAME "slic_mpi_${example_name}"
                      COMMAND "slic_mpi_${example_name}" )
    endforeach()
endif()

