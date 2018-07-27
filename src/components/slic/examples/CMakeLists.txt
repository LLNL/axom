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

set(example_sources
    basic/cpplogger.cpp
    basic/logging.cpp
    multicode/driver.cpp
    )

set(mpi_example_sources
    basic/parallel_logging.cpp
    basic/separate_file_per_rank.cpp
    )

if (AXOM_ENABLE_LUMBERJACK)
    list(APPEND mpi_example_sources basic/lumberjack_logging.cpp )
endif()

foreach( example_source ${example_sources} )
    get_filename_component( example_name ${example_source} NAME_WE )
    blt_add_executable( NAME "slic_${example_name}_ex" 
                        SOURCES ${example_source} 
                        DEPENDS_ON slic
                        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY} 
                        FOLDER axom/slic/examples )
    blt_add_test( NAME "slic_${example_name}"
                  COMMAND "slic_${example_name}_ex" )
endforeach()


if ( ENABLE_MPI )
    foreach( example_source ${mpi_example_sources} )
        get_filename_component( example_name ${example_source} NAME_WE )
        blt_add_executable( NAME "slic_mpi_${example_name}_ex" 
                            SOURCES ${example_source} 
                            OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
                            DEPENDS_ON slic 
                            FOLDER      axom/slic/examples )
        blt_add_test( NAME "slic_mpi_${example_name}"
                      COMMAND "slic_mpi_${example_name}_ex" 
                      NUM_MPI_TASKS 4)
    endforeach()
endif()

