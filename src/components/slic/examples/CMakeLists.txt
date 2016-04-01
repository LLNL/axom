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
    basic/lumberjack_logging_example.cc
    basic/separate_file_per_rank.cc
    )


foreach(example_source ${example_sources})
    get_filename_component(example_name ${example_source} NAME_WE)
    make_executable(EXECUTABLE_NAME "slic_${example_name}" 
                    EXECUTABLE_SOURCE ${example_source} 
                    DEPENDS_ON slic
                    ADD_CTEST IS_EXAMPLE
                    )
endforeach()


if (ENABLE_MPI)
    foreach(example_source ${mpi_example_sources})
        get_filename_component(example_name ${example_source} NAME_WE)
        make_executable(EXECUTABLE_NAME "slic_mpi_${example_name}" 
                        EXECUTABLE_SOURCE ${example_source} 
                        DEPENDS_ON slic
                        ADD_CTEST IS_EXAMPLE
                        )
    endforeach()
endif()

