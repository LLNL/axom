# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Spio examples
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# List of single source file examples
#------------------------------------------------------------------------------
set(example_sources
    IORead.cpp
    IOWrite.cpp
    )

if (SCR_FOUND)
    list(APPEND example_sources IO_SCR_Read.cpp)
    list(APPEND example_sources IO_SCR_Write.cpp)
    list(APPEND example_sources IO_SCR_Checkpoint.cpp)
    list(APPEND example_sources IO_SCR_Output.cpp)
endif() 

set(spio_example_depends axom gtest ${EXTRA_LIBS})
blt_list_append(TO spio_example_depends ELEMENTS hdf5 IF HDF5_FOUND)
blt_list_append(TO spio_example_depends ELEMENTS scr  IF SCR_FOUND)

#------------------------------------------------------------------------------
# Add targets and tests for spio examples
#------------------------------------------------------------------------------
foreach(src ${example_sources})
    get_filename_component(exe_name ${src} NAME_WE)
    blt_add_executable(
        NAME spio_${exe_name}_ex
        SOURCES ${src}
        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON ${spio_example_depends}
        FOLDER axom/sidre/examples
        )
endforeach()

    if(AXOM_ENABLE_TESTS)
        if(SCR_FOUND)
            if(ENABLE_MPI)
                axom_add_test(
                    NAME spio_IO_SCR_Checkpoint
                    COMMAND spio_IO_SCR_Checkpoint_ex 
                    NUM_MPI_TASKS 4 )
            else()
                axom_add_test(
                    NAME spio_IO_SCR_Checkpoint
                    COMMAND spio_IO_SCR_Checkpoint_ex )
            endif()
        endif()
    endif()
