# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
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
    list(APPEND example_sources IOSCRRead.cpp)
    list(APPEND example_sources IOSCRWrite.cpp)
    list(APPEND example_sources IOSCRCkpt.cpp)
    list(APPEND example_sources IOSCROutput.cpp)
endif() 

set(spio_example_depends sidre conduit::conduit slic ${EXTRA_LIBS})
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

