# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#
# List of single source file examples
#
set(example_sources
    IORead.cpp
    IOWrite.cpp
    )

if (SCR_FOUND)
    list(APPEND example_sources IOSCRRead.cpp)
    list(APPEND example_sources IOSCRWrite.cpp)
endif() 

set(spio_example_depends sidre conduit conduit_relay slic ${EXTRA_LIBS})

if (HDF5_FOUND)
    list(APPEND spio_example_depends hdf5)
endif()

if (SCR_FOUND)
    list(APPEND spio_example_depends scr)
endif()


#
# Build the example sources
#
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

