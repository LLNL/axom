# Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Sina examples
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# List of single source file examples
#------------------------------------------------------------------------------
set(sina_example_sources
    sina_basic.cpp
    sina_curve_set.cpp
    sina_tutorial.cpp
)

set(sina_example_depends sina conduit::conduit)

#------------------------------------------------------------------------------
# Add targets for Sina examples
#------------------------------------------------------------------------------
foreach(src ${sina_example_sources})
    get_filename_component(exe_name ${src} NAME_WE)
    axom_add_executable(
        NAME       ${exe_name}_ex
        SOURCES    ${src}
        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON ${sina_example_depends}
        FOLDER     axom/sina/examples
        )
endforeach()
