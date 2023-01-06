# Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Primal examples
#------------------------------------------------------------------------------

set( primal_examples 
    primal_introduction.cpp
   )

set(primal_example_depends
        axom
        ${axom_device_depends}
        fmt
        )

foreach ( example ${primal_examples} )

    get_filename_component( example_name ${example} NAME_WE )

    blt_add_executable(
        NAME        ${example_name}_ex
        SOURCES     ${example}
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${primal_example_depends}
        FOLDER      axom/primal/examples
        )

endforeach()
