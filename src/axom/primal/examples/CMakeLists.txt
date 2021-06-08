# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Primal examples
#------------------------------------------------------------------------------

set(primal_example_depends
        axom
        fmt
        )

blt_list_append( TO primal_example_depends ELEMENTS cuda IF ${ENABLE_CUDA} )

blt_add_executable(
    NAME        primal_introduction_ex
    SOURCES     primal_introduction.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${primal_example_depends}
    FOLDER      axom/primal/examples
    )
