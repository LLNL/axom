# Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

set( mir_example_dependencies
     core
     slic
     mir
   )

axom_add_executable(
    NAME       mir_tutorial_simple_ex
    SOURCES    mir_tutorial_simple.cpp
    OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON ${mir_example_dependencies}
    FOLDER     axom/mir/examples
    )
