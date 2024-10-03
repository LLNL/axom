# Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

set( mir_example_dependencies
     core
     slic
     mir
   )

# Speed up compilation by separating out the execspaces into different runMIR_xxx.cpp files.
axom_add_executable(
    NAME       mir_concentric_circles
    SOURCES    mir_concentric_circles.cpp runMIR_seq.cpp runMIR_omp.cpp runMIR_cuda.cpp runMIR_hip.cpp
    OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON ${mir_example_dependencies}
    FOLDER     axom/mir/examples
    )
