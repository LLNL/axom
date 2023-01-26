# Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Examples for Spin component
#------------------------------------------------------------------------------

set(spin_example_depends
        axom
        ${axom_device_depends}
        fmt
        )

blt_add_executable(
    NAME        spin_introduction_ex
    SOURCES     spin_introduction.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${spin_example_depends}
    FOLDER      axom/spin/examples
    )

# Note: The following example uses quest to read in an STL file
#       and mint for meshing
if (AXOM_ENABLE_MINT AND AXOM_ENABLE_QUEST AND RAJA_FOUND AND UMPIRE_FOUND)

    blt_add_executable(
        NAME        spin_bvh_two_pass_ex
        SOURCES     spin_bvh_two_pass.cpp
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${spin_example_depends}
        FOLDER      axom/spin/examples
        )

endif()
