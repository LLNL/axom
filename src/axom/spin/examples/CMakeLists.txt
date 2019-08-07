# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

set(spin_example_depends
        core
        slam
        primal
        spin
        fmt
        )

    #blt_list_append( TO spin_example_depends ELEMENTS raja IF ${RAJA_FOUND} )
    #blt_list_append( TO spin_example_depends ELEMENTS umpire IF ${UMPIRE_FOUND} )
    #blt_list_append( TO spin_example_depends ELEMENTS openmp  IF ${ENABLE_OPENMP} )
blt_list_append( TO spin_example_depends ELEMENTS cuda IF ${ENABLE_CUDA} )

blt_add_executable(
    NAME        spin_introduction_ex
    SOURCES     spin_introduction.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${spin_example_depends}
    FOLDER      axom/spin/examples
    )

