# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# MultiMat examples
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Single source standalone examples
#------------------------------------------------------------------------------
set(example_sources
    calculate.cpp
	traversal.cpp
    )

set(multimat_tests_depends_on slic slam multimat)

set(multimat_tests_extra_depends)

foreach(example_source ${example_sources})
    get_filename_component(example_name ${example_source} NAME_WE)
    blt_add_executable(
        NAME        "multimat_${example_name}_ex"
        SOURCES     ${example_source}
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${multimat_tests_depends_on} ${multimat_tests_extra_depends}
        FOLDER      axom/multimat/examples
        )

    blt_add_test(
        NAME        multimat_${example_name}
        COMMAND     multimat_${example_name}_ex
        )
endforeach()
