#------------------------------------------------------------------------------
# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
#------------------------------------------------------------------------------

###################
# MultiMat examples
###################

#
# Single source standalone examples
#
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
