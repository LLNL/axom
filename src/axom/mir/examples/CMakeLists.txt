# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

set( mir_examples
     mir_tutorial_simple.cpp
     mirConcentricCircles.cpp
   )

set( mir_example_dependencies
     #core
     mir
     #slic
   )

blt_list_append( TO mir_example_dependencies ELEMENTS sidre conduit IF AXOM_MIR_USE_SIDRE )
blt_list_append( TO mir_example_dependencies ELEMENTS raja IF RAJA_FOUND )
blt_list_append( TO mir_example_dependencies ELEMENTS openmp IF ENABLE_OPENMP )
blt_list_append( TO mir_example_dependencies ELEMENTS cuda IF ENABLE_CUDA )

foreach( example ${mir_examples} )

   get_filename_component( example_name ${example} NAME_WE )

   blt_add_executable(
        NAME       ${example_name}_ex
        SOURCES    ${example}
        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON ${mir_example_dependencies}
        FOLDER     axom/mir/examples
        )

endforeach()
