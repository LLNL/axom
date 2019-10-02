# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

set( curves_examples
     curves_basic.cpp
   )

set( curves_example_dependencies
     core
     slic
     curves
   )

foreach( example ${curves_examples} )

   get_filename_component( example_name ${example} NAME_WE )

   blt_add_executable(
        NAME       ${example_name}_ex
        SOURCES    ${example}
        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON ${curves_example_dependencies}
        FOLDER     axom/curves/examples
        )

endforeach()