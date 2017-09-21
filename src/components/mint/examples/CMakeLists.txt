#-------------------------------------------------------------------------------
# Copyright (c) 2017, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------

set( mint_examples
     mint_heat_equation.cpp
   )

foreach( example ${mint_examples} )

   get_filename_component( example_name ${example} NAME_WE )

   blt_add_executable(
        NAME       ${example_name}_ex
        SOURCES    ${example}
        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON axom_utils mint slic
        )

endforeach()
