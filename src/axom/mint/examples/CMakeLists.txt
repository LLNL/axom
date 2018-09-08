#------------------------------------------------------------------------------
# Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
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

set( mint_examples
     mint_curvilinear_mesh.cpp
     mint_heat_equation.cpp
     mint_nbody_solver.cpp
     mint_particle_mesh.cpp
     mint_rectilinear_mesh.cpp
     mint_uniform_mesh.cpp
     mint_unstructured_mixed_topology_mesh.cpp
     mint_unstructured_single_topology_mesh.cpp
   )

set( mint_dependencies
     core
     mint
     slic
   )

if ( ${RAJA_FOUND} )
    list(APPEND mint_dependencies raja)
endif()

if ( ${ENABLE_OPENMP} )
    list(APPEND mint_dependencies openmp)
endif()

foreach( example ${mint_examples} )

   get_filename_component( example_name ${example} NAME_WE )

   blt_add_executable(
        NAME       ${example_name}_ex
        SOURCES    ${example}
        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON ${mint_dependencies}
        FOLDER     axom/mint/examples
        )

endforeach()
