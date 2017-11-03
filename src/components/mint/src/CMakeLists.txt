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

## mint headers
set( mint_headers

    ## FEM classes
    fem/FEBasis.hpp
    fem/FEBasisTypes.hpp
    fem/FiniteElement.hpp
    fem/shape_functions/Lagrange.hpp
    fem/shape_functions/ShapeFunction.hpp

    ## Lagrange Basis
    fem/shape_functions/Lagrange.hpp
    fem/shape_functions/lagrange/lagrange_hexa_27.hpp
    fem/shape_functions/lagrange/lagrange_hexa_8.hpp
    fem/shape_functions/lagrange/lagrange_prism_6.hpp
    fem/shape_functions/lagrange/lagrange_pyra_5.hpp
    fem/shape_functions/lagrange/lagrange_quad_4.hpp
    fem/shape_functions/lagrange/lagrange_quad_9.hpp
    fem/shape_functions/lagrange/lagrange_tetra_4.hpp
    fem/shape_functions/lagrange/lagrange_tri_3.hpp

    ## Mesh data-model
    mesh/CellConnectivity.hpp
    mesh/CellType.hpp
    mesh/CurvilinearMesh.hpp
    mesh/Extent.hpp
    mesh/Field.hpp
    mesh/FieldData.hpp
    mesh/FieldTypes.hpp
    mesh/FieldVariable.hpp
    mesh/Mesh.hpp
    mesh/MeshCoordinates.hpp
    mesh/MeshType.hpp
    mesh/ParticleMesh.hpp
    mesh/RectilinearMesh.hpp
    mesh/StructuredMesh.hpp
    mesh/UniformMesh.hpp
    mesh/UnstructuredMesh.hpp

    ## utils
    utils/vtk_utils.hpp
    utils/Vector.hpp
    utils/DataTypes.hpp
   )

## mint sources
set( mint_sources

    ## FEM classes
    fem/FiniteElement.cpp

    ## Mesh data-model
    mesh/CurvilinearMesh.cpp
    mesh/Extent.cpp
    mesh/Mesh.cpp
    mesh/MeshCoordinates.cpp
    mesh/ParticleMesh.cpp
    mesh/RectilinearMesh.cpp
    mesh/StructuredMesh.cpp
    mesh/UniformMesh.cpp

    ## utils
    utils/vtk_utils.cpp
   )

## create the mint library
blt_add_library(
    NAME mint
    SOURCES ${mint_sources}
    HEADERS ${mint_headers}
    HEADERS_OUTPUT_SUBDIR mint
    DEPENDS_ON axom_utils slic
    FOLDER axom/mint
   )

install(TARGETS              mint
        EXPORT               mint-targets
        DESTINATION          lib
        INCLUDES DESTINATION include/mint
        )
install(FILES       ${mint_headers}
        DESTINATION include/mint
        )
