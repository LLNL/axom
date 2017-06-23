#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------

## list of mint dependencies
set(mint_depends mint)

## load mint cmake macros
include(${PROJECT_SOURCE_DIR}/cmake/MintMacros.cmake)

## mint Lagrange basis functions
mint_fe_basis( NAME LAGRANGE
               HEADER_FILES
                  hexa_8.hpp
                  quad_4.hpp
                  tetra_4.hpp
                  tri_3.hpp
               )

## custom target to put the mint basis functions under mint/lagrange
## mint headers
set( mint_headers

    ## FEM classes
    fem/FEBasis.hpp
    fem/FEBasisTypes.hpp
    fem/FiniteElement.hpp
    fem/Lagrange.hpp
    fem/ShapeFunction.hpp

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
   )

## mint sources
set( mint_sources

    ## FEM classes
    fem/FiniteElement.cpp

    ## Mesh data-model
    mesh/CurvilinearMesh.cpp
    mesh/Field.cpp
    mesh/FieldData.cpp
    mesh/Mesh.cpp
    mesh/MeshCoordinates.cpp
    mesh/ParticleMesh.cpp
    mesh/RectilinearMesh.cpp
    mesh/StructuredMesh.cpp
    mesh/UniformMesh.cpp
   )

## create the mint library
blt_add_library(
    NAME mint
    SOURCES ${mint_sources}
    HEADERS ${mint_headers}
    HEADERS_OUTPUT_SUBDIR mint
   )

add_dependencies( ${mint_depends} )