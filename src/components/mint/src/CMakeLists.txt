#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------

## mint headers
set( mint_headers
    CellConnectivity.hpp
    CellType.hpp
    CurvilinearMesh.hpp
    Extent.hpp
    Field.hpp
    FieldData.hpp
    FieldTypes.hpp
    FieldVariable.hpp
    Mesh.hpp
    MeshCoordinates.hpp
    MeshType.hpp
    ParticleMesh.hpp
    RectilinearMesh.hpp
    StructuredMesh.hpp
    UniformMesh.hpp
    UnstructuredMesh.hpp
   )

## mint sources
set( mint_sources
    CurvilinearMesh.cpp
    Field.cpp
    FieldData.cpp
    Mesh.cpp
    MeshCoordinates.cpp
    ParticleMesh.cpp
    RectilinearMesh.cpp
    StructuredMesh.cpp
    UniformMesh.cpp
   )

## create the mint library 
blt_add_library(
    NAME mint
    SOURCES ${mint_sources}
    HEADERS ${mint_headers}
    HEADERS_OUTPUT_SUBDIR mint
    DEPENDS_ON common
   )