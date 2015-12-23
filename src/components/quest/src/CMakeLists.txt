#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------

set( quest_headers
   
    ## Geometry Headers
    geom/BoundingBox.hpp
    geom/Determinants.hpp
    geom/HyperSphere.hpp
    geom/Orientation.hpp
    geom/Point.hpp
    geom/Segment.hpp
    geom/SquaredDistance.hpp
    geom/Triangle.hpp
    geom/Vector.hpp
    geom/fuzzy_compare.hpp
    
    ## Mesh Library Headers
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
    
    ## STL library headers
    stl/stla_io.hpp
    stl/STLReader.hpp
    
   )
   
set( quest_sources
    
    ## CG Sources
    # <nothing yet>
        
    ## Mesh Library Sources
    mesh/CurvilinearMesh.cpp
    mesh/Field.cpp
    mesh/FieldData.cpp
    mesh/Mesh.cpp
    mesh/MeshCoordinates.cpp
    mesh/ParticleMesh.cpp
    mesh/RectilinearMesh.cpp
    mesh/StructuredMesh.cpp
    mesh/UniformMesh.cpp
    
    ## STL library sources
    stl/stla_io.cpp
    stl/STLReader.cpp
    
   )
   
copy_headers_target( ${PROJECT_NAME} "${quest_headers}"
                     ${HEADER_INCLUDES_DIRECTORY}/${PROJECT_NAME} )
                     
make_library( LIBRARY_NAME 
                quest
                
              LIBRARY_SOURCES
                "${quest_sources}"
                "${quest_headers}"
              
              DEPENDS_ON
                common
                slic
                #slam
             )
             
install( TARGETS quest DESTINATION lib EXPORT ${PROJECT_NAME}-targets )