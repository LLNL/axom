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

    SignedDistance.hpp

    ## Geometry Headers
    geom/BoundingBox.hpp
    geom/BVHTree.hpp
    geom/Determinants.hpp
    geom/HyperSphere.hpp
    geom/Intersection.hpp
    geom/InOutOctree.hpp
    geom/MortonIndex.hpp
    geom/NumericArray.hpp
    geom/OctreeBase.hpp
    geom/Orientation.hpp
    geom/Point.hpp
    geom/Ray.hpp
    geom/Segment.hpp
    geom/SpatialOctree.hpp
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
    stl/STLReader.hpp
    
    ## quest interface headers
    interface/quest.hpp
    
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
    stl/STLReader.cpp
    
    ## quest interface sources
    interface/quest.cpp
    
   )

set( quest_depends_on 
     common
     slic
     slam
     )

if (ENABLE_MPI)
  list(APPEND quest_headers stl/PSTLReader.hpp)
  list(APPEND quest_sources stl/PSTLReader.cpp)
endif()

if (ENABLE_FORTRAN)
  
  ## append shroud generated headers/sources to corresponding lists
  list(APPEND quest_headers interface/shroud/wrapQUEST.h)
  list(APPEND quest_sources interface/shroud/wrapQUEST.cpp)
  list(APPEND quest_sources interface/shroud/shroudrt.cpp)
  list(APPEND quest_sources interface/shroud/wrapfquest.f)
  
  add_subdirectory(interface/shroud)
  
  include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})
endif()

blt_add_library( NAME 
                    quest
                
                 SOURCES
                    "${quest_sources}"

                 HEADERS
                    "${quest_headers}"
              
                 DEPENDS_ON
                    "${quest_depends_on}"
                 )
target_include_directories( quest PRIVATE "${CMAKE_SOURCE_DIR}/thirdparty/stla_io" )
             
install( TARGETS quest DESTINATION lib EXPORT ${PROJECT_NAME}-targets )
