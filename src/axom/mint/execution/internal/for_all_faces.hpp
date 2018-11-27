/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifndef MINT_FOR_ALL_FACES_HPP_
#define MINT_FOR_ALL_FACES_HPP_

// mint includes
#include "axom/mint/execution/xargs.hpp"        // for xargs

#include "axom/mint/config.hpp"                 // for compile-time definitions
#include "axom/mint/mesh/Mesh.hpp"              // for Mesh
#include "axom/mint/mesh/StructuredMesh.hpp"    // for StructuredMesh
#include "axom/mint/mesh/UniformMesh.hpp"       // for UniformMesh
#include "axom/mint/mesh/RectilinearMesh.hpp"   // for RectilinearMesh
#include "axom/mint/mesh/CurvilinearMesh.hpp"   // for CurvilinearMesh
#include "axom/mint/execution/policy.hpp"       // execution policies/traits

#ifdef AXOM_USE_RAJA
#include "RAJA/RAJA.hpp"
#endif

namespace axom
{
namespace mint
{
namespace internal
{

template < typename ExecPolicy, typename KernelType >
inline void for_all_faces( xargs::index, const Mesh* m,
                           KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  const IndexType numFaces = m->getNumberOfFaces();

#ifdef AXOM_USE_RAJA

  using exec_pol = typename policy_traits< ExecPolicy >::raja_exec_policy;
  RAJA::forall< exec_pol >( RAJA::RangeSegment( 0, numFaces ), kernel );

#else

  constexpr bool is_serial = std::is_same< ExecPolicy, policy::serial >::value;
  AXOM_STATIC_ASSERT( is_serial );

  for ( IndexType cellID = 0 ; cellID < numFaces ; ++cellID )
  {
    kernel( cellID );
  }

#endif
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_I_faces_2D( const Mesh* m, 
                                KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ERROR_IF( !m->isStructured(), "Mesh must be a 2D StructuredMesh." );
  SLIC_ERROR_IF( m->getDimension() != 2, "Mesh must be a 2D StructuredMesh." );

  const StructuredMesh* sm =
    static_cast< const StructuredMesh* >( m );

  const IndexType INodeResolution = sm->getNodeResolution( I_DIRECTION );
  const IndexType Ni = INodeResolution;
  const IndexType Nj = sm->getCellResolution( J_DIRECTION );

#ifdef AXOM_USE_RAJA

  RAJA::RangeSegment i_range(0,Ni);
  RAJA::RangeSegment j_range(0,Nj);

  using exec_pol = typename policy_traits< ExecPolicy >::raja_2d_exec;
  RAJA::kernel< exec_pol >( RAJA::make_tuple( i_range, j_range ),
    AXOM_LAMBDA( IndexType i, IndexType j )
    {
      const IndexType faceID = i + j * INodeResolution;
      kernel( faceID, i, j );
    }
  );

#else

  constexpr bool is_serial = std::is_same< ExecPolicy, policy::serial >::value;
  AXOM_STATIC_ASSERT( is_serial );

  for ( IndexType j = 0; j < Nj; ++j )
  {
    const IndexType offset = j * INodeResolution;
    for ( IndexType i = 0; i < Ni; ++i )
    {
      const IndexType faceID = i + offset;
      kernel( faceID, i, j );
    }
  }

#endif
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_J_faces_2D( const Mesh* m, 
                                KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ERROR_IF( !m->isStructured(), "Mesh must be a 2D StructuredMesh." );
  SLIC_ERROR_IF( m->getDimension() != 2, "Mesh must be a 2D StructuredMesh." );

  const StructuredMesh* sm =
    static_cast< const StructuredMesh* >( m );

  const IndexType ICellResolution = sm->getCellResolution( I_DIRECTION );
  const IndexType numIFaces       = sm->getTotalNumFaces( I_DIRECTION );
  const IndexType Ni              = ICellResolution;
  const IndexType Nj              = sm->getNodeResolution( J_DIRECTION );

#ifdef AXOM_USE_RAJA

  RAJA::RangeSegment i_range(0,Ni);
  RAJA::RangeSegment j_range(0,Nj);

  using exec_pol = typename policy_traits< ExecPolicy >::raja_2d_exec;
  RAJA::kernel< exec_pol >( RAJA::make_tuple( i_range, j_range ),
    AXOM_LAMBDA( IndexType i, IndexType j )
    {
      const IndexType faceID = numIFaces + i + j * ICellResolution;
      kernel( faceID, i, j );
    }
  );

#else

  constexpr bool is_serial = std::is_same< ExecPolicy, policy::serial >::value;
  AXOM_STATIC_ASSERT( is_serial );

  for ( IndexType j = 0; j < Nj; ++j )
  {
    const IndexType offset = numIFaces + j * ICellResolution;
    for ( IndexType i = 0; i < Ni; ++i )
    {
      const IndexType faceID = i + offset;
      kernel( faceID, i, j );
    }
  }

#endif
}


//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_facenodes_structured_2D( const Mesh* m,
                                              KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ERROR_IF( !m->isStructured(), "Mesh must be a 2D StructuredMesh." );
  SLIC_ERROR_IF( m->getDimension() != 2, "Mesh must be a 2D StructuredMesh." );

  const StructuredMesh* sm =
    static_cast< const StructuredMesh* >( m );

  const IndexType numIFaces = sm->getTotalNumFaces( I_DIRECTION );
  const IndexType* offsets  = sm->getCellNodeOffsetsArray();
  const IndexType cellNodeOffset3 = offsets[ 3 ];

  for_all_I_faces_2D< ExecPolicy >( m,
    AXOM_LAMBDA( IndexType faceID, IndexType AXOM_NOT_USED(i),
                 IndexType AXOM_NOT_USED(j) )
    {
      IndexType nodes[ 2 ];
      nodes[ 0 ] = faceID;
      nodes[ 1 ] = nodes[ 0 ] + cellNodeOffset3;
      kernel( faceID, nodes, 2 );
    }
  );

  for_all_J_faces_2D< ExecPolicy >( m,
    AXOM_LAMBDA( IndexType faceID, IndexType AXOM_NOT_USED(i), IndexType j )
    {
      const IndexType shiftedID = faceID - numIFaces;
      IndexType nodes[ 2 ];
      nodes[ 0 ] = shiftedID + j;
      nodes[ 1 ] = nodes[ 0 ] + 1;
      kernel( faceID, nodes, 2 );
    }
  );
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_I_faces_3D( const Mesh* m, 
                                KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ERROR_IF( !m->isStructured(), "Mesh must be a 3D StructuredMesh." );
  SLIC_ERROR_IF( m->getDimension() != 3, "Mesh must be a 3D StructuredMesh." );

  const StructuredMesh* sm =
    static_cast< const StructuredMesh* >( m );

  const IndexType INodeResolution = sm->getNodeResolution( I_DIRECTION );
  const IndexType numIFacesInKSlice =
                        INodeResolution * sm->getCellResolution( J_DIRECTION );
  const IndexType Ni = INodeResolution;
  const IndexType Nj = sm->getCellResolution( J_DIRECTION );
  const IndexType Nk = sm->getCellResolution( K_DIRECTION );

#ifdef AXOM_USE_RAJA

  RAJA::RangeSegment i_range(0,Ni);
  RAJA::RangeSegment j_range(0,Nj);
  RAJA::RangeSegment k_range(0,Nk);

  using exec_pol = typename policy_traits< ExecPolicy >::raja_3d_exec;
  RAJA::kernel< exec_pol >( RAJA::make_tuple( i_range, j_range, k_range ),
    AXOM_LAMBDA( IndexType i, IndexType j, IndexType k )
    {
      const IndexType faceID = i + j * INodeResolution + k * numIFacesInKSlice;
      kernel( faceID, i, j, k );
    }
  );

#else

  constexpr bool is_serial = std::is_same< ExecPolicy, policy::serial >::value;
  AXOM_STATIC_ASSERT( is_serial );

  for ( IndexType k = 0; k < Nk; ++k )
  {
    const IndexType k_offset = k * numIFacesInKSlice;
    for ( IndexType j = 0; j < Nj; ++j )
    {
      const IndexType offset = j * INodeResolution + k_offset;
      for ( IndexType i = 0; i < Ni; ++i )
      {
        const IndexType faceID = i + offset;
        kernel( faceID, i, j, k );
      }
    }
  }

#endif
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_J_faces_3D( const Mesh* m, 
                                KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ERROR_IF( !m->isStructured(), "Mesh must be a 3D StructuredMesh." );
  SLIC_ERROR_IF( m->getDimension() != 3, "Mesh must be a 3D StructuredMesh." );

  const StructuredMesh* sm =
    static_cast< const StructuredMesh* >( m );

  const IndexType numIFaces = sm->getTotalNumFaces( I_DIRECTION );
  const IndexType ICellResolution = sm->getCellResolution( I_DIRECTION );
  const IndexType numJFacesInKSlice =
                        ICellResolution * sm->getNodeResolution( J_DIRECTION );
  const IndexType Ni = ICellResolution;
  const IndexType Nj = sm->getNodeResolution( J_DIRECTION );
  const IndexType Nk = sm->getCellResolution( K_DIRECTION );

#ifdef AXOM_USE_RAJA

  RAJA::RangeSegment i_range(0,Ni);
  RAJA::RangeSegment j_range(0,Nj);
  RAJA::RangeSegment k_range(0,Nk);

  using exec_pol = typename policy_traits< ExecPolicy >::raja_3d_exec;
  RAJA::kernel< exec_pol >( RAJA::make_tuple( i_range, j_range, k_range ),
    AXOM_LAMBDA( IndexType i, IndexType j, IndexType k )
    {
      const IndexType jp     = j * ICellResolution;
      const IndexType kp     = k * numJFacesInKSlice;
      const IndexType faceID = numIFaces + i + jp + kp;
      kernel( faceID, i, j, k );
    }
  );

#else

  constexpr bool is_serial = std::is_same< ExecPolicy, policy::serial >::value;
  AXOM_STATIC_ASSERT( is_serial );

  for ( IndexType k = 0; k < Nk; ++k )
  {
    const IndexType k_offset = k * numJFacesInKSlice + numIFaces;
    for ( IndexType j = 0; j < Nj; ++j )
    {
      const IndexType offset = j * ICellResolution + k_offset;
      for ( IndexType i = 0; i < Ni; ++i )
      {
        const IndexType faceID = i + offset;
        kernel( faceID, i, j, k );
      }
    }
  }

#endif
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_K_faces_3D( const Mesh* m, 
                                KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ERROR_IF( !m->isStructured(), "Mesh must be a 3D StructuredMesh." );
  SLIC_ERROR_IF( m->getDimension() != 3, "Mesh must be a 3D StructuredMesh." );

  const StructuredMesh* sm =
    static_cast< const StructuredMesh* >( m );

  const IndexType numIJFaces = sm->getTotalNumFaces( I_DIRECTION ) +
                               sm->getTotalNumFaces( J_DIRECTION );
  const IndexType ICellResolution = sm->getCellResolution( I_DIRECTION );
  const IndexType cellKp = sm->cellKp();
  const IndexType Ni = ICellResolution;
  const IndexType Nj = sm->getCellResolution( J_DIRECTION );
  const IndexType Nk = sm->getNodeResolution( K_DIRECTION );

#ifdef AXOM_USE_RAJA

  RAJA::RangeSegment i_range(0,Ni);
  RAJA::RangeSegment j_range(0,Nj);
  RAJA::RangeSegment k_range(0,Nk);

  using exec_pol = typename policy_traits< ExecPolicy >::raja_3d_exec;
  RAJA::kernel< exec_pol >( RAJA::make_tuple( i_range, j_range, k_range ),
    AXOM_LAMBDA( IndexType i, IndexType j, IndexType k )
    {
      const IndexType jp = j * ICellResolution;
      const IndexType kp = k * cellKp;
      const IndexType faceID = numIJFaces + i + jp + kp;
      kernel( faceID, i, j, k );
    }
  );

#else

  constexpr bool is_serial = std::is_same< ExecPolicy, policy::serial >::value;
  AXOM_STATIC_ASSERT( is_serial );

  for ( IndexType k = 0; k < Nk; ++k )
  {
    const IndexType k_offset = k * cellKp + numIJFaces;
    for ( IndexType j = 0; j < Nj; ++j )
    {
      const IndexType offset = j * ICellResolution + k_offset;
      for ( IndexType i = 0; i < Ni; ++i )
      {
        const IndexType faceID = i + offset;
        kernel( faceID, i, j, k );
      }
    }
  }

#endif
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_facenodes_structured_3D( const Mesh* m,
                                              KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ERROR_IF( !m->isStructured(), "Mesh must be a 3D StructuredMesh." );
  SLIC_ERROR_IF( m->getDimension() != 3, "Mesh must be a 3D StructuredMesh." );

  const StructuredMesh* sm =
    static_cast< const StructuredMesh* >( m );

  const IndexType numIFaces = sm->getTotalNumFaces( I_DIRECTION );
  const IndexType numIJFaces = numIFaces + sm->getTotalNumFaces( J_DIRECTION );
  const IndexType INodeResolution = sm->getNodeResolution( I_DIRECTION );
  const IndexType JNodeResolution = sm->getNodeResolution( J_DIRECTION );
  const IndexType KFaceNodeStride = sm->getCellResolution( I_DIRECTION ) +
                                    sm->getCellResolution( J_DIRECTION ) + 1;

  const IndexType* offsets  = sm->getCellNodeOffsetsArray();
  const IndexType cellNodeOffset2 = offsets[ 2 ];
  const IndexType cellNodeOffset3 = offsets[ 3 ];
  const IndexType cellNodeOffset4 = offsets[ 4 ];
  const IndexType cellNodeOffset5 = offsets[ 5 ];
  const IndexType cellNodeOffset7 = offsets[ 7 ];

  for_all_I_faces_3D< ExecPolicy >( m,
    AXOM_LAMBDA( IndexType faceID, IndexType AXOM_NOT_USED(i),
                 IndexType AXOM_NOT_USED(j), IndexType k )
    {
      IndexType nodes[ 4 ];
      nodes[ 0 ] = faceID + k * INodeResolution;
      nodes[ 1 ] = nodes[ 0 ] + cellNodeOffset4;
      nodes[ 2 ] = nodes[ 0 ] + cellNodeOffset7;
      nodes[ 3 ] = nodes[ 0 ] + cellNodeOffset3;
      kernel( faceID, nodes, 4 );
    }
  );

  for_all_J_faces_3D< ExecPolicy >( m,
    AXOM_LAMBDA( IndexType faceID, IndexType AXOM_NOT_USED(i), IndexType j,
                 IndexType k )
    {
      const IndexType shiftedID = faceID - numIFaces;
      IndexType nodes[ 4 ];
      nodes[ 0 ] = shiftedID + j + k * JNodeResolution;
      nodes[ 1 ] = nodes[ 0 ] + 1;
      nodes[ 2 ] = nodes[ 0 ] + cellNodeOffset5;
      nodes[ 3 ] = nodes[ 0 ] + cellNodeOffset4;
      kernel( faceID, nodes, 4 );
    }
  );

  for_all_K_faces_3D< ExecPolicy >( m,
    AXOM_LAMBDA( IndexType faceID, IndexType AXOM_NOT_USED(i), IndexType j,
                 IndexType k )
    {
      const IndexType shiftedID = faceID - numIJFaces;
      IndexType nodes[ 4 ];
      nodes[ 0 ] = shiftedID + j + k * KFaceNodeStride;
      nodes[ 1 ] = nodes[ 0 ] + 1;
      nodes[ 2 ] = nodes[ 0 ] + cellNodeOffset2;
      nodes[ 3 ] = nodes[ 0 ] + cellNodeOffset3;
      kernel( faceID, nodes, 4 );
    }
  );
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_facescells_structured_2D( const Mesh* m,
                                              KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ERROR_IF( !m->isStructured(), "Mesh must be a 2D StructuredMesh." );
  SLIC_ERROR_IF( m->getDimension() != 2, "Mesh must be a 2D StructuredMesh." );

  const StructuredMesh* sm = 
    static_cast< const StructuredMesh* >( m );

  const IndexType ICellResolution = sm->getCellResolution( I_DIRECTION );
  const IndexType JCellResolution = sm->getCellResolution( J_DIRECTION );

  const IndexType cellJp = sm->cellJp();

  for_all_I_faces_2D< ExecPolicy >( m,
    AXOM_LAMBDA( IndexType faceID, IndexType i, IndexType j )
    {
      IndexType cellIDTwo = i + j * cellJp;
      IndexType cellIDOne = cellIDTwo - 1;
      if ( i == 0 )
      {
        cellIDOne = cellIDTwo;
        cellIDTwo = -1;
      }
      else if ( i == ICellResolution )
      {
        cellIDTwo = -1;
      }

      kernel( faceID, cellIDOne, cellIDTwo );
    }
  );

  for_all_J_faces_2D< ExecPolicy >( m,
    AXOM_LAMBDA( IndexType faceID, IndexType i, IndexType j )
    {
      IndexType cellIDTwo = i + j * cellJp;
      IndexType cellIDOne = cellIDTwo - cellJp;
      if ( j == 0 )
      {
        cellIDOne = cellIDTwo;
        cellIDTwo = -1;
      }
      else if ( j == JCellResolution )
      {
        cellIDTwo = -1;
      }

      kernel( faceID, cellIDOne, cellIDTwo );
    }
  );
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_facescells_structured_3D( const Mesh* m,
                                              KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ERROR_IF( !m->isStructured(), "Mesh must be a 3D StructuredMesh." );
  SLIC_ERROR_IF( m->getDimension() != 3, "Mesh must be a 3D StructuredMesh." );

  const StructuredMesh* sm = 
    static_cast< const StructuredMesh* >( m );

  const IndexType ICellResolution = sm->getCellResolution( I_DIRECTION );
  const IndexType JCellResolution = sm->getCellResolution( J_DIRECTION );
  const IndexType KCellResolution = sm->getCellResolution( K_DIRECTION );

  const IndexType cellJp = sm->cellJp();
  const IndexType cellKp = sm->cellKp();

  for_all_I_faces_3D< ExecPolicy >( m,
    AXOM_LAMBDA( IndexType faceID, IndexType i, IndexType j, IndexType k )
    {
      IndexType cellIDTwo = i + j * cellJp + k * cellKp;
      IndexType cellIDOne = cellIDTwo - 1;
      if ( i == 0 )
      {
        cellIDOne = cellIDTwo;
        cellIDTwo = -1;
      }
      else if ( i == ICellResolution )
      {
        cellIDTwo = -1;
      }

      kernel( faceID, cellIDOne, cellIDTwo );
    }
  );

  for_all_J_faces_3D< ExecPolicy >( m,
    AXOM_LAMBDA( IndexType faceID, IndexType i, IndexType j, IndexType k )
    {
      IndexType cellIDTwo = i + j * cellJp + k * cellKp;
      IndexType cellIDOne = cellIDTwo - cellJp;
      if ( j == 0 )
      {
        cellIDOne = cellIDTwo;
        cellIDTwo = -1;
      }
      else if ( j == JCellResolution )
      {
        cellIDTwo = -1;
      }

      kernel( faceID, cellIDOne, cellIDTwo );
    }
  );

  for_all_K_faces_3D< ExecPolicy >( m,
    AXOM_LAMBDA( IndexType faceID, IndexType i, IndexType j, IndexType k )
    {
      IndexType cellIDTwo = i + j * cellJp + k * cellKp;
      IndexType cellIDOne = cellIDTwo - cellKp;
      if ( k == 0 )
      {
        cellIDOne = cellIDTwo;
        cellIDTwo = -1;
      }
      else if ( k == KCellResolution )
      {
        cellIDTwo = -1;
      }

      kernel( faceID, cellIDOne, cellIDTwo );
    }
  );
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_facenodes_unstructured_single( const Mesh* m,
                                                   KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ERROR_IF( m->getNumberOfFaces() <= 0,
                 "No faces in the mesh, perhaps you meant to call " <<
                 "UnstructuredMesh::initializeFaceConnectivity first." );
  SLIC_ERROR_IF( m->isStructured(),
                 "Mesh must be an UnstructuredMesh<SINGLE_SHAPE>." );
  SLIC_ERROR_IF( m->hasMixedCellTypes(),
                 "Mesh must be an UnstructuredMesh<SINGLE_SHAPE>." );

  using UnstructuredMeshType = UnstructuredMesh< SINGLE_SHAPE >;

  const UnstructuredMeshType* um =
                                static_cast< const UnstructuredMeshType* >( m );

  const IndexType* faces_to_nodes = um->getFaceNodesArray();
  const IndexType num_nodes = um->getNumberOfFaceNodes();

  for_all_faces< ExecPolicy, xargs::index >( m,
    AXOM_LAMBDA( IndexType faceID )
    {
      kernel( faceID, faces_to_nodes + faceID * num_nodes, num_nodes );
    }
  );
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_facenodes_unstructured_mixed( const Mesh* m,
                                                  KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ERROR_IF( m->getNumberOfFaces() <= 0,
                 "No faces in the mesh, perhaps you meant to call " <<
                 "UnstructuredMesh::initializeFaceConnectivity first." );
  SLIC_ERROR_IF( m->isStructured(),
                 "Mesh must be an UnstructuredMesh<SINGLE_SHAPE>." );
  SLIC_ERROR_IF( !m->hasMixedCellTypes(),
                 "Mesh must be an UnstructuredMesh<MIXED_SHAPE>." );

  using UnstructuredMeshType = UnstructuredMesh< MIXED_SHAPE >;

  const UnstructuredMeshType* um =
                                static_cast< const UnstructuredMeshType* >( m );

  const IndexType* faces_to_nodes = um->getFaceNodesArray();
  const IndexType* offsets        = um->getFaceNodesOffsetsArray();

  for_all_faces< ExecPolicy, xargs::index >( m,
    AXOM_LAMBDA( IndexType faceID )
    {
      const IndexType num_nodes = offsets[ faceID + 1 ] - offsets[ faceID ];
      kernel( faceID, faces_to_nodes + offsets[ faceID ], num_nodes );
    }
  );
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, Topology TOPO, typename KernelType >
inline void for_all_facecells_unstructured( const Mesh* m,
                                            KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ERROR_IF( m->getNumberOfFaces() <= 0,
                 "No faces in the mesh, perhaps you meant to call " <<
                 "UnstructuredMesh::initializeFaceConnectivity first." );
  SLIC_ERROR_IF( m->isStructured(), "Mesh must be an UnstructuredMesh." );

  using UnstructuredMeshType = UnstructuredMesh< TOPO >;

  const UnstructuredMeshType* um =
                                static_cast< const UnstructuredMeshType* >( m );

  const IndexType* faces_to_cells = um->getFaceCellsArray();

  for_all_faces< ExecPolicy, xargs::index >( m,
    AXOM_LAMBDA( IndexType faceID )
    {
      const IndexType offset = 2 * faceID;
      kernel( faceID, faces_to_cells[ offset ], faces_to_cells[ offset + 1 ] );
    }
  );
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_facecoords_uniform( const Mesh* m, KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ASSERT( m->getMeshType() == STRUCTURED_UNIFORM_MESH );
  SLIC_ASSERT( m->getDimension() > 1 && m->getDimension() <= 3 );

  constexpr bool NO_COPY = true;

  const UniformMesh* um  = static_cast< const UniformMesh* >( m );
  const int dimension    = um->getDimension();
  const double * x0      = um->getOrigin( );
  const double * h       = um->getSpacing( );
  const IndexType nodeJp = um->nodeJp();
  const IndexType nodeKp = um->nodeKp();

  if ( dimension == 2 )
  {
    for_all_I_faces_2D< ExecPolicy >( um,
      AXOM_LAMBDA( IndexType faceID, IndexType i, IndexType j )
      {
        const IndexType n0 = i + j * nodeJp;
        const IndexType nodeIDs[2] = { n0, n0 + nodeJp };

        double coords[4] = { x0[0] + i * h[0], x0[1] +  j      * h[1],
                             x0[0] + i * h[0], x0[1] + (j + 1) * h[1] };
        
        numerics::Matrix<double> coordsMatrix( dimension, 2, coords, NO_COPY );
        kernel( faceID, coordsMatrix, nodeIDs );
      }
    );

    for_all_J_faces_2D< ExecPolicy >( um,
      AXOM_LAMBDA( IndexType faceID, IndexType i, IndexType j )
      {
        const IndexType n0 = i + j * nodeJp;
        const IndexType nodeIDs[2] = { n0, n0 + 1 };

        double coords[4] = { x0[0] +  i      * h[0], x0[1] + j * h[1],
                             x0[0] + (i + 1) * h[0], x0[1] + j * h[1] };
        
        numerics::Matrix<double> coordsMatrix( dimension, 2, coords, NO_COPY );
        kernel( faceID, coordsMatrix, nodeIDs );
      }
    );
  }
  else
  {
    for_all_I_faces_3D< ExecPolicy >( um,
      AXOM_LAMBDA( IndexType faceID, IndexType i, IndexType j, IndexType k )
      {
        const IndexType n0 = i + j * nodeJp + k * nodeKp;
        const IndexType nodeIDs[4] = { n0,
                                       n0 + nodeKp,
                                       n0 + nodeJp + nodeKp,
                                       n0 + nodeJp };

        double coords[12] = { 
          x0[0] + i * h[0], x0[1] +  j      * h[1], x0[2] +  k      * h[2],
          x0[0] + i * h[0], x0[1] +  j      * h[1], x0[2] + (k + 1) * h[2],
          x0[0] + i * h[0], x0[1] + (j + 1) * h[1], x0[2] + (k + 1) * h[2],
          x0[0] + i * h[0], x0[1] + (j + 1) * h[1], x0[2] +  k      * h[2] };
        
        numerics::Matrix<double> coordsMatrix( dimension, 4, coords, NO_COPY );
        kernel( faceID, coordsMatrix, nodeIDs );
      }
    );

    for_all_J_faces_3D< ExecPolicy >( um,
      AXOM_LAMBDA( IndexType faceID, IndexType i, IndexType j, IndexType k )
      {
        const IndexType n0 = i + j * nodeJp + k * nodeKp;
        const IndexType nodeIDs[4] = { n0,
                                       n0 + 1,
                                       n0 + 1 + nodeKp,
                                       n0 + nodeKp };

        double coords[12] = { 
          x0[0] +  i      * h[0], x0[1] + j * h[1], x0[2] +  k      * h[2],
          x0[0] + (i + 1) * h[0], x0[1] + j * h[1], x0[2] +  k      * h[2],
          x0[0] + (i + 1) * h[0], x0[1] + j * h[1], x0[2] + (k + 1) * h[2],
          x0[0] +  i      * h[0], x0[1] + j * h[1], x0[2] + (k + 1) * h[2] };
        
        numerics::Matrix<double> coordsMatrix( dimension, 4, coords, NO_COPY );
        kernel( faceID, coordsMatrix, nodeIDs );
      }
    );

    for_all_K_faces_3D< ExecPolicy >( um,
      AXOM_LAMBDA( IndexType faceID, IndexType i, IndexType j, IndexType k )
      {
        const IndexType n0 = i + j * nodeJp + k * nodeKp;
        const IndexType nodeIDs[4] = { n0,
                                       n0 + 1,
                                       n0 + 1 + nodeJp,
                                       n0 + nodeJp };

        double coords[12] = { 
          x0[0] +  i      * h[0], x0[1] +  j      * h[1], x0[2] + k * h[2],
          x0[0] + (i + 1) * h[0], x0[1] +  j      * h[1], x0[2] + k * h[2],
          x0[0] + (i + 1) * h[0], x0[1] + (j + 1) * h[1], x0[2] + k * h[2],
          x0[0] +  i      * h[0], x0[1] + (j + 1) * h[1], x0[2] + k * h[2] };
        
        numerics::Matrix<double> coordsMatrix( dimension, 4, coords, NO_COPY );
        kernel( faceID, coordsMatrix, nodeIDs );
      }
    );
  }
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_facecoords_rectilinear( const Mesh* m, KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ASSERT( m->getMeshType() == STRUCTURED_RECTILINEAR_MESH );
  SLIC_ASSERT( m->getDimension() > 1 && m->getDimension() <= 3 );

  constexpr bool NO_COPY = true;

  const RectilinearMesh* rm = static_cast< const RectilinearMesh* >( m );
  const int dimension       = rm->getDimension();
  const IndexType nodeJp    = rm->nodeJp();
  const IndexType nodeKp    = rm->nodeKp();
  const double * x          = rm->getCoordinateArray( X_COORDINATE );
  const double * y          = rm->getCoordinateArray( Y_COORDINATE );

  if ( dimension == 2 )
  {
    for_all_I_faces_2D< ExecPolicy >( rm,
      AXOM_LAMBDA( IndexType faceID, IndexType i, IndexType j )
      {
        const IndexType n0 = i + j * nodeJp;
        const IndexType nodeIDs[2] = { n0, n0 + nodeJp };

        double coords[4] = { x[ i ], y[ j ],
                             x[ i ], y[ j + 1 ] };
        
        numerics::Matrix<double> coordsMatrix( dimension, 2, coords, NO_COPY );
        kernel( faceID, coordsMatrix, nodeIDs );
      }
    );

    for_all_J_faces_2D< ExecPolicy >( rm,
      AXOM_LAMBDA( IndexType faceID, IndexType i, IndexType j )
      {
        const IndexType n0 = i + j * nodeJp;
        const IndexType nodeIDs[2] = { n0, n0 + 1 };

        double coords[4] = { x[ i ]   , y[ j ],
                             x[ i + 1], y[ j ] };
        
        numerics::Matrix<double> coordsMatrix( dimension, 2, coords, NO_COPY );
        kernel( faceID, coordsMatrix, nodeIDs );
      }
    );
  }
  else
  {
    const double * z = rm->getCoordinateArray( Z_COORDINATE );
    for_all_I_faces_3D< ExecPolicy >( rm,
      AXOM_LAMBDA( IndexType faceID, IndexType i, IndexType j, IndexType k )
      {
        const IndexType n0 = i + j * nodeJp + k * nodeKp;
        const IndexType nodeIDs[4] = { n0,
                                       n0 + nodeKp,
                                       n0 + nodeJp + nodeKp,
                                       n0 + nodeJp };

        double coords[12] = { x[ i ], y[ j     ], z[ k     ],
                              x[ i ], y[ j     ], z[ k + 1 ],
                              x[ i ], y[ j + 1 ], z[ k + 1 ],
                              x[ i ], y[ j + 1 ], z[ k     ] };
        
        numerics::Matrix<double> coordsMatrix( dimension, 4, coords, NO_COPY );
        kernel( faceID, coordsMatrix, nodeIDs );
      }
    );

    for_all_J_faces_3D< ExecPolicy >( rm,
      AXOM_LAMBDA( IndexType faceID, IndexType i, IndexType j, IndexType k )
      {
        const IndexType n0 = i + j * nodeJp + k * nodeKp;
        const IndexType nodeIDs[4] = { n0,
                                       n0 + 1,
                                       n0 + 1 + nodeKp,
                                       n0 + nodeKp };

        double coords[12] = { x[ i     ], y[ j ], z[ k     ],
                              x[ i + 1 ], y[ j ], z[ k     ],
                              x[ i + 1 ], y[ j ], z[ k + 1 ],
                              x[ i     ], y[ j ], z[ k + 1 ] };
        
        numerics::Matrix<double> coordsMatrix( dimension, 4, coords, NO_COPY );
        kernel( faceID, coordsMatrix, nodeIDs );
      }
    );

    for_all_K_faces_3D< ExecPolicy >( rm,
      AXOM_LAMBDA( IndexType faceID, IndexType i, IndexType j, IndexType k )
      {
        const IndexType n0 = i + j * nodeJp + k * nodeKp;
        const IndexType nodeIDs[4] = { n0,
                                       n0 + 1,
                                       n0 + 1 + nodeJp,
                                       n0 + nodeJp };

        double coords[12] = { x[ i     ], y[ j     ], z[ k ],
                              x[ i + 1 ], y[ j     ], z[ k ],
                              x[ i + 1 ], y[ j + 1 ], z[ k ],
                              x[ i     ], y[ j + 1 ], z[ k ] };
        
        numerics::Matrix<double> coordsMatrix( dimension, 4, coords, NO_COPY );
        kernel( faceID, coordsMatrix, nodeIDs );
      }
    );
  }
}

//------------------------------------------------------------------------------
struct for_all_face_nodes_functor
{
  template < typename ExecPolicy, typename KernelType >
  inline void operator()( ExecPolicy AXOM_NOT_USED(policy), const Mesh* m,
                          KernelType&& kernel ) const
  {
    for_all_faces< ExecPolicy >( xargs::nodeids(), m,
                                 std::forward< KernelType >( kernel ) );
  }
};

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_facecoords_curvilinear( const Mesh* m, KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ASSERT( m->getMeshType() == STRUCTURED_CURVILINEAR_MESH );
  SLIC_ASSERT( m->getDimension() > 1 && m->getDimension() <= 3 );

  const int dimension = m->getDimension();
  if ( dimension == 2 )
  {
    for_all_coords< ExecPolicy, 2, 2 >( for_all_face_nodes_functor(), m,
                                        std::forward< KernelType >( kernel ) );
  }
  else
  {
    for_all_coords< ExecPolicy, 3, 4 >( for_all_face_nodes_functor(), m,
                                        std::forward< KernelType >( kernel ) );
  }
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_facecoords_unstructured( const Mesh* m,
                                             KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ASSERT( m->getMeshType() == UNSTRUCTURED_MESH );
  SLIC_ASSERT( m->getDimension() > 1 && m->getDimension() <= 3 );

  constexpr bool NO_COPY = true;

  const int dimension = m->getDimension();
  const double * x = m->getCoordinateArray( X_COORDINATE );
  const double * y = m->getCoordinateArray( Y_COORDINATE );

  if ( dimension == 2 )
  {
    for_all_faces< ExecPolicy >( xargs::nodeids(), m, 
      AXOM_LAMBDA( IndexType faceID, const IndexType * nodeIDs,
                   IndexType numNodes )
      {
        double coords[ 2 * MAX_FACE_NODES ];
        for ( int i = 0; i < numNodes; ++ i )
        {
          const IndexType nodeID = nodeIDs[ i ];
          coords[ 2 * i     ] = x[ nodeID ];
          coords[ 2 * i + 1 ] = y[ nodeID ];
        }
        
        numerics::Matrix<double> coordsMatrix( dimension, 2, coords, NO_COPY );
        kernel( faceID, coordsMatrix, nodeIDs );
      }
    );
  }
  else
  {
    const double * z = m->getCoordinateArray( Z_COORDINATE );
    for_all_faces< ExecPolicy >( xargs::nodeids(), m, 
      AXOM_LAMBDA( IndexType faceID, const IndexType * nodeIDs,
                   IndexType numNodes )
      {
        double coords[ 3 * MAX_FACE_NODES ];
        for ( int i = 0; i < numNodes; ++ i )
        {
          const IndexType nodeID = nodeIDs[ i ];
          coords[ 3 * i     ] = x[ nodeID ];
          coords[ 3 * i + 1 ] = y[ nodeID ];
          coords[ 3 * i + 2 ] = z[ nodeID ];
        }
        
        numerics::Matrix<double> coordsMatrix( dimension, 4, coords, NO_COPY );
        kernel( faceID, coordsMatrix, nodeIDs );
      }
    );
  }
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_faces( xargs::nodeids, const Mesh* m,
                           KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ASSERT( m->getDimension() > 1 && m->getDimension() <= 3 );

  if ( m->isStructured() )
  {
    if ( m->getDimension() == 2 )
    {
      for_all_facenodes_structured_2D< ExecPolicy >( m,
        std::forward< KernelType >( kernel ) );
    }
    else
    {
      for_all_facenodes_structured_3D< ExecPolicy >( m,
        std::forward< KernelType >( kernel ) );
    }
  }
  else if ( m->hasMixedCellTypes() )
  {
    for_all_facenodes_unstructured_mixed< ExecPolicy >(
      m, std::forward< KernelType >( kernel ) );
  }
  else
  {
    for_all_facenodes_unstructured_single< ExecPolicy >(
      m, std::forward< KernelType >( kernel ) );
  }
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_faces( xargs::cellids, const Mesh* m,
                           KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ASSERT( m->getDimension() > 1 && m->getDimension() <= 3 );

  if ( m->isStructured() )
  {
    if ( m->getDimension() == 2 )
    {
      for_all_facescells_structured_2D< ExecPolicy >( m,
        std::forward< KernelType >( kernel ) );
    }
    else
    {
      for_all_facescells_structured_3D< ExecPolicy >( m,
        std::forward< KernelType >( kernel ) );
    }
  }
  else if ( m->hasMixedCellTypes() )
  {
    for_all_facecells_unstructured< ExecPolicy, MIXED_SHAPE >(
      m, std::forward< KernelType >( kernel ) );
  }
  else
  {
    for_all_facecells_unstructured< ExecPolicy, SINGLE_SHAPE >(
      m, std::forward< KernelType >( kernel ) );
  }
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_faces( xargs::coords, const Mesh* m, KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ASSERT( m->getDimension() > 1 && m->getDimension() <= 3 );

  if ( m->getMeshType() == STRUCTURED_UNIFORM_MESH )
  {
    for_all_facecoords_uniform< ExecPolicy >(
      m, std::forward< KernelType >( kernel ) );
  }
  else if ( m->getMeshType() == STRUCTURED_RECTILINEAR_MESH )
  {
    for_all_facecoords_rectilinear< ExecPolicy >( 
      m, std::forward< KernelType >( kernel ) );
  }
  else if ( m->getMeshType() == STRUCTURED_CURVILINEAR_MESH )
  {
    for_all_facecoords_curvilinear< ExecPolicy >( 
      m, std::forward< KernelType >( kernel ) );
  }
  else if ( m->getMeshType() == UNSTRUCTURED_MESH )
  {
    for_all_facecoords_unstructured< ExecPolicy >(
      m, std::forward< KernelType >( kernel ) );
  }
  else
  {
    SLIC_ERROR( "Unknown mesh type." );
  }
}

} /* namespace internal */
} /* namespace mint     */
} /* namespace axom     */

#endif /* MINT_FOR_ALL_FACES_HPP_ */
