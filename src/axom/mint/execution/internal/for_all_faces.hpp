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
#include "axom/mint/mesh/Mesh.hpp"              // for mint::Mesh
#include "axom/mint/mesh/StructuredMesh.hpp"    // for mint::StructuredMesh
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
inline void for_all_faces( xargs::index, const mint::Mesh* m,
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
inline void for_all_I_faces_2D( const mint::Mesh* m, 
                                KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ASSERT( m->isStructured() );
  SLIC_ASSERT( m->getDimension() == 2 );

  const mint::StructuredMesh* sm =
    static_cast< const mint::StructuredMesh* >( m );

  const IndexType INodeResolution = sm->getNodeResolution( I_DIRECTION );
  const IndexType i_range = INodeResolution;
  const IndexType j_range = sm->getCellResolution( J_DIRECTION );

#ifdef AXOM_USE_RAJA
  
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

  for ( IndexType j = 0; j < j_range; ++j )
  {
    const IndexType offset = j * INodeResolution;
    for ( IndexType i = 0; i < i_range; ++i )
    {
      const IndexType faceID = i + offset;
      kernel( faceID, i, j );
    }
  }

#endif
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_J_faces_2D( const mint::Mesh* m, 
                                KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ASSERT( m->isStructured() );
  SLIC_ASSERT( m->getDimension() == 2 );

  const mint::StructuredMesh* sm =
    static_cast< const mint::StructuredMesh* >( m );

  const IndexType ICellResolution = sm->getCellResolution( I_DIRECTION );
  const IndexType numIFaces = sm->getTotalNumFaces( I_DIRECTION );
  const IndexType i_range = ICellResolution;
  const IndexType j_range = sm->getNodeResolution( J_DIRECTION );

#ifdef AXOM_USE_RAJA
  
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

  for ( IndexType j = 0; j < j_range; ++j )
  {
    const IndexType offset = numIFaces + j * ICellResolution;
    for ( IndexType i = 0; i < i_range; ++i )
    {
      const IndexType faceID = i + offset;
      kernel( faceID, i, j );
    }
  }

#endif
}


//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_facenodes_structured_2D( const mint::Mesh* m,
                                              KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ASSERT( m->isStructured() );
  SLIC_ASSERT( m->getDimension() == 2 );

  const mint::StructuredMesh* sm =
    static_cast< const mint::StructuredMesh* >( m );

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
inline void for_all_I_faces_3D( const mint::Mesh* m, 
                                KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ASSERT( m->isStructured() );
  SLIC_ASSERT( m->getDimension() == 3 );

  const mint::StructuredMesh* sm =
    static_cast< const mint::StructuredMesh* >( m );

  const IndexType INodeResolution = sm->getNodeResolution( I_DIRECTION );
  const IndexType numIFacesInKSlice = 
                        INodeResolution * sm->getCellResolution( J_DIRECTION );
  const IndexType i_range = INodeResolution;
  const IndexType j_range = sm->getCellResolution( J_DIRECTION );
  const IndexType k_range = sm->getCellResolution( K_DIRECTION );

#ifdef AXOM_USE_RAJA
  
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

  for ( IndexType k = 0; k < k_range; ++k )
  {
    const IndexType k_offset = k * numIFacesInKSlice;
    for ( IndexType j = 0; j < j_range; ++j )
    {
      const IndexType offset = j * INodeResolution + k_offset;
      for ( IndexType i = 0; i < i_range; ++i )
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
inline void for_all_J_faces_3D( const mint::Mesh* m, 
                                KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ASSERT( m->isStructured() );
  SLIC_ASSERT( m->getDimension() == 3 );

  const mint::StructuredMesh* sm =
    static_cast< const mint::StructuredMesh* >( m );

  const IndexType numIFaces = sm->getTotalNumFaces( I_DIRECTION );
  const IndexType ICellResolution = sm->getCellResolution( I_DIRECTION );
  const IndexType numJFacesInKSlice = 
                        ICellResolution * sm->getNodeResolution( J_DIRECTION );
  const IndexType i_range = ICellResolution;
  const IndexType j_range = sm->getNodeResolution( J_DIRECTION );
  const IndexType k_range = sm->getCellResolution( K_DIRECTION );

#ifdef AXOM_USE_RAJA
  
  using exec_pol = typename policy_traits< ExecPolicy >::raja_3d_exec;
  RAJA::kernel< exec_pol >( RAJA::make_tuple( i_range, j_range, k_range ),
    AXOM_LAMBDA( IndexType i, IndexType j, IndexType k )
    {
      const IndexType faceID = numIFaces + i + j * ICellResolution + k * numJFacesInKSlice;
      kernel( faceID, i, j, k ); 
    }
  );

#else

  constexpr bool is_serial = std::is_same< ExecPolicy, policy::serial >::value;
  AXOM_STATIC_ASSERT( is_serial );

  for ( IndexType k = 0; k < k_range; ++k )
  {
    const IndexType k_offset = k * numJFacesInKSlice + numIFaces;
    for ( IndexType j = 0; j < j_range; ++j )
    {
      const IndexType offset = j * ICellResolution + k_offset;
      for ( IndexType i = 0; i < i_range; ++i )
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
inline void for_all_K_faces_3D( const mint::Mesh* m, 
                                KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ASSERT( m->isStructured() );
  SLIC_ASSERT( m->getDimension() == 3 );

  const mint::StructuredMesh* sm =
    static_cast< const mint::StructuredMesh* >( m );

  const IndexType numIJFaces = sm->getTotalNumFaces( I_DIRECTION ) + 
                               sm->getTotalNumFaces( J_DIRECTION );
  const IndexType ICellResolution = sm->getCellResolution( I_DIRECTION );
  const IndexType cellKp = sm->cellKp();
  const IndexType i_range = ICellResolution;
  const IndexType j_range = sm->getCellResolution( J_DIRECTION );
  const IndexType k_range = sm->getNodeResolution( K_DIRECTION );

#ifdef AXOM_USE_RAJA
  
  using exec_pol = typename policy_traits< ExecPolicy >::raja_3d_exec;
  RAJA::kernel< exec_pol >( RAJA::make_tuple( i_range, j_range, k_range ),
    AXOM_LAMBDA( IndexType i, IndexType j, IndexType k )
    {
      const IndexType faceID = numIJFaces + i + j * ICellResolution + k * cellKp;
      kernel( faceID, i, j, k ); 
    }
  );

#else

  constexpr bool is_serial = std::is_same< ExecPolicy, policy::serial >::value;
  AXOM_STATIC_ASSERT( is_serial );

  for ( IndexType k = 0; k < k_range; ++k )
  {
    const IndexType k_offset = k * cellKp + numIJFaces;
    for ( IndexType j = 0; j < j_range; ++j )
    {
      const IndexType offset = j * ICellResolution + k_offset;
      for ( IndexType i = 0; i < i_range; ++i )
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
inline void for_all_facenodes_structured_3D( const mint::Mesh* m,
                                              KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ASSERT( m->isStructured() );
  SLIC_ASSERT( m->getDimension() == 3 );

  const mint::StructuredMesh* sm =
    static_cast< const mint::StructuredMesh* >( m );

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
inline void for_all_facescells_structured_2D( const mint::Mesh* m,
                                              KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ASSERT( m->isStructured() );
  SLIC_ASSERT( m->getDimension() == 2 );

  const mint::StructuredMesh* sm = 
    static_cast< const mint::StructuredMesh* >( m );

  const IndexType ICellResolution = sm->getCellResolution( I_DIRECTION );
  const IndexType JCellResolution = sm->getCellResolution( J_DIRECTION );
  
  const IndexType cellJp = sm->cellJp();

  for_all_I_faces_2D< ExecPolicy >( m,
    AXOM_LAMBDA( IndexType faceID, IndexType i, IndexType j )
    {
      IndexType cellIDOne = -1;
      IndexType cellIDTwo = -1;

      const IndexType cellID = i + j * cellJp;
      if ( i != 0 ) cellIDOne = cellID - 1;
      if ( i != ICellResolution ) cellIDTwo = cellID; 
      
      kernel( faceID, cellIDOne, cellIDTwo );
    }
  );

  for_all_J_faces_2D< ExecPolicy >( m,
    AXOM_LAMBDA( IndexType faceID, IndexType i, IndexType j )
    {
      IndexType cellIDOne = -1;
      IndexType cellIDTwo = -1;

      const IndexType cellID = i + j * cellJp;
            if ( j != 0 ) cellIDOne = cellID - cellJp;
      if ( j != JCellResolution ) cellIDTwo = cellID; 
      
      kernel( faceID, cellIDOne, cellIDTwo );
    }
  );
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_facescells_structured_3D( const mint::Mesh* m,
                                              KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ASSERT( m->isStructured() );
  SLIC_ASSERT( m->getDimension() == 3 );

  const mint::StructuredMesh* sm = 
    static_cast< const mint::StructuredMesh* >( m );

  const IndexType ICellResolution = sm->getCellResolution( I_DIRECTION );
  const IndexType JCellResolution = sm->getCellResolution( J_DIRECTION );
  const IndexType KCellResolution = sm->getCellResolution( K_DIRECTION );

  const IndexType cellJp = sm->cellJp();
  const IndexType cellKp = sm->cellKp();

  for_all_I_faces_3D< ExecPolicy >( m,
    AXOM_LAMBDA( IndexType faceID, IndexType i, IndexType j, IndexType k )
    {
      IndexType cellIDOne = -1;
      IndexType cellIDTwo = -1;

      const IndexType cellID = i + j * cellJp + k * cellKp;
      if ( i != 0 ) cellIDOne = cellID - 1;
      if ( i != ICellResolution ) cellIDTwo = cellID; 
      
      kernel( faceID, cellIDOne, cellIDTwo );
    }
  );

  for_all_J_faces_3D< ExecPolicy >( m,
    AXOM_LAMBDA( IndexType faceID, IndexType i, IndexType j, IndexType k )
    {
      IndexType cellIDOne = -1;
      IndexType cellIDTwo = -1;

      const IndexType cellID = i + j * cellJp + k * cellKp;
      if ( j != 0 ) cellIDOne = cellID - cellJp;
      if ( j != JCellResolution ) cellIDTwo = cellID; 
      
      kernel( faceID, cellIDOne, cellIDTwo );
    }
  );

  for_all_K_faces_3D< ExecPolicy >( m,
    AXOM_LAMBDA( IndexType faceID, IndexType i, IndexType j, IndexType k )
    {
      IndexType cellIDOne = -1;
      IndexType cellIDTwo = -1;

      const IndexType cellID = i + j * cellJp + k * cellKp;
      if ( k != 0 ) cellIDOne = cellID - cellKp;
      if ( k != KCellResolution ) cellIDTwo = cellID; 
      
      kernel( faceID, cellIDOne, cellIDTwo );
    }
  );
}


//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_faces( xargs::nodeids, const mint::Mesh* m,
                           KernelType&& kernel )
{
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
  else
  {
    SLIC_ERROR("For all faces only supported for structured meshes.");
  }
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_faces( xargs::cellids, const mint::Mesh* m,
                           KernelType&& kernel )
{
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
  else
  {
    SLIC_ERROR("For all faces only supported for structured meshes.");
  }
}

} /* namespace internal */
} /* namespace mint     */
} /* namespace axom     */

#endif /* MINT_FOR_ALL_FACES_HPP_ */
