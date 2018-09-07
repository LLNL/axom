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

#ifndef MINT_FOR_ALL_CELLS_HPP_
#define MINT_FOR_ALL_CELLS_HPP_

// mint includes
#include "axom/mint/execution/xargs.hpp"        // for xargs

#include "axom/mint/config.hpp"                 // for compile-time definitions
#include "axom/mint/mesh/Mesh.hpp"              // for mint::Mesh
#include "axom/mint/mesh/StructuredMesh.hpp"    // for mint::StructuredMesh
#include "axom/mint/mesh/UnstructuredMesh.hpp"  // for mint::UnstructuredMesh
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
inline void for_all_cells( xargs::index, const mint::Mesh* m,
                           KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  const IndexType numCells = m->getNumberOfCells();

#ifdef AXOM_USE_RAJA

  using exec_pol = typename policy_traits< ExecPolicy >::raja_exec_policy;
  RAJA::forall< exec_pol >( RAJA::RangeSegment(0,numCells), kernel );

#else

  constexpr bool is_serial = std::is_same< ExecPolicy, policy::serial >::value;
  AXOM_STATIC_ASSERT( is_serial );

  for ( IndexType cellIdx=0; cellIdx < numCells; ++cellIdx )
  {
    kernel( cellIdx );
  }

#endif

}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_cells( xargs::ij, const mint::Mesh* m,
                           KernelType&& kernel )
{
  // run-time checks
  SLIC_ASSERT( m != nullptr );
  SLIC_ERROR_IF( !m->isStructured(),
                 "xargs::ij is only valid on Structured meshes!" );
  SLIC_ERROR_IF( m->getDimension() != 2,
                 "xargs::ij is only valid for 2-D structured meshes!" );

  const mint::StructuredMesh* sm =
      static_cast< const mint::StructuredMesh* >( m );

  const IndexType jp = sm->cellJp();
  const IndexType Ni = sm->getCellResolution( I_DIRECTION );
  const IndexType Nj = sm->getCellResolution( J_DIRECTION );

#ifdef AXOM_USE_RAJA

  RAJA::RangeSegment i_range(0,Ni);
  RAJA::RangeSegment j_range(0,Nj);
  using exec_pol = typename policy_traits< ExecPolicy >::raja_2d_exec;

  RAJA::kernel< exec_pol >( RAJA::make_tuple(i_range,j_range),
                            AXOM_LAMBDA(IndexType i, IndexType j) {
    const IndexType cellIdx = i + j*jp;
    kernel( cellIdx, i, j );
  } );

#else

  constexpr bool is_serial = std::is_same< ExecPolicy, policy::serial >::value;
  AXOM_STATIC_ASSERT( is_serial );

  for( IndexType j=0; j < Nj; ++j )
  {
    const IndexType j_offset = j * jp;
    for ( IndexType i=0; i < Ni; ++i )
    {
      const IndexType cellIdx = i + j_offset;
      kernel( cellIdx, i, j );
    } // END for all i
  } // END for all j

#endif

}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_cells( xargs::ijk, const mint::Mesh* m,
                           KernelType&& kernel )
{
  // run-time checks
  SLIC_ASSERT( m != nullptr );
  SLIC_ERROR_IF( !m->isStructured(),
                 "xargs::ijk is only valid on structured meshes!" );
  SLIC_ERROR_IF( m->getDimension() != 3,
                 "xargs::ijk is only valid for 3-D structured meshes!" );

  const mint::StructuredMesh* sm =
       static_cast< const mint::StructuredMesh* >( m );

  const IndexType Ni = sm->getCellResolution( I_DIRECTION );
  const IndexType Nj = sm->getCellResolution( J_DIRECTION );
  const IndexType Nk = sm->getCellResolution( K_DIRECTION );

  const IndexType jp = sm->cellJp();
  const IndexType kp = sm->cellKp();

#ifdef AXOM_USE_RAJA

  RAJA::RangeSegment i_range(0,Ni);
  RAJA::RangeSegment j_range(0,Nj);
  RAJA::RangeSegment k_range(0,Nk);

  using exec_pol = typename policy_traits< ExecPolicy >::raja_3d_exec;

  RAJA::kernel< exec_pol >( RAJA::make_tuple( i_range, j_range, k_range ),
                          AXOM_LAMBDA(IndexType i, IndexType j, IndexType k) {

    const IndexType cellIdx = i + j*jp + k*kp;
    kernel( cellIdx, i, j, k );
  } );

#else

  constexpr bool is_serial = std::is_same< ExecPolicy, policy::serial >::value;
  AXOM_STATIC_ASSERT( is_serial );

  for ( IndexType k=0; k < Nk; ++k )
  {
    const IndexType k_offset = k * kp;
    for ( IndexType j=0; j < Nj; ++j )
    {
      const IndexType j_offset = j * jp;
      for ( IndexType i=0; i < Ni; ++i )
      {
        const IndexType cellIdx = i + j_offset + k_offset;
        kernel( cellIdx, i, j, k );
      } // END for all i
    } // END for all j
  } // END for all k

#endif

}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_cellnodes_structured( const mint::Mesh* m,
                                          KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );

  const mint::StructuredMesh* sm =
      static_cast< const mint::StructuredMesh* >( m );

  const IndexType dimension = sm->getDimension();
  const IndexType nodeJp    = sm->nodeJp();
  const IndexType nodeKp    = sm->nodeKp();
  const IndexType* offsets  = sm->getCellNodeOffsetsArray();
  const IndexType numCells  = sm->getNumberOfCells();

#ifdef AXOM_USE_RAJA

  switch( dimension )
  {
  case 1:

    using exec_pol = typename policy_traits< ExecPolicy >::raja_exec_policy;
    RAJA::forall< exec_pol >( RAJA::RangeSegment(0,numCells),
        AXOM_LAMBDA(IndexType cellIdx) {

      IndexType cell_connectivity[ 2 ] = { cellIdx,  cellIdx+1 };
      kernel( cellIdx, cell_connectivity, 2 );
    } );

    break;
  case 2:

    for_all_cells< ExecPolicy, xargs::ij >(
        m, AXOM_LAMBDA(IndexType cellIdx, IndexType i, IndexType j ) {

      const IndexType n0 = i + j * nodeJp;
      IndexType cell_connectivity[ 4 ];

      cell_connectivity[ 0 ] = n0;
      cell_connectivity[ 1 ] = n0 + offsets[ 1 ];
      cell_connectivity[ 2 ] = n0 + offsets[ 2 ];
      cell_connectivity[ 3 ] = n0 + offsets[ 3 ];

      kernel( cellIdx, cell_connectivity, 4 );

    } );

    break;
  default:
    SLIC_ASSERT( dimension==3 );

    for_all_cells< ExecPolicy, xargs::ijk >(
     m, AXOM_LAMBDA(IndexType cellIdx, IndexType i, IndexType j, IndexType k) {

      const IndexType n0 = i + j * nodeJp + k * nodeKp;
      IndexType cell_connectivity[ 8 ];

      cell_connectivity[ 0 ] = n0;
      cell_connectivity[ 1 ] = n0 + offsets[ 1 ];
      cell_connectivity[ 2 ] = n0 + offsets[ 2 ];
      cell_connectivity[ 3 ] = n0 + offsets[ 3 ];

      cell_connectivity[ 4 ] = n0 + offsets[ 4 ];
      cell_connectivity[ 5 ] = n0 + offsets[ 5 ];
      cell_connectivity[ 6 ] = n0 + offsets[ 6 ];
      cell_connectivity[ 7 ] = n0 + offsets[ 7 ];

      kernel( cellIdx, cell_connectivity, 8 );

    } );

  } // END switch

#else

  switch ( dimension )
  {
  case 1:
    {
      IndexType cell_connectivity[ 2 ];
      for ( IndexType icell=0; icell < numCells; ++icell )
      {
        cell_connectivity[ 0 ] = icell;
        cell_connectivity[ 1 ] = icell+1;
        kernel( icell, cell_connectivity, 2 );
      } // END for all cells

    } // END 3D
    break;
  case 2:
    {
      IndexType cell_connectivity[ 4 ];

      IndexType icell    = 0;
      const IndexType Ni = sm->getCellResolution( I_DIRECTION );
      const IndexType Nj = sm->getCellResolution( J_DIRECTION );
      for ( IndexType j=0; j < Nj; ++j )
      {
        const IndexType j_offset = j * nodeJp;
        for ( IndexType i=0; i < Ni; ++i )
        {
          const IndexType n0 = i + j_offset;

          cell_connectivity[ 0 ] = n0;
          cell_connectivity[ 1 ] = n0 + offsets[ 1 ];
          cell_connectivity[ 2 ] = n0 + offsets[ 2 ];
          cell_connectivity[ 3 ] = n0 + offsets[ 3 ];

          kernel( icell, cell_connectivity, 4 );

          ++icell;

        } // END for all i
      } // END for all j

    } // END 2D
    break;
  default:
    SLIC_ASSERT( dimension == 3 );
    {

      IndexType cell_connectivity[ 8 ];

      IndexType icell    = 0;
      const IndexType Ni = sm->getCellResolution( I_DIRECTION );
      const IndexType Nj = sm->getCellResolution( J_DIRECTION );
      const IndexType Nk = sm->getCellResolution( K_DIRECTION );

      for ( IndexType k=0; k < Nk; ++k )
      {
        const IndexType k_offset = k * nodeKp;
        for ( IndexType j=0; j < Nj; ++j )
        {
          const IndexType j_offset = j * nodeJp;
          for ( IndexType i=0; i < Ni; ++i )
          {
            const IndexType n0 = i + j_offset + k_offset;

            cell_connectivity[ 0 ] = n0;
            cell_connectivity[ 1 ] = n0 + offsets[ 1 ];
            cell_connectivity[ 2 ] = n0 + offsets[ 2 ];
            cell_connectivity[ 3 ] = n0 + offsets[ 3 ];

            cell_connectivity[ 4 ] = n0 + offsets[ 4 ];
            cell_connectivity[ 5 ] = n0 + offsets[ 5 ];
            cell_connectivity[ 6 ] = n0 + offsets[ 6 ];
            cell_connectivity[ 7 ] = n0 + offsets[ 7 ];

            kernel( icell, cell_connectivity, 8 );

            ++icell;

          } // END for all i
        } // END for all j
      } // END for all k

    } // END default

  } // END switch

#endif

}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_cellnodes_mixed( const mint::Mesh* m, KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );

  using UnstructuredMeshType = mint::UnstructuredMesh< mint::MIXED_SHAPE >;
  const UnstructuredMeshType* um =
     static_cast< const UnstructuredMeshType* >( m );

  const IndexType* cell_connectivity = um->getCellConnectivityArray( );
  const IndexType* cell_offsets      = um->getCellOffsetsArray( );

#ifdef AXOM_USE_RAJA

  for_all_cells< ExecPolicy >( m, AXOM_LAMBDA( IndexType cellIdx ) {
    const IndexType N = cell_offsets[ cellIdx+1 ] - cell_offsets[ cellIdx ];
    kernel( cellIdx, &cell_connectivity[ cell_offsets[ cellIdx ] ], N );
  } );

#else

  constexpr bool is_serial = std::is_same< ExecPolicy, policy::serial >::value;
  AXOM_STATIC_ASSERT( is_serial );

  const IndexType  numCells = um->getNumberOfCells();
  for ( IndexType icell=0; icell < numCells; ++icell )
  {
    const IndexType N = cell_offsets[ icell+1 ] - cell_offsets[ icell ];
    kernel( icell, &cell_connectivity[ cell_offsets[ icell ] ], N );
  } // END for all cells

#endif

}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_cellnodes_unstructured( const mint::Mesh* m,
                                            KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );

  using UnstructuredMeshType = mint::UnstructuredMesh< mint::SINGLE_SHAPE >;
  const UnstructuredMeshType* um =
      static_cast< const UnstructuredMeshType* >( m );

  const IndexType* cell_connectivity = um->getCellConnectivityArray();
  const IndexType stride             = um->getNumberOfCellNodes();

#ifdef AXOM_USE_RAJA

  for_all_cells< ExecPolicy >( m, AXOM_LAMBDA( IndexType cellIdx ) {
    kernel( cellIdx, &cell_connectivity[ cellIdx*stride ], stride );
  } );

#else

  constexpr bool is_serial = std::is_same< ExecPolicy, policy::serial >::value;
  AXOM_STATIC_ASSERT( is_serial );

  const IndexType numCells = um->getNumberOfCells();
  for ( IndexType icell=0; icell < numCells; ++icell )
  {
    kernel( icell, &cell_connectivity[ icell*stride ], stride );
  } // END for all cells

#endif
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_cells( xargs::nodeids, const mint::Mesh* m,
                           KernelType&& kernel )
{
  if ( m->isStructured() )
  {
    for_all_cellnodes_structured< ExecPolicy >(
        m, std::forward< KernelType >( kernel ) );
  }
  else if ( m->hasMixedCellTypes() )
  {
    for_all_cellnodes_mixed< ExecPolicy >(
        m, std::forward< KernelType >( kernel ) );
  }
  else
  {
    for_all_cellnodes_unstructured< ExecPolicy >(
        m, std::forward< KernelType >( kernel ) );
  }

}

} /* namespace internal */
} /* namespace mint     */
} /* namespace axom     */

#endif /* MINT_FOR_ALL_CELLS_HPP_ */
