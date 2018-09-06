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

#ifndef MINT_FOR_ALL_NODES_HPP_
#define MINT_FOR_ALL_NODES_HPP_

// mint includes
#include "axom/mint/execution/xargs.hpp"        // for xargs

#include "axom/mint/config.hpp"                 // for compile-time definitions
#include "axom/mint/mesh/Mesh.hpp"              // for mint::Mesh
#include "axom/mint/mesh/RectilinearMesh.hpp"   // for mint::RectilinearMesh
#include "axom/mint/mesh/StructuredMesh.hpp"    // for mint::StructuredMesh
#include "axom/mint/mesh/UniformMesh.hpp"       // for mint::UniformMesh

#include "axom/mint/execution/policy.hpp"       // for mint execution policies

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
inline void for_all_nodes( xargs::index,
                           const mint::Mesh* m,
                           KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  const IndexType numNodes = m->getNumberOfNodes();

#ifdef AXOM_USE_RAJA

  using exec_pol = typename policy_traits< ExecPolicy >::raja_exec_policy;
  RAJA::forall< exec_pol >( RAJA::RangeSegment(0,numNodes), kernel );

#else

  constexpr bool is_serial = std::is_same< ExecPolicy, policy::serial >::value;
  AXOM_STATIC_ASSERT( is_serial );

  for ( int nodeIdx=0; nodeIdx < numNodes; ++nodeIdx )
  {
    kernel( nodeIdx );
  }

#endif

}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_nodes( xargs::ij,
                           const mint::Mesh* m,
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

  const IndexType jp = sm->nodeJp();
  const IndexType Ni = sm->getNodeResolution( I_DIRECTION );
  const IndexType Nj = sm->getNodeResolution( J_DIRECTION );

#ifdef AXOM_USE_RAJA

  RAJA::RangeSegment i_range(0,Ni);
  RAJA::RangeSegment j_range(0,Nj);
  using exec_pol = typename policy_traits< ExecPolicy >::raja_2d_exec;

  RAJA::kernel< exec_pol >( RAJA::make_tuple(i_range, j_range),
                            AXOM_LAMBDA( IndexType i, IndexType j ) {
    const IndexType nodeIdx = i + j * jp;
    kernel( nodeIdx, i, j );
  } );

#else

  constexpr bool is_serial = std::is_same< ExecPolicy, policy::serial >::value;
  AXOM_STATIC_ASSERT( is_serial );

  for( IndexType j=0; j < Nj; ++j )
  {
    const IndexType j_offset = j * jp;
    for ( IndexType i=0; i < Ni; ++i )
    {
      const IndexType nodeIdx = i + j_offset;
      kernel( nodeIdx, i, j );
    } // END for all i
  } // END for all j
#endif

}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_nodes( xargs::ijk,
                           const mint::Mesh* m,
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

  const IndexType jp = sm->nodeJp();
  const IndexType kp = sm->nodeKp();
  const IndexType Ni = sm->getNodeResolution( I_DIRECTION );
  const IndexType Nj = sm->getNodeResolution( J_DIRECTION );
  const IndexType Nk = sm->getNodeResolution( K_DIRECTION );

#ifdef AXOM_USE_RAJA

  RAJA::RangeSegment i_range(0,Ni);
  RAJA::RangeSegment j_range(0,Nj);
  RAJA::RangeSegment k_range(0,Nk);
  using exec_pol = typename policy_traits< ExecPolicy >::raja_3d_exec;

  RAJA::kernel< exec_pol >( RAJA::make_tuple(i_range, j_range, k_range),
                        AXOM_LAMBDA( IndexType i, IndexType j, IndexType k ) {
    const IndexType nodeIdx = i + j*jp + k*kp;
    kernel( nodeIdx, i, j, k );
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
        const IndexType nodeIdx = i + j_offset + k_offset;
        kernel( nodeIdx, i, j, k );
      } // END for all i
    } // END for all j
  } // END for all k

#endif

}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_nodes_x( const mint::Mesh* m, KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );

  const double* x = m->getCoordinateArray( X_COORDINATE );
  SLIC_ASSERT( x != nullptr );

  const IndexType numNodes = m->getNumberOfNodes();

#ifdef AXOM_USE_RAJA

  using exec_pol = typename policy_traits< ExecPolicy >::raja_exec_policy;
  RAJA::forall< exec_pol >( RAJA::RangeSegment(0,numNodes),
                            AXOM_LAMBDA(IndexType nodeIdx) {
    kernel( nodeIdx, x[ nodeIdx ] );
  } );

#else

  constexpr bool is_serial = std::is_same< ExecPolicy, policy::serial >::value;
  AXOM_STATIC_ASSERT( is_serial );

  for ( int inode=0; inode < numNodes; ++inode )
  {
    kernel( inode, x[ inode ] );
  } // END for all nodes

#endif

}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_nodes_x_uniform( const mint::Mesh* m, KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ASSERT( m->getMeshType()==STRUCTURED_UNIFORM_MESH );

  const mint::UniformMesh* um = static_cast< const mint::UniformMesh* >( m );
  const double* x0            = um->getOrigin( );
  const double* h             = um->getSpacing( );
  const IndexType numNodes    = um->getNumberOfNodes();

#ifdef AXOM_USE_RAJA

  using exec_pol = typename policy_traits< ExecPolicy >::raja_exec_policy;
  RAJA::forall< exec_pol >( RAJA::RangeSegment(0,numNodes),
                            AXOM_LAMBDA(IndexType nodeIdx) {
    const double x = x0[ X_COORDINATE ] + nodeIdx*h[ I_DIRECTION ];
    kernel( nodeIdx, x );
  } );

#else

  constexpr bool is_serial = std::is_same< ExecPolicy, policy::serial >::value;
  AXOM_STATIC_ASSERT( is_serial );

  for ( int inode=0; inode < numNodes; ++inode )
  {
    const double x = x0[ X_COORDINATE ] + inode*h[ I_DIRECTION ];
    kernel( inode, x );
  } // END for all nodes

#endif

}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_nodes( xargs::x,
                           const mint::Mesh* m,
                           KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ERROR_IF( m->getDimension() != 1,
                 "xargs::xy is only valid for 1-D meshes" );

  const int mesh_type = m->getMeshType();
  switch( mesh_type )
  {
  case STRUCTURED_UNIFORM_MESH:
    for_all_nodes_x_uniform< ExecPolicy >(
        m, std::forward< KernelType >( kernel ) );
    break;
  default:
    for_all_nodes_x< ExecPolicy >( m, std::forward< KernelType >( kernel ) );
  } // END switch

}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_nodes_xy( const mint::Mesh* m, KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );

  const double* x = m->getCoordinateArray( X_COORDINATE );
  const double* y = m->getCoordinateArray( Y_COORDINATE );
  SLIC_ASSERT( x != nullptr );
  SLIC_ASSERT( y != nullptr );

  const IndexType numNodes = m->getNumberOfNodes();

#ifdef AXOM_USE_RAJA

  using exec_pol = typename policy_traits< ExecPolicy >::raja_exec_policy;
  RAJA::forall< exec_pol >( RAJA::RangeSegment(0,numNodes),
                            AXOM_LAMBDA(IndexType nodeIdx) {
      kernel( nodeIdx, x[ nodeIdx ], y[ nodeIdx ] );
  } );

#else

  constexpr bool is_serial = std::is_same< ExecPolicy, policy::serial >::value;
  AXOM_STATIC_ASSERT( is_serial );

  for ( int nodeIdx=0; nodeIdx < numNodes; ++nodeIdx )
  {
    kernel( nodeIdx, x[ nodeIdx ], y[ nodeIdx ] );
  }

#endif

}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_nodes_xy_uniform( const mint::Mesh* m, KernelType&& kernel)
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ASSERT( m->getMeshType()==STRUCTURED_UNIFORM_MESH );

  const mint::UniformMesh* um = static_cast< const mint::UniformMesh* >( m );
  const double* x0 = um->getOrigin( );
  const double* h  = um->getSpacing( );

#ifdef AXOM_USE_RAJA

  for_all_nodes< ExecPolicy, xargs::ij >(
      m, AXOM_LAMBDA(IndexType nodeIdx, IndexType i, IndexType j ) {
      const double x = x0[ X_COORDINATE ] + i*h[ I_DIRECTION ];
      const double y = x0[ Y_COORDINATE ] + j*h[ J_DIRECTION ];
      kernel( nodeIdx, x, y );
  } );

#else

  constexpr bool is_serial = std::is_same< ExecPolicy, policy::serial >::value;
  AXOM_STATIC_ASSERT( is_serial );

  const IndexType jp = um->jp();
  const IndexType Ni = um->getNumberOfNodesAlongDim( I_DIRECTION );
  const IndexType Nj = um->getNumberOfNodesAlongDim( J_DIRECTION );

  for( IndexType j=0; j < Nj; ++j )
  {
    const IndexType j_offset = j * jp;                // logical offset
    const double j_h         = j * h[ J_DIRECTION ];  // physical offset

    for ( IndexType i=0; i < Ni; ++i )
    {
      const IndexType nodeIdx = i + j_offset;
      const double x = x0[ X_COORDINATE ] + i * h[ I_DIRECTION ];
      const double y = x0[ Y_COORDINATE ] + j_h;

      kernel( nodeIdx, x, y );
    } // END for all i
  } // END for all j
#endif

}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_nodes_xy_rectilinear( const mint::Mesh* m,
                                          KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ASSERT( m->getMeshType()==STRUCTURED_RECTILINEAR_MESH );

  const mint::RectilinearMesh* rm =
      static_cast< const mint::RectilinearMesh* >( m );
  const double* x = rm->getCoordinateArray( X_COORDINATE );
  const double* y = rm->getCoordinateArray( Y_COORDINATE );

#ifdef AXOM_USE_RAJA

  for_all_nodes< ExecPolicy, xargs::ij >(
      m, AXOM_LAMBDA( IndexType nodeIdx, IndexType i, IndexType j ) {
      kernel( nodeIdx, x[ i ], y[ j ] );
  } );

#else

  constexpr bool is_serial = std::is_same< ExecPolicy, policy::serial >::value;
  AXOM_STATIC_ASSERT( is_serial );

  const IndexType jp = rm->jp();
  const IndexType Ni = rm->getNumberOfNodesAlongDim( I_DIRECTION );
  const IndexType Nj = rm->getNumberOfNodesAlongDim( J_DIRECTION );

  for ( IndexType j=0; j < Nj; ++j )
  {
    const IndexType j_offset = j * jp;
    for ( IndexType i=0; i < Ni; ++i )
    {
      const IndexType nodeIdx = i + j_offset;
      kernel( nodeIdx, x[ i ], y[ j ] );
    } // END for all i
  } // END for all j

#endif

}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_nodes( xargs::xy, const mint::Mesh* m,
                           KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ERROR_IF( m->getDimension() != 2,
                 "xargs::xy is only valid for 2-D meshes" );

  const int mesh_type = m->getMeshType();
  switch( mesh_type )
  {
  case STRUCTURED_RECTILINEAR_MESH:
    for_all_nodes_xy_rectilinear< ExecPolicy >(
        m, std::forward< KernelType >( kernel ) );
    break;
  case STRUCTURED_UNIFORM_MESH:
    for_all_nodes_xy_uniform< ExecPolicy >(
        m, std::forward< KernelType >( kernel ) );
    break;
  default:
    for_all_nodes_xy< ExecPolicy >( m, std::forward< KernelType >( kernel ) );
  } // END switch

}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_nodes_xyz( const mint::Mesh* m, KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );

  const double* x = m->getCoordinateArray( X_COORDINATE );
  const double* y = m->getCoordinateArray( Y_COORDINATE );
  const double* z = m->getCoordinateArray( Z_COORDINATE );
  SLIC_ASSERT( x != nullptr );
  SLIC_ASSERT( y != nullptr );
  SLIC_ASSERT( z != nullptr );

  const IndexType numNodes = m->getNumberOfNodes();

#ifdef AXOM_USE_RAJA

  using exec_pol = typename policy_traits< ExecPolicy >::raja_exec_policy;
  RAJA::forall< exec_pol >( RAJA::RangeSegment(0,numNodes),
                            AXOM_LAMBDA(IndexType nodeIdx) {
      kernel( nodeIdx, x[ nodeIdx ], y[ nodeIdx ], z[ nodeIdx ] );
  } );

#else

  constexpr bool is_serial = std::is_same< ExecPolicy, policy::serial >::value;
  AXOM_STATIC_ASSERT( is_serial );

  for ( int nodeIdx=0; nodeIdx < numNodes; ++nodeIdx )
  {
    kernel( nodeIdx, x[ nodeIdx ], y[ nodeIdx ], z[ nodeIdx ] );
  }

#endif

}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_nodes_xyz_uniform( const mint::Mesh* m,
                                       KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ASSERT( m->getMeshType()==STRUCTURED_UNIFORM_MESH );

  const mint::UniformMesh* um = static_cast< const mint::UniformMesh* >( m );
  const double* x0 = um->getOrigin( );
  const double* h  = um->getSpacing( );

#ifdef AXOM_USE_RAJA

  for_all_nodes< ExecPolicy, xargs::ijk >(
     m, AXOM_LAMBDA(IndexType nodeIdx, IndexType i, IndexType j, IndexType k) {
     const double x = x0[ X_COORDINATE ] + i*h[ I_DIRECTION ];
     const double y = x0[ Y_COORDINATE ] + j*h[ J_DIRECTION ];
     const double z = x0[ Z_COORDINATE ] + k*h[ K_DIRECTION ];
     kernel( nodeIdx, x, y, z );
  } );

#else

  constexpr bool is_serial = std::is_same< ExecPolicy, policy::serial >::value;
  AXOM_STATIC_ASSERT( is_serial );

  const IndexType jp = um->jp();
  const IndexType kp = um->kp();
  const IndexType Ni = um->getNumberOfNodesAlongDim( I_DIRECTION );
  const IndexType Nj = um->getNumberOfNodesAlongDim( J_DIRECTION );
  const IndexType Nk = um->getNumberOfNodesAlongDim( K_DIRECTION );

  for ( IndexType k=0; k < Nk; ++k )
  {
    const IndexType k_offset = k * kp;               // logical offset (k)
    const double k_h         = k * h[ K_DIRECTION ]; // physical offset (k)
    for ( IndexType j=0; j < Nj; ++j )
    {
      const IndexType j_offset = j * jp;               // logical offset (j)
      const double j_h         = j * h[ J_DIRECTION ]; // physical offset (j)
      for ( IndexType i=0; i < Ni; ++i )
      {
        const IndexType nodeIdx = i + j_offset + k_offset;
        const double x = x0[ X_COORDINATE ] + i * h[ I_DIRECTION ];
        const double y = x0[ Y_COORDINATE ] + j_h;
        const double z = x0[ Z_COORDINATE ] + k_h;

        kernel( nodeIdx, x, y, z );
      } // END for all i
    } // END for all j
  } // END for all k
#endif

}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_nodes_xyz_rectilinear( const mint::Mesh* m,
                                           KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ASSERT( m->getMeshType()==STRUCTURED_RECTILINEAR_MESH );

  const mint::RectilinearMesh* rm =
      static_cast< const mint::RectilinearMesh* >( m );

  const double* x = rm->getCoordinateArray( X_COORDINATE );
  const double* y = rm->getCoordinateArray( Y_COORDINATE );
  const double* z = rm->getCoordinateArray( Z_COORDINATE );

#ifdef AXOM_USE_RAJA

  for_all_nodes< ExecPolicy, xargs::ijk >(
     m, AXOM_LAMBDA(IndexType nodeIdx, IndexType i, IndexType j, IndexType k) {
     kernel( nodeIdx, x[ i ], y[ j ], z[ k ] );
  } );

#else

  constexpr bool is_serial = std::is_same< ExecPolicy, policy::serial >::value;
  AXOM_STATIC_ASSERT( is_serial );

  const IndexType jp = rm->jp();
  const IndexType kp = rm->kp();
  const IndexType Ni = rm->getNumberOfNodesAlongDim( I_DIRECTION );
  const IndexType Nj = rm->getNumberOfNodesAlongDim( J_DIRECTION );
  const IndexType Nk = rm->getNumberOfNodesAlongDim( K_DIRECTION );

  for ( IndexType k=0; k < Nk; ++k )
  {
    const IndexType k_offset = k * kp;
    for ( IndexType j=0; j < Nj; ++j )
    {
      const IndexType j_offset = j * jp;
      for ( IndexType i=0; i < Ni; ++ i )
      {

        const IndexType nodeIdx = i + j_offset + k_offset;
        kernel( nodeIdx, x[ i ], y[ j ], z[ k ] );

      } // END for all i
    } // END for all j
  } // END for all k
#endif

}

//------------------------------------------------------------------------------
template < typename ExecPolicy, typename KernelType >
inline void for_all_nodes( xargs::xyz, const mint::Mesh* m,
                           KernelType&& kernel )
{
  SLIC_ASSERT( m != nullptr );
  SLIC_ERROR_IF( m->getDimension() != 3,
                 "xargs::xyz is only valid for 3-D meshes" );

  const int mesh_type = m->getMeshType();
  switch( mesh_type )
  {
  case STRUCTURED_RECTILINEAR_MESH:
    for_all_nodes_xyz_rectilinear< ExecPolicy >(
        m, std::forward< KernelType >( kernel ) );
    break;
  case STRUCTURED_UNIFORM_MESH:
    for_all_nodes_xyz_uniform< ExecPolicy >(
        m, std::forward< KernelType >( kernel ) );
    break;
  default:
    for_all_nodes_xyz< ExecPolicy >( m, std::forward< KernelType >( kernel ) );
  } // END switch
}

} /* namespace internal */
} /* namespace mint     */
} /* namespace axom     */

#endif /* MINTFOR_ALL_NODES_STRUCTURED_HPP_ */
