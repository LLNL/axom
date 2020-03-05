// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_BVH_IMPL_HPP_
#define AXOM_SPIN_BVH_IMPL_HPP_

#include "axom/core/memory_management.hpp"          // for memory functions
#include "axom/core/Types.hpp"                      // for fixed bitwidth types

#include "axom/core/execution/for_all.hpp"          // for generic for_all()

// slic includes
#include "axom/slic/interface/slic.hpp"             // for SLIC macros

// linear bvh includes
#include "axom/spin/internal/linear_bvh/aabb.hpp"
#include "axom/spin/internal/linear_bvh/build_radix_tree.hpp"
#include "axom/spin/internal/linear_bvh/bvh_traverse.hpp"
#include "axom/spin/internal/linear_bvh/bvh_vtkio.hpp"
#include "axom/spin/internal/linear_bvh/BVHData.hpp"
#include "axom/spin/internal/linear_bvh/emit_bvh.hpp"
#include "axom/spin/internal/linear_bvh/QueryAccessor.hpp"
#include "axom/spin/internal/linear_bvh/TraversalPredicates.hpp"
#include "axom/spin/internal/linear_bvh/vec.hpp"

// RAJA includes
#include "RAJA/RAJA.hpp"


// C/C++ includes
#include <fstream>  // for std::ofstream
#include <sstream>  // for std::ostringstream
#include <string>   // for std::string
#include <cstring>  // for memcpy

namespace axom
{
namespace spin
{

//------------------------------------------------------------------------------
//  BVH IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//  PRIVATE HELPER METHOD IMPLEMENTATION
//------------------------------------------------------------------------------
namespace lbvh = internal::linear_bvh;

template < typename FloatType >
using vec4_t = internal::linear_bvh::Vec< FloatType, 4 >;

template < typename FloatType, int NDIMS >
using point_t = internal::linear_bvh::Vec< FloatType, NDIMS >;

/*!
 * \def BVH_PREDICATE
 *
 * \brief Macro that defines a traversal predicate functor.
 *
 * \param _predicateName the name of the predicate, e.g., `leftPredicate`
 * \param _p the primitive type
 * \param _s1 the 1st BVH segment consisting of the BVH bin information
 * \param _s2 the 2nd BVH segment consisting of the BVH bin information
 *
 * \note The BVH_PREDICATE may be instantiated within or outside a kernel.
 *
 * \note This macro is intended to be used internally by the BVH implementation.
 */
#define BVH_PREDICATE(_predicateName, _p, _s1, _s2 )    \
  auto _predicateName = [] AXOM_HOST_DEVICE( _p, _s1, _s2 )->bool

/*!
 * \def BVH_LEAF_ACTION
 *
 * \brief Macro that defines the leaf action functor for a BVH traversal.
 *
 * \param _funcName the nem of the functor, .e.g, `leafAction`
 * \param _node the BVH node ID
 * \param _leafNodes the leafNodes
 *
 * \note The BVH_LEAF_ACTION macro must be called within a kernel.
 *
 * \note This macro is intended to be used internally by the BVH implementation.
 */
#define BVH_LEAF_ACTION( _funcName, _node, _leafNodes ) \
  auto _funcName = [&]( _node, _leafNodes )->void

namespace
{

/*!
 * \brief Performs a traversal to count the candidates for each query point.
 *
 * \param [in] leftCheck functor for left bin predicate check.
 * \param [in] rightCheck functor for right bin predicate check.
 * \param [in] inner_nodes array of vec4s for the BVH inner nodes.
 * \param [in] leaf_nodes array of BVH leaf node indices
 * \param [in] N the number of user-supplied query points
 * \param [out] counts array of candidate counts for each query point.
 * \param [in] x user-supplied array of x-coordinates
 * \param [in] y user-supplied array of y-coordinates
 * \param [in] z user-supplied array of z-coordinates
 *
 * \return total_count the total count of candidates for all query points.
 */
template < int NDIMS, typename ExecSpace,
           typename LeftPredicate,
           typename RightPredicate,
           typename FloatType >
IndexType bvh_get_counts(
  LeftPredicate&& leftCheck,
  RightPredicate&& rightCheck,
  const vec4_t< FloatType >* inner_nodes,
  const int32* leaf_nodes,
  IndexType N,
  IndexType* counts,
  const FloatType* x,
  const FloatType* y,
  const FloatType* z ) noexcept
{
  // sanity checks
  SLIC_ASSERT( inner_nodes != nullptr );
  SLIC_ASSERT( leaf_nodes != nullptr );
  SLIC_ERROR_IF( counts == nullptr, "supplied null pointer for counts!" );
  SLIC_ERROR_IF( x == nullptr, "supplied null pointer for x-coordinates!" );
  SLIC_ERROR_IF( y == nullptr, "supplied null pointer for y-coordinates!" );
  SLIC_ERROR_IF( (z==nullptr && NDIMS==3),
                 "supplied null pointer for z-coordinates!" );

  // STEP 1: count number of candidates for each query point
  using reduce_pol = typename axom::execution_space< ExecSpace >::reduce_policy;
  RAJA::ReduceSum< reduce_pol, IndexType > total_count( 0 );

  using QueryAccessor = lbvh::QueryAccessor< NDIMS, FloatType >;
  for_all< ExecSpace >( N, AXOM_LAMBDA(IndexType i)
  {
    int32 count = 0;
    point_t< FloatType, NDIMS > point;
    QueryAccessor::getPoint( point, i, x, y, z );


    BVH_LEAF_ACTION( leafAction,
                      int32 AXOM_NOT_USED(current_node),
                     const int32* AXOM_NOT_USED(leaf_nodes) ) {
      count ++;
    };

    lbvh::bvh_traverse( inner_nodes,
                        leaf_nodes,
                        point,
                        leftCheck,
                        rightCheck,
                        leafAction );

    counts[ i ]  = count;
    total_count += count;

  } );

  return ( total_count.get() );
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
//  PUBLIC API IMPLEMENTATION
//------------------------------------------------------------------------------
template< int NDIMS, typename ExecSpace, typename FloatType >
BVH< NDIMS, ExecSpace, FloatType >::BVH( const FloatType* boxes,
                                         IndexType numItems ) :
  m_scaleFactor( DEFAULT_SCALE_FACTOR ),
  m_numItems( numItems ),
  m_boxes( boxes )
{}

//------------------------------------------------------------------------------
template< int NDIMS, typename ExecSpace, typename FloatType >
BVH< NDIMS, ExecSpace, FloatType >::~BVH()
{
  m_bvh.deallocate();
}

//------------------------------------------------------------------------------
template< int NDIMS, typename ExecSpace, typename FloatType >
int BVH< NDIMS, ExecSpace, FloatType >::build()
{
  // STEP 0: set the default memory allocator to use for the execution space.
  const int currentAllocatorID = axom::getDefaultAllocatorID();
  const int allocatorID = axom::execution_space< ExecSpace >::allocatorID();
  axom::setDefaultAllocator( allocatorID );

  // STEP 1: Handle case when user supplied a single bounding box
  int numBoxes        = m_numItems;
  FloatType* boxesptr = nullptr;
  if ( m_numItems == 1 )
  {
    numBoxes          = 2;
    constexpr int32 M = NDIMS * 2;      // number of entries for one box
    const int N       = numBoxes * M;   // number of entries for N boxes
    boxesptr          = axom::allocate< FloatType >( N );

    const FloatType* myboxes = m_boxes;

    // copy first box and add a fake 2nd box
    for_all< ExecSpace >( N, AXOM_LAMBDA(IndexType i)
    {
      boxesptr[ i ] = ( i < M ) ? myboxes[ i ] : 0.0;
    } );

  } // END if single item
  else
  {
    boxesptr = const_cast< FloatType* >( m_boxes );
  }

  // STEP 2: Build a RadixTree consisting of the bounding boxes, sorted
  // by their corresponding morton code.
  lbvh::RadixTree< FloatType, NDIMS > radix_tree;
  lbvh::AABB< FloatType, NDIMS > global_bounds;
  lbvh::build_radix_tree< ExecSpace >(
    boxesptr, numBoxes, global_bounds, radix_tree, m_scaleFactor );

  // STEP 3: emit the BVH data-structure from the radix tree
  m_bvh.m_bounds = global_bounds;
  m_bvh.allocate( numBoxes );

  // STEP 4: emit the BVH
  lbvh::emit_bvh< ExecSpace >( radix_tree, m_bvh );

  radix_tree.deallocate();

  // STEP 5: deallocate boxesptr if user supplied a single box
  if ( m_numItems == 1 )
  {
    SLIC_ASSERT( boxesptr != m_boxes );
    axom::deallocate( boxesptr );
  }

  // STEP 6: restore default allocator
  axom::setDefaultAllocator( currentAllocatorID );
  return BVH_BUILD_OK;
}

//------------------------------------------------------------------------------
template< int NDIMS, typename ExecSpace, typename FloatType >
void BVH< NDIMS, ExecSpace, FloatType >::getBounds( FloatType* min,
                                                    FloatType* max ) const
{
  SLIC_ASSERT( min != nullptr );
  SLIC_ASSERT( max != nullptr );
  m_bvh.m_bounds.min( min );
  m_bvh.m_bounds.max( max );
}

//------------------------------------------------------------------------------
template< int NDIMS, typename ExecSpace, typename FloatType >
void BVH< NDIMS, ExecSpace, FloatType >::find( IndexType* offsets,
                                               IndexType* counts,
                                               IndexType*& candidates,
                                               IndexType numPts,
                                               const FloatType* x,
                                               const FloatType* y,
                                               const FloatType* z ) const
{
  SLIC_ASSERT( offsets != nullptr );
  SLIC_ASSERT( counts != nullptr );
  SLIC_ASSERT( candidates == nullptr );
  SLIC_ASSERT( x != nullptr );
  SLIC_ASSERT( y != nullptr );

  // STEP 0: set the default memory allocator to use for the execution space.
  const int currentAllocatorID = axom::getDefaultAllocatorID();
  const int allocatorID = axom::execution_space< ExecSpace >::allocatorID();
  axom::setDefaultAllocator( allocatorID );

  using PointType           = point_t< FloatType, NDIMS >;
  using TraversalPredicates = lbvh::TraversalPredicates< NDIMS, FloatType >;
  using QueryAccessor       = lbvh::QueryAccessor< NDIMS, FloatType >;

  // STEP 1: count number of candidates for each query point
  const vec4_t< FloatType >* inner_nodes = m_bvh.m_inner_nodes;
  const int32* leaf_nodes  = m_bvh.m_leaf_nodes;
  SLIC_ASSERT( inner_nodes != nullptr );
  SLIC_ASSERT( leaf_nodes != nullptr );

  // STEP 2: define traversal predicates
  BVH_PREDICATE( leftPredicate,
                 const PointType &p,
                 const vec4_t< FloatType >&s1,
                 const vec4_t< FloatType >&s2 )
  {
    return TraversalPredicates::pointInLeftBin( p, s1, s2 );
  };

  BVH_PREDICATE( rightPredicate,
                 const PointType &p,
                 const vec4_t< FloatType >&s2,
                 const vec4_t< FloatType >&s3 )
  {
    return TraversalPredicates::pointInRightBin( p, s2, s3 );
  };

  // STEP 3: get counts
  int total_count = bvh_get_counts< NDIMS,ExecSpace >(
    leftPredicate, rightPredicate, inner_nodes, leaf_nodes,
    numPts, counts, x, y, z );

  using exec_policy = typename axom::execution_space< ExecSpace >::loop_policy;
  RAJA::exclusive_scan< exec_policy >(
    counts, counts+numPts, offsets, RAJA::operators::plus<IndexType>{} );

  IndexType total_candidates = static_cast< IndexType >( total_count );
  candidates = axom::allocate< IndexType >( total_candidates);

  // STEP 4: fill in candidates for each point
  for_all< ExecSpace >( numPts, AXOM_LAMBDA (IndexType i)
  {
    int32 offset = offsets[ i ];

    PointType point;
    QueryAccessor::getPoint( point, i, x, y, z );

    BVH_LEAF_ACTION( leafAction, int32 current_node, const int32* leaf_nodes ) {
      candidates[offset] = leaf_nodes[current_node];
      offset++;
    };

    lbvh::bvh_traverse( inner_nodes,
                        leaf_nodes,
                        point,
                        leftPredicate,
                        rightPredicate,
                        leafAction );

  } );

  // STEP 3: restore default allocator
  axom::setDefaultAllocator( currentAllocatorID );
}

//------------------------------------------------------------------------------
template < int NDIMS, typename ExecSpace, typename FloatType >
void BVH< NDIMS, ExecSpace, FloatType >::writeVtkFile(
  const std::string& fileName ) const
{
  std::ostringstream nodes;
  std::ostringstream cells;
  std::ostringstream levels;

  // STEP 0: Write VTK header
  std::ofstream ofs;
  ofs.open( fileName.c_str() );
  ofs << "# vtk DataFile Version 3.0\n";
  ofs << " BVHTree \n";
  ofs << "ASCII\n";
  ofs << "DATASET UNSTRUCTURED_GRID\n";

  // STEP 1: write root
  int32 numPoints = 0;
  int32 numBins   = 0;
  lbvh::write_root( m_bvh.m_bounds, numPoints, numBins,nodes,cells,levels );


  // STEP 2: traverse the BVH and dump each bin
  constexpr int32 ROOT = 0;
  lbvh::write_recursive< FloatType, NDIMS >(
    m_bvh.m_inner_nodes, ROOT, 1, numPoints, numBins, nodes, cells, levels );

  // STEP 3: write nodes
  ofs << "POINTS " << numPoints << " double\n";
  ofs << nodes.str() << std::endl;

  // STEP 4: write cells
  const int32 nnodes = (NDIMS==2) ? 4 : 8;
  ofs << "CELLS " << numBins << " " << numBins*(nnodes+1) << std::endl;
  ofs << cells.str() << std::endl;

  // STEP 5: write cell types
  ofs << "CELL_TYPES " << numBins << std::endl;
  const int32 cellType = (NDIMS==2) ? 9 : 12;
  for ( int32 i=0 ; i < numBins ; ++i )
  {
    ofs << cellType << std::endl;
  }

  // STEP 6: dump level information
  ofs << "CELL_DATA " << numBins << std::endl;
  ofs << "SCALARS level int\n";
  ofs << "LOOKUP_TABLE default\n";
  ofs << levels.str() << std::endl;
  ofs << std::endl;

  // STEP 7: close file
  ofs.close();
}

#undef BVH_PREDICATE
#undef BVH_LEAF_ACTION


} /* namespace spin */
} /* namespace axom */

#endif /* AXOM_SPIN_BVH_IMPL_HPP_ */
