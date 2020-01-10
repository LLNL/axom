// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_BVH_H_
#define AXOM_SPIN_BVH_H_

// axom core includes
#include "axom/config.hpp"                 // for Axom compile-time definitions
#include "axom/core/Macros.hpp"            // for Axom macros
#include "axom/core/memory_management.hpp" // for memory functions
#include "axom/core/Types.hpp"             // for fixed bitwidth types

#include "axom/core/execution/execution_space.hpp" // for execution spaces
#include "axom/core/execution/for_all.hpp"         // for generic for_all()

// slic includes
#include "axom/slic/interface/slic.hpp"    // for SLIC macros

// spin includes
#include "axom/spin/internal/linear_bvh/aabb.hpp"
#include "axom/spin/internal/linear_bvh/vec.hpp"
#include "axom/spin/internal/linear_bvh/bvh_vtkio.hpp"
#include "axom/spin/internal/linear_bvh/BVHData.hpp"
#include "axom/spin/internal/linear_bvh/emit_bvh.hpp"
#include "axom/spin/internal/linear_bvh/build_radix_tree.hpp"

// C/C++ includes
#include <fstream>  // for std::ofstream
#include <sstream>  // for std::ostringstream
#include <string>   // for std::string
#include <cstring>  // for memcpy

#if !defined(AXOM_USE_RAJA) || !defined(AXOM_USE_UMPIRE)
#error *** The spin::BVH class requires RAJA and Umpire ***
#endif

// RAJA includes
#include "RAJA/RAJA.hpp"

namespace axom
{
namespace spin
{

/*!
 * \brief Enumerates the list of return codes for various BVH operations.
 */
enum BVHReturnCodes
{
  BVH_BUILD_FAILED=-1, //!< indicates that generation of the BVH failed
  BVH_BUILD_OK,        //!< indicates that the BVH was generated successfully
};


/*!
 * \class BVH
 *
 * \brief Defines a Bounding Volume Hierarchy (BVH) spatial acceleration
 *  data structure over a set of geometric entities.
 *
 * The BVH class provides functionality for generating a hierarchical spatial
 * partitioning over a set of geometric entities. Each entity in the BVH is
 * represented by a bounding volume, in this case an axis-aligned bounding box.
 * Once the BVH structure is generated, it is used to accelerate various spatial
 * queries, such as, collision detection, ray tracing, etc., by reducing the
 * search space for a given operation to an abbreviated list of candidate
 * geometric entities to check for a particular query.
 *
 * \tparam NDIMS the number of dimensions, e.g., 2 or 3.
 * \tparam ExecSpace the execution space to use, e.g. SEQ_EXEC, CUDA_EXEC, etc.
 * \tparam FloatType floating precision, e.g., `double` or `float`. Optional.
 *
 * \note The last template parameter is optional. Defaults to double precision
 *  if not specified.
 *
 * \pre The spin::BVH class requires RAJA and Umpire. For a CPU-only, sequential
 *  implementation, see the spin::BVHTree class.
 *
 * \note The Execution Space, supplied as the 2nd template argument, specifies
 *
 *  1. Where and how the BVH is generated and stored
 *  2. Where and how subsequent queries are performed
 *  3. The default memory space, bound to the corresponding execution space
 *
 * \see axom::execution_space for more details.
 *
 *  A simple example illustrating how to use the BVH class is given below:
 *  \code
 *
 *     namespace spin = axom::spin;
 *     constexpr int DIMENSION = 3;
 *
 *     // get a list of axis-aligned bounding boxes in a flat array
 *     const double* aabbs = ...
 *
 *     // create a 3D BVH instance in parallel on the CPU using OpenMP
 *     spin::BVH< DIMENSION, axom::OMP_EXEC > bvh( aabbs, numItems );
 *     bvh.build();
 *
 *     // query points supplied in arrays, qx, qy, qz,
 *     const axom::IndexType numPoints = ...
 *     const double* qx = ...
 *     const double* qy = ...
 *     const double* qz = ...
 *
 *     // output array buffers
 *     axom::IndexType* offsets    = axom::allocate< IndexType >( numPoints );
 *     axom::IndexType* counts     = axom::allocate< IndexType >( numPoints );
 *     axom::IndexType* candidates = nullptr;
 *
 *     // find candidates in parallel, allocates and populates the supplied
 *     // candidates array
 *     bvh.find( offsets, counts, candidates, numPoints, qx, qy, qz );
 *     SLIC_ASSERT( candidates != nullptr );
 *
 *     ...
 *
 *     // caller is responsible for properly de-allocating the candidates array
 *     axom::deallocate( candidates );
 *
 *  \endcode
 *
 */
template < int NDIMS, typename ExecSpace, typename FloatType = double >
class BVH
{
public:

  AXOM_STATIC_ASSERT_MSG( ( (NDIMS==2) || (NDIMS==3) ),
                          "The BVH class may be used only in 2D or 3D." );
  AXOM_STATIC_ASSERT_MSG( std::is_floating_point< FloatType >::value,
                          "A valid FloatingType must be used for the BVH." );
  AXOM_STATIC_ASSERT_MSG( axom::execution_space< ExecSpace >::valid(),
      "A valid execution space must be supplied to the BVH." );


  /*!
   * \brief Default constructor. Disabled.
   */
  BVH() = delete;

  /*!
   * \brief Creates a BVH instance, of specified dimension, over a given set
   *  of geometric entities, each represented by its corresponding axis-aligned
   *  bounding box.
   *
   * \param [in] boxes buffer consisting of bounding boxes for each entity.
   * \param [in] numItems the total number of items to store in the BVH.
   *
   * \note boxes is an array of length 2*dimension*numItems, that stores the
   *  two corners of the axis-aligned bounding box corresponding to a given
   *  geometric entity. For example, in 3D, the two corners of the ith bounding
   *  box are given by:
   *  \code
   *    const int offset = i*6;
   *
   *    double xmin = boxes[ offset   ];
   *    double ymin = boxes[ offset+1 ];
   *    double zmin = boxes[ offset+2 ];
   *
   *    double xmax = boxes[ offset+3 ];
   *    double ymax = boxes[ offset+4 ];
   *    double zmax = boxes[ offset+5 ];
   *  \endcode
   *
   * \warning The supplied boxes array must point to a buffer in a memory space
   *  that is compatible with the execution space. For example, when using
   *  CUDA_EXEC, boxes must be in unified memory or GPU memory. The code
   *  currently does not check for that.
   *
   * \pre boxes != nullptr
   * \pre numItems > 0
   */
  BVH( const FloatType* boxes, IndexType numItems );

  /*!
   * \brief Destructor.
   */
  ~BVH();

  /*!
   * \brief Sets the scale factor for scaling the supplied bounding boxes.
   * \param [in] scale_factor the scale factor
   *
   * \note The default scale factor is set to 1.001
   */
  void setScaleFactor( FloatType scale_factor )
  { m_scaleFactor = scale_factor; };

  /*!
   * \brief Generates the BVH
   * \return status set to BVH_BUILD_OK on success.
   */
  int build( );

  /*!
   * \brief Returns the bounds of the BVH, given by the the root bounding box.
   *
   * \param [out] min buffer to store the lower corner of the root bounding box.
   * \param [out] max buffer to store the upper corner of the root bounding box.
   *
   * \note min/max point to arrays that are at least NDIMS long.
   *
   * \pre min != nullptr
   * \pre max != nullptr
   */
  void getBounds( FloatType* min, FloatType* max ) const;

  /*!
   * \brief Finds the candidate geometric entities that contain each of the
   *  given query points.
   *
   * \param [out] offsets offset to the candidates array for each query point
   * \param [out] counts stores the number of candidates per query point
   * \param [out] candidates array of the candidate IDs for each query point
   * \param [in]  numPts the total number of query points supplied
   * \param [in]  x array of x-coordinates
   * \param [in]  y array of y-coordinates
   * \param [in]  z array of z-coordinates, may be nullptr if 2D
   *
   * \note offsets and counts are pointers to arrays of size numPts that are
   *  pre-allocated by the caller before calling find().
   *
   * \note The candidates array is allocated internally by the method and
   *  ownership of the memory is transferred to the caller. Consequently, the
   *  caller is responsible for properly deallocating the candidates buffer.
   *
   * \note Upon completion, the ith query point has:
   *  * counts[ i ] candidates
   *  * Stored in the candidates array in the following range:
   *    [ offsets[ i ], offsets[ i ]+counts[ i ] ]
   *
   * \pre offsets != nullptr
   * \pre counts  != nullptr
   * \pre candidates == nullptr
   * \pre x != nullptr
   * \pre y != nullptr if dimension==2 || dimension==3
   * \pre z != nullptr if dimension==3
   */
  /// @{
  void find( IndexType* offsets, IndexType* counts, IndexType*& candidates,
             IndexType numPts, const FloatType* x, const FloatType* y ) const;
  void find( IndexType* offsets, IndexType* counts,
             IndexType*& candidates, IndexType numPts,
             const FloatType* x, const FloatType* y, const FloatType* z ) const;
  /// @}

  /*!
   * \brief Writes the BVH to the specified VTK file for visualization.
   * \param [in] fileName the name of VTK file.
   * \note Primarily used for debugging.
   */
  void writeVtkFile( const std::string& fileName ) const;

private:

/// \name Private Members
/// @{

  FloatType m_scaleFactor;
  IndexType m_numItems;
  const FloatType* m_boxes;
  internal::linear_bvh::BVHData< FloatType,NDIMS > m_bvh;

  static constexpr FloatType DEFAULT_SCALE_FACTOR = 1.001;
/// @}

  DISABLE_COPY_AND_ASSIGNMENT(BVH);
  DISABLE_MOVE_AND_ASSIGNMENT(BVH);
};

//------------------------------------------------------------------------------
//  BVH IMPLEMENTATION HELPER METHODS
//------------------------------------------------------------------------------
namespace
{

//------------------------------------------------------------------------------
template< typename FloatType >
AXOM_HOST_DEVICE
bool InLeft( const internal::linear_bvh::Vec< FloatType, 3>& point,
             const internal::linear_bvh::Vec< FloatType, 4>& s1,
             const internal::linear_bvh::Vec< FloatType, 4>& s2 )
{
  bool in_left = true;
  if ( point[0]  < s1[0] ) in_left = false;
  if ( point[1]  < s1[1] ) in_left = false;
  if ( point[2]  < s1[2] ) in_left = false;

  if ( point[0]  > s1[3] ) in_left = false;
  if ( point[1]  > s2[0] ) in_left = false;
  if ( point[2]  > s2[1] ) in_left = false;

  return in_left;
}

//------------------------------------------------------------------------------
template< typename FloatType >
AXOM_HOST_DEVICE
bool InLeft( const internal::linear_bvh::Vec< FloatType, 2>& point,
             const internal::linear_bvh::Vec< FloatType, 4>& s1,
             const internal::linear_bvh::Vec< FloatType, 4>& s2 )
{
  bool in_left = true;
  if ( point[0]  < s1[0] ) in_left = false;
  if ( point[1]  < s1[1] ) in_left = false;

  if ( point[0]  > s1[3] ) in_left = false;
  if ( point[1]  > s2[0] ) in_left = false;

  return in_left;
}

//------------------------------------------------------------------------------
template< typename FloatType >
AXOM_HOST_DEVICE
bool InRight( const internal::linear_bvh::Vec< FloatType, 3>& point,
              const internal::linear_bvh::Vec< FloatType, 4>& s2,
              const internal::linear_bvh::Vec< FloatType, 4>& s3 )
{
  bool in_right = true;

  if ( point[0] < s2[2] ) in_right = false;
  if ( point[1] < s2[3] ) in_right = false;
  if ( point[2] < s3[0] ) in_right = false;

  if ( point[0] > s3[1] ) in_right = false;
  if ( point[1] > s3[2] ) in_right = false;
  if ( point[2] > s3[3] ) in_right = false;

  return in_right;
}

//------------------------------------------------------------------------------
template< typename FloatType >
AXOM_HOST_DEVICE
bool InRight( const internal::linear_bvh::Vec< FloatType, 2>& point,
              const internal::linear_bvh::Vec< FloatType, 4>& s2,
              const internal::linear_bvh::Vec< FloatType, 4>& s3 )
{
  bool in_right = true;

  if ( point[0] < s2[2] ) in_right = false;
  if ( point[1] < s2[3] ) in_right = false;

  if ( point[0] > s3[1] ) in_right = false;
  if ( point[1] > s3[2] ) in_right = false;

  return in_right;
}

}

//------------------------------------------------------------------------------
//  BVH IMPLEMENTATION
//------------------------------------------------------------------------------
template< int NDIMS, typename ExecSpace, typename FloatType >
BVH< NDIMS, ExecSpace, FloatType >::BVH( const FloatType* boxes,
                                         IndexType numItems ) :
  m_scaleFactor( DEFAULT_SCALE_FACTOR ),
  m_numItems( numItems ),
  m_boxes( boxes )
{

}

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
  umpire::Allocator current_allocator = axom::getDefaultAllocator();
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
  internal::linear_bvh::RadixTree< FloatType, NDIMS > radix_tree;
  internal::linear_bvh::AABB< FloatType, NDIMS > global_bounds;
  internal::linear_bvh::build_radix_tree< ExecSpace >(
      boxesptr, numBoxes, global_bounds, radix_tree, m_scaleFactor );

  // STEP 3: emit the BVH data-structure from the radix tree
  m_bvh.m_bounds = global_bounds;
  m_bvh.allocate( numBoxes );

  // STEP 4: emit the BVH
  internal::linear_bvh::emit_bvh< ExecSpace >( radix_tree, m_bvh );

  radix_tree.deallocate();

  // STEP 5: deallocate boxesptr if user supplied a single box
  if ( m_numItems == 1 )
  {
    SLIC_ASSERT( boxesptr != m_boxes );
    axom::deallocate( boxesptr );
  }

  // STEP 6: restore default allocator
  axom::setDefaultAllocator( current_allocator );
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
  AXOM_STATIC_ASSERT_MSG( NDIMS==3,
     "The 3D version of find() must be called on a 3D BVH" );

  SLIC_ASSERT( offsets != nullptr );
  SLIC_ASSERT( counts != nullptr );
  SLIC_ASSERT( candidates == nullptr );
  SLIC_ASSERT( x != nullptr );
  SLIC_ASSERT( y != nullptr );
  SLIC_ASSERT( z != nullptr );

  // STEP 0: set the default memory allocator to use for the execution space.
  umpire::Allocator current_allocator = axom::getDefaultAllocator();
  const int allocatorID = axom::execution_space< ExecSpace >::allocatorID();
  axom::setDefaultAllocator( allocatorID );

  // STEP 1: count number of candidates for each query point
  using vec4_t = internal::linear_bvh::Vec< FloatType, 4 >;

  const vec4_t* inner_nodes = m_bvh.m_inner_nodes;
  const int32*  leaf_nodes  = m_bvh.m_leaf_nodes;
  SLIC_ASSERT( inner_nodes != nullptr );
  SLIC_ASSERT( leaf_nodes != nullptr );

  using reduce_policy =
      typename axom::execution_space< ExecSpace >::reduce_policy;
  RAJA::ReduceSum< reduce_policy, IndexType > total_count( 0 );

  for_all< ExecSpace >( numPts, AXOM_LAMBDA(IndexType i)
  {
    int32 count = 0;
    internal::linear_bvh::Vec< FloatType, NDIMS > point;
    point[0] = x[i];
    point[1] = y[i];
    point[2] = z[i];

   int32 current_node = 0;
   int32 todo[64];
   int32 stackptr = 0;

   constexpr int32 barrier = -2000000000;
   todo[stackptr] = barrier;
   while (current_node != barrier)
   {
     if (current_node > -1)
     {
       const vec4_t first4  = inner_nodes[current_node + 0];
       const vec4_t second4 = inner_nodes[current_node + 1];
       const vec4_t third4  = inner_nodes[current_node + 2];

       const bool in_left  = InLeft( point, first4, second4 );
       const bool in_right = InRight( point, second4, third4 );

       if (!in_left && !in_right)
       {
         // pop the stack and continue
         current_node = todo[stackptr];
         stackptr--;
       }
       else
       {
         vec4_t children = inner_nodes[current_node + 3];
         int32 l_child;
         constexpr int32 isize = sizeof(int32);
         // memcpy the int bits hidden in the floats
         memcpy(&l_child, &children[0], isize);
         int32 r_child;
         memcpy(&r_child, &children[1], isize);

         current_node = (in_left) ? l_child : r_child;

          if (in_left && in_right)
          {
            stackptr++;
            todo[stackptr] = r_child;
            // TODO: if we are in both children we could
            // go down the "closer" first by perhaps the distance
            // from the point to the center of the aabb
          }

       } // END else

     } // END if
     else
     {
       // leaf node
       count++;
       current_node = todo[stackptr];
       stackptr--;
     }

   } // while

  counts[ i ]  = count;
  total_count += count;

  } );

  using exec_policy = typename axom::execution_space< ExecSpace >::loop_policy;
  RAJA::exclusive_scan< exec_policy >(
      counts, counts+numPts, offsets, RAJA::operators::plus<IndexType>{} );

  IndexType total_candidates = static_cast< IndexType >( total_count.get() );
  candidates = axom::allocate< IndexType >( total_candidates);

  // STEP 2: fill in candidates for each point
  for_all< ExecSpace >( numPts, AXOM_LAMBDA (IndexType i)
  {
    int32 offset = offsets[ i ];

    internal::linear_bvh::Vec< FloatType,NDIMS > point;
    point[0] = x[i];
    point[1] = y[i];
    point[2] = z[i];

    int32 current_node = 0;
    int32 todo[64];
    int32 stackptr = 0;

    constexpr int32 barrier = -2000000000;
    todo[stackptr] = barrier;
    while (current_node != barrier)
    {
      if (current_node > -1)
      {
        const vec4_t first4  = inner_nodes[current_node + 0];
        const vec4_t second4 = inner_nodes[current_node + 1];
        const vec4_t third4  = inner_nodes[current_node + 2];

        const bool in_left  = InLeft( point, first4, second4 );
        const bool in_right = InRight( point, second4, third4 );

        if (!in_left && !in_right)
        {
          // pop the stack and continue
          current_node = todo[stackptr];
          stackptr--;
        }
        else
        {
          vec4_t children = inner_nodes[current_node + 3];
          int32 l_child;
          constexpr int32 isize = sizeof(int32);
          // memcpy the int bits hidden in the floats
          memcpy(&l_child, &children[0], isize);
          int32 r_child;
          memcpy(&r_child, &children[1], isize);

          current_node = (in_left) ? l_child : r_child;

          if (in_left && in_right)
          {
            stackptr++;
            todo[stackptr] = r_child;
            // TODO: if we are in both children we could
            // go down the "closer" first by perhaps the distance
            // from the point to the center of the aabb
          }
        }
      }
      else
      {
        current_node = -current_node - 1; //swap the neg address
        candidates[offset] = leaf_nodes[current_node];
        offset++;
        current_node = todo[stackptr];
        stackptr--;
      }

    } // while

  } );

  // STEP 3: restore default allocator
  axom::setDefaultAllocator( current_allocator );
}

//------------------------------------------------------------------------------
template< int NDIMS, typename ExecSpace, typename FloatType >
void BVH< NDIMS, ExecSpace, FloatType >::find( IndexType* offsets,
                                               IndexType* counts,
                                               IndexType*& candidates,
                                               IndexType numPts,
                                               const FloatType* x,
                                               const FloatType* y ) const
{
  AXOM_STATIC_ASSERT_MSG( NDIMS==2,
     "The 2D version of find() must be called on a 2D BVH" );

  SLIC_ASSERT( offsets != nullptr );
  SLIC_ASSERT( counts != nullptr );
  SLIC_ASSERT( candidates == nullptr );
  SLIC_ASSERT( x != nullptr );
  SLIC_ASSERT( y != nullptr );

  // STEP 0: set the default memory allocator to use for the execution space.
  umpire::Allocator current_allocator = axom::getDefaultAllocator();
  const int allocatorID = axom::execution_space< ExecSpace >::allocatorID();
  axom::setDefaultAllocator( allocatorID );

  // STEP 1: count number of candidates for each query point
  const internal::linear_bvh::Vec< FloatType, 4>  *inner_nodes = m_bvh.m_inner_nodes;
  const int32 *leaf_nodes = m_bvh.m_leaf_nodes;
  SLIC_ASSERT( inner_nodes != nullptr );
  SLIC_ASSERT( leaf_nodes != nullptr );

  using vec4_t      = internal::linear_bvh::Vec< FloatType, 4 >;

  using reduce_policy =
        typename axom::execution_space< ExecSpace >::reduce_policy;
  RAJA::ReduceSum< reduce_policy, IndexType > total_count( 0 );

  for_all< ExecSpace >( numPts, AXOM_LAMBDA (IndexType i)
  {
    int32 count = 0;
    internal::linear_bvh::Vec< FloatType, NDIMS > point;
    point[0] = x[i];
    point[1] = y[i];

   int32 current_node = 0;
   int32 todo[64];
   int32 stackptr = 0;

   constexpr int32 barrier = -2000000000;
   todo[stackptr] = barrier;
   while (current_node != barrier)
   {
     if (current_node > -1)
     {
       const vec4_t first4  = inner_nodes[current_node + 0];
       const vec4_t second4 = inner_nodes[current_node + 1];
       const vec4_t third4  = inner_nodes[current_node + 2];

       const bool in_left  = InLeft( point, first4, second4 );
       const bool in_right = InRight( point, second4, third4 );

       if (!in_left && !in_right)
       {
         // pop the stack and continue
         current_node = todo[stackptr];
         stackptr--;
       }
       else
       {
         vec4_t children = inner_nodes[current_node + 3];
         int32 l_child;
         constexpr int32 isize = sizeof(int32);
         // memcpy the int bits hidden in the floats
         memcpy(&l_child, &children[0], isize);
         int32 r_child;
         memcpy(&r_child, &children[1], isize);

         current_node = (in_left) ? l_child : r_child;

          if (in_left && in_right)
          {
            stackptr++;
            todo[stackptr] = r_child;
            // TODO: if we are in both children we could
            // go down the "closer" first by perhaps the distance
            // from the point to the center of the aabb
          }

       } // END else

     } // END if
     else
     {
       // leaf node
       count++;
       current_node = todo[stackptr];
       stackptr--;
     }

   } // while

   counts[ i ]  = count;
   total_count += count;

  } );

  using exec_policy = typename axom::execution_space< ExecSpace >::loop_policy;
  RAJA::exclusive_scan< exec_policy >(
      counts, counts+numPts, offsets, RAJA::operators::plus<IndexType>{} );

  IndexType total_candidates = static_cast< IndexType >( total_count.get() );

  candidates = axom::allocate< IndexType >( total_candidates);

  // STEP 2: fill in candidates for each point
  for_all< ExecSpace >( numPts, AXOM_LAMBDA (IndexType i)
  {
    int32 offset = offsets[ i ];

    internal::linear_bvh::Vec< FloatType,NDIMS > point;
    point[0] = x[i];
    point[1] = y[i];

    int32 current_node = 0;
    int32 todo[64];
    int32 stackptr = 0;

    constexpr int32 barrier = -2000000000;
    todo[stackptr] = barrier;
    while (current_node != barrier)
    {
      if (current_node > -1)
      {
        const vec4_t first4  = inner_nodes[current_node + 0];
        const vec4_t second4 = inner_nodes[current_node + 1];
        const vec4_t third4  = inner_nodes[current_node + 2];

        const bool in_left  = InLeft( point, first4, second4 );
        const bool in_right = InRight( point, second4, third4 );

        if (!in_left && !in_right)
        {
          // pop the stack and continue
          current_node = todo[stackptr];
          stackptr--;
        }
        else
        {
          vec4_t children = inner_nodes[current_node + 3];
          int32 l_child;
          constexpr int32 isize = sizeof(int32);
          // memcpy the int bits hidden in the floats
          memcpy(&l_child, &children[0], isize);
          int32 r_child;
          memcpy(&r_child, &children[1], isize);

          current_node = (in_left) ? l_child : r_child;

          if (in_left && in_right)
          {
            stackptr++;
            todo[stackptr] = r_child;
            // TODO: if we are in both children we could
            // go down the "closer" first by perhaps the distance
            // from the point to the center of the aabb
          }
        }
      }
      else
      {
        current_node = -current_node - 1; //swap the neg address
        candidates[offset] = leaf_nodes[current_node];
        offset++;
        current_node = todo[stackptr];
        stackptr--;
      }

    } // while

  } );

  // STEP 3: restore default allocator
  axom::setDefaultAllocator( current_allocator );
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
  internal::linear_bvh::write_root(
      m_bvh.m_bounds, numPoints, numBins,nodes,cells,levels );


  // STEP 2: traverse the BVH and dump each bin
  constexpr int32 ROOT = 0;
  internal::linear_bvh::write_recursive< FloatType, NDIMS >(
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
  for ( int32 i=0; i < numBins; ++i )
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

} /* namespace primal */
} /* namespace axom */


#endif /* AXOM_PRIMAL_BVH_H_ */
