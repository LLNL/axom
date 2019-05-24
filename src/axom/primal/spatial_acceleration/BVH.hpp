// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_BVH_H_
#define AXOM_PRIMAL_BVH_H_

#include "axom/config.hpp"                 // for Axom compile-time definitions
#include "axom/core/Macros.hpp"            // for Axom macros
#include "axom/core/memory_management.hpp" // for memory functions
#include "axom/core/Types.hpp"             // for fixed bitwidth types
#include "axom/slic/interface/slic.hpp"    // for SLIC macros

#include "axom/primal/spatial_acceleration/linear_bvh/bvh_builder.hpp"
#include "axom/primal/spatial_acceleration/linear_bvh/bvh_traverse.hpp"
#include "axom/primal/spatial_acceleration/linear_bvh/bvh_vtkio.hpp"

// C/C++ includes
#include <fstream>  // for std::ofstream
#include <sstream>  // for std::ostringstream
#include <string>   // for std::string

#ifdef AXOM_USE_RAJA
#include "RAJA/RAJA.hpp"
#endif

namespace axom
{
namespace primal
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
 * geometric entities to check for a given query.
 *
 */
template < int NDIMS, typename FloatType = double >
class BVH
{
public:

  AXOM_STATIC_ASSERT_MSG( ( (NDIMS==2) || (NDIMS==3) ),
                          "The BVH class may be used only in 2D or 3D." );
  AXOM_STATIC_ASSERT_MSG( std::is_floating_point< FloatType >::value,
                          "A valid FloatingType must be used for the BVH." );

  /*!
   * \brief Default constructor. Disabled.
   */
  BVH() = delete;

  /*!
   * \brief Creates a BVH instance, of specified dimension, over a given set
   *  of geometric entities, each represented by its corresponding axis-aligned
   *  bounding box.
   *
   * \param [in] dimension the spatial dimension for the BVH
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
   * \pre dimension >= 1 && dimension <= 3
   * \pre boxes != nullptr
   * \pre numItems > 0
   */
  BVH( const FloatType* boxes, IndexType numItems );

  /*!
   * \brief Destructor.
   */
  ~BVH();

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
   * \param [out] candidates array of the candidate IDs for each query point
   * \param [in]  numPts the total number of query points supplied
   * \param [in]  x array of x-coordinates
   * \param [in]  y array of y-coordinates, may be nullptr if 1D
   * \param [in]  z array of z-coordinates, may be nullptr if 1D or 2D
   *
   * \note the offsets array must be pre-allocated by the caller and
   *
   * \note The candidates array is allocated internally by the method and
   *  ownership of the memory is transferred to the caller. Consequently, the
   *  caller is responsible for properly deallocating the candidates buffer.
   *
   * \pre offsets != nullptr
   * \pre candidates == nullptr
   * \pre x != nullptr
   * \pre y != nullptr if dimension==2 || dimension==3
   * \pre z != nullptr if dimension==3
   */
  /// @{
  void find( IndexType* offsets, IndexType*& candidates, IndexType numPts,
             const FloatType* x, const FloatType* y ) const;
  void find( IndexType* offsets, IndexType*& candidates, IndexType numPts,
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

  IndexType m_numItems;
  const FloatType* m_boxes;
  bvh::BVH< FloatType,NDIMS > m_bvh;

/// @}

  DISABLE_COPY_AND_ASSIGNMENT(BVH);
  DISABLE_MOVE_AND_ASSIGNMENT(BVH);
};

//------------------------------------------------------------------------------
//  BVH IMPLEMENTATION
//------------------------------------------------------------------------------
template< int NDIMS, typename FloatType >
BVH< NDIMS, FloatType >::BVH( const FloatType* boxes, IndexType numItems ) :
  m_numItems( numItems ),
  m_boxes( boxes )
{

}

//------------------------------------------------------------------------------
template< int NDIMS, typename FloatType >
BVH< NDIMS, FloatType >::~BVH()
{
  m_bvh.free();
}

//------------------------------------------------------------------------------
template< int NDIMS, typename FloatType >
int BVH< NDIMS, FloatType >::build()
{
  bvh::LinearBVHBuilder builder;
  m_bvh = builder.construct< FloatType, NDIMS >( m_boxes, m_numItems );
  return BVH_BUILD_OK;
}

//------------------------------------------------------------------------------
template< int NDIMS, typename FloatType >
void BVH< NDIMS, FloatType >::getBounds( FloatType* min, FloatType* max ) const
{
  SLIC_ASSERT( min != nullptr );
  SLIC_ASSERT( max != nullptr );
  m_bvh.m_bounds.min( min );
  m_bvh.m_bounds.max( max );
}

//------------------------------------------------------------------------------
template< int NDIMS, typename FloatType >
void BVH< NDIMS, FloatType >::find( IndexType* offsets,
                                    IndexType*& candidates,
                                    IndexType numPts,
                                    const FloatType* x,
                                    const FloatType* y,
                                    const FloatType* z ) const
{
  AXOM_STATIC_ASSERT_MSG( NDIMS==3,
      "The 3D version of find() must be called on a 3D BVH" );

  SLIC_ASSERT( offsets != nullptr );
  SLIC_ASSERT( candidates == nullptr );

  // create local reference to member to capture by value on the device,
  // otherwise, the capture would attempt to capture (this)
  auto const& mybvh = this->m_bvh;

  // STEP 0: count candidates
  IndexType* candidate_counts = axom::allocate< IndexType >( numPts );

  using exec_pol = bvh::raja_for_policy;
  RAJA::forall< exec_pol >(
      RAJA::RangeSegment(0,numPts), AXOM_LAMBDA(IndexType i)
  {

    int32 count = 0;

    bvh::Vec< FloatType, NDIMS > point;
    point[ 0 ] = x[ i ];
    point[ 1 ] = y[ i ];
    point[ 2 ] = z[ i ];

    auto inLeft = [&] AXOM_DEVICE ( const bvh::Vec< FloatType,4 >& s1,
                                    const bvh::Vec< FloatType,4 >& s2  )
    {
      bool in_left = true;
      if(point[0]  < s1[0]) in_left = false; // L.x min
      if(point[1]  < s1[1]) in_left = false; // L.y min
      if(point[2]  < s1[2]) in_left = false; // L.z min

      if(point[0]  > s2[3]) in_left = false; // L.x max
      if(point[1]  > s2[0]) in_left = false; // L.y max
      if(point[2]  > s2[1]) in_left = false; // L.z max
      return in_left;
    };

    auto inRight = [&] AXOM_DEVICE ( const bvh::Vec< FloatType,4 >& s2,
                                     const bvh::Vec< FloatType,4 >& s3 )
    {
      bool in_right = true;
      if(point[0]  < s2[2]) in_right = false; // R.x min
      if(point[1]  < s2[3]) in_right = false; // R.y min
      if(point[2]  < s3[0]) in_right = false; // R.z min

      if(point[0]  > s3[1]) in_right = false; // R.x max
      if(point[1]  > s3[2]) in_right = false; // R.y max
      if(point[2]  > s3[3]) in_right = false; // R.z max
      return in_right;
    };

    auto leafKernel = [&] AXOM_DEVICE (
                        int32 AXOM_NOT_USED(current_node),
                        const bvh::BVH< FloatType,NDIMS >& AXOM_NOT_USED(bvh) )
    {
      ++count;
    };

    bvh_traverse< NDIMS, FloatType >( mybvh, inLeft, inRight, leafKernel );

    candidate_counts[ i ] = count;
  } );

  // STEP 1: prefix sum of candidate counts
  RAJA::exclusive_scan< bvh::raja_for_policy >(
      candidate_counts, candidate_counts + numPts, offsets,
      RAJA::operators::plus<IndexType>{} );

  // STEP 2: populate candidates
  // TODO: this will segault with raw(unmanaged) cuda pointers
  IndexType total_candidates = offsets[numPts-1] + candidate_counts[numPts-1];

  candidates = axom::allocate< IndexType >( total_candidates );

  RAJA::forall< exec_pol >(
        RAJA::RangeSegment(0,numPts), AXOM_LAMBDA(IndexType i)
    {

      int32 offset = offsets[ i ];

      bvh::Vec< FloatType, NDIMS > point;
      point[ 0 ] = x[ i ];
      point[ 1 ] = y[ i ];
      point[ 2 ] = z[ i ];

      auto inLeft = [&] AXOM_DEVICE ( const bvh::Vec< FloatType,4 >& s1,
                                      const bvh::Vec< FloatType,4 >& s2  )
      {
        bool in_left = true;
        if(point[0]  < s1[0]) in_left = false; // L.x min
        if(point[1]  < s1[1]) in_left = false; // L.y min
        if(point[2]  < s1[2]) in_left = false; // L.z min

        if(point[0]  > s2[3]) in_left = false; // L.x max
        if(point[1]  > s2[0]) in_left = false; // L.y max
        if(point[2]  > s2[1]) in_left = false; // L.z max
        return in_left;
      };

      auto inRight = [&] AXOM_DEVICE ( const bvh::Vec< FloatType,4 >& s2,
                                       const bvh::Vec< FloatType,4 >& s3 )
      {
        bool in_right = true;
        if(point[0]  < s2[2]) in_right = false; // R.x min
        if(point[1]  < s2[3]) in_right = false; // R.y min
        if(point[2]  < s3[0]) in_right = false; // R.z min

        if(point[0]  > s3[1]) in_right = false; // R.x max
        if(point[1]  > s3[2]) in_right = false; // R.y max
        if(point[2]  > s3[3]) in_right = false; // R.z max
        return in_right;
      };

      auto leafKernel = [&] AXOM_DEVICE (
                          int32 current_node,
                          const bvh::BVH< FloatType,NDIMS >& mybvh )
      {
        candidates[ offset ] = mybvh.m_leaf_nodes[ current_node ];
        ++offset;
      };

      bvh_traverse< NDIMS, FloatType >( mybvh, inLeft, inRight, leafKernel );

    } );

}

//------------------------------------------------------------------------------
template< int NDIMS, typename FloatType >
void BVH< NDIMS, FloatType >::find( IndexType* offsets,
                                    IndexType*& candidates,
                                    IndexType numPts,
                                    const FloatType* x,
                                    const FloatType* y ) const
{
  AXOM_STATIC_ASSERT_MSG( NDIMS==2,
        "The 2D version of find() must be called on a 2D BVH" );

  SLIC_ASSERT( offsets != nullptr );
  SLIC_ASSERT( candidates == nullptr );

  // create local reference to member to capture by value on the device,
  // otherwise, the capture would attempt to capture (this)
  auto const& mybvh = this->m_bvh;

  // STEP 0: count candidates
  IndexType* candidate_counts = axom::allocate< IndexType >( numPts );

  using exec_pol = bvh::raja_for_policy;
  RAJA::forall< exec_pol >(
      RAJA::RangeSegment(0,numPts), AXOM_LAMBDA(IndexType i)
  {

    int32 count = 0;

    bvh::Vec< FloatType, NDIMS > point;
    point[ 0 ] = x[ i ];
    point[ 1 ] = y[ i ];

    auto inLeft = [&] AXOM_DEVICE ( const bvh::Vec< FloatType,4 >& s1,
                                    const bvh::Vec< FloatType,4 >& s2  )
    {
      bool in_left = true;
      if(point[0]  < s1[0]) in_left = false; // L.x min
      if(point[1]  < s1[1]) in_left = false; // L.y min

      if(point[0]  > s2[3]) in_left = false; // L.x max
      if(point[1]  > s2[0]) in_left = false; // L.y max
      return in_left;
    };

    auto inRight = [&] AXOM_DEVICE ( const bvh::Vec< FloatType,4 >& s2,
                                     const bvh::Vec< FloatType,4 >& s3 )
    {
      bool in_right = true;
      if(point[0]  < s2[2]) in_right = false; // R.x min
      if(point[1]  < s2[3]) in_right = false; // R.y min

      if(point[0]  > s3[1]) in_right = false; // R.x max
      if(point[1]  > s3[2]) in_right = false; // R.y max
      return in_right;
    };

    auto leafKernel = [&] AXOM_DEVICE (
                       int32 AXOM_NOT_USED(current_node),
                       const bvh::BVH< FloatType,NDIMS >& AXOM_NOT_USED(mybvh))
    {
      ++count;
    };

    bvh_traverse< NDIMS, FloatType >( mybvh, inLeft, inRight, leafKernel );

    candidate_counts[ i ] = count;
  } );

  // STEP 1: prefix sum of candidate counts
  RAJA::exclusive_scan< bvh::raja_for_policy >(
      candidate_counts, candidate_counts + numPts, offsets,
      RAJA::operators::plus<IndexType>{} );

  // STEP 2: populate candidates
  // TODO: this will segault with raw(unmanaged) cuda pointers
  IndexType total_candidates = offsets[numPts-1] + candidate_counts[numPts-1];

  candidates = axom::allocate< IndexType >( total_candidates );

  RAJA::forall< exec_pol >(
        RAJA::RangeSegment(0,numPts), AXOM_LAMBDA(IndexType i)
  {

    int32 offset = offsets[ i ];

    bvh::Vec< FloatType, NDIMS > point;
    point[ 0 ] = x[ i ];
    point[ 1 ] = y[ i ];

    auto inLeft = [&] AXOM_DEVICE ( const bvh::Vec< FloatType,4 >& s1,
                                    const bvh::Vec< FloatType,4 >& s2  )
    {
      bool in_left = true;
      if(point[0]  < s1[0]) in_left = false; // L.x min
      if(point[1]  < s1[1]) in_left = false; // L.y min

      if(point[0]  > s2[3]) in_left = false; // L.x max
      if(point[1]  > s2[0]) in_left = false; // L.y max
      return in_left;
    };

    auto inRight = [&] AXOM_DEVICE ( const bvh::Vec< FloatType,4 >& s2,
                                     const bvh::Vec< FloatType,4 >& s3 )
    {
      bool in_right = true;
      if(point[0]  < s2[2]) in_right = false; // R.x min
      if(point[1]  < s2[3]) in_right = false; // R.y min

      if(point[0]  > s3[1]) in_right = false; // R.x max
      if(point[1]  > s3[2]) in_right = false; // R.y max
      return in_right;
    };

    auto leafKernel = [&] AXOM_DEVICE (
                        int32 current_node,
                        const bvh::BVH< FloatType,NDIMS >& mybvh )
    {
      candidates[ offset ] = mybvh.m_leaf_nodes[ current_node ];
      ++offset;
    };

    bvh_traverse< NDIMS, FloatType >( mybvh, inLeft, inRight, leafKernel );

  } );

}

//------------------------------------------------------------------------------
template < int NDIMS, typename FloatType >
void BVH< NDIMS, FloatType >::writeVtkFile( const std::string& fileName ) const
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
  bvh::write_root( m_bvh.m_bounds, numPoints, numBins,nodes,cells,levels );


  // STEP 2: traverse the BVH and dump each bin
  constexpr int32 ROOT = 0;
  bvh::write_recursive< FloatType, NDIMS >(
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
