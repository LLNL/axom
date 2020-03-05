// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_BVH_H_
#define AXOM_SPIN_BVH_H_

// axom core includes
#include "axom/config.hpp"       // for Axom compile-time definitions
#include "axom/core/Macros.hpp"  // for Axom macros
#include "axom/core/Types.hpp"   // for axom::IndexType

#include "axom/core/execution/execution_space.hpp"  // for execution spaces

// C/C++ includes
#include <type_traits> // for std::is_floating_point()

#if !defined(AXOM_USE_RAJA) || !defined(AXOM_USE_UMPIRE)
#error *** The spin::BVH class requires RAJA and Umpire ***
#endif

namespace axom
{
namespace spin
{

// forward declarations
namespace internal
{
namespace linear_bvh
{
template < typename FloatType, int NDIMS >
struct BVHData;
}
}

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
  void find( IndexType* offsets, IndexType* counts,
             IndexType*& candidates, IndexType numPts,
             const FloatType* x,
             const FloatType* y,
             const FloatType* z=nullptr ) const;

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

}
}

#include "axom/spin/internal/linear_bvh/BVH_impl.hpp"

#endif /* AXOM_PRIMAL_BVH_H_ */
