// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
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
#include <type_traits>  // for std::is_floating_point(), std::is_same()

#if !defined(AXOM_USE_RAJA) || !defined(AXOM_USE_UMPIRE)
  #error*** The spin::BVH class requires RAJA and Umpire ***
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
template <typename FloatType, int NDIMS>
struct BVHData;
}
}  // namespace internal

/*!
 * \brief Enumerates the list of return codes for various BVH operations.
 */
enum BVHReturnCodes
{
  BVH_BUILD_FAILED = -1,  //!< indicates that generation of the BVH failed
  BVH_BUILD_OK,           //!< indicates that the BVH was generated successfully
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
 *     bvh.findPoints( offsets, counts, candidates, numPoints, qx, qy, qz );
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
template <int NDIMS, typename ExecSpace = axom::SEQ_EXEC, typename FloatType = double>
class BVH
{
private:
  // compile time checks
  AXOM_STATIC_ASSERT_MSG(((NDIMS == 2) || (NDIMS == 3)),
                         "The BVH class may be used only in 2D or 3D.");
  AXOM_STATIC_ASSERT_MSG(std::is_floating_point<FloatType>::value,
                         "A valid FloatingType must be used for the BVH.");
  AXOM_STATIC_ASSERT_MSG(
    axom::execution_space<ExecSpace>::valid(),
    "A valid execution space must be supplied to the BVH.");

public:
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
   * \param [in] allocatorID Umpire allocator ID to use (optional)
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
   * \note If an allocatorID is not specified, the code will use the default
   *  allocator ID for the execution space specified via the template argument
   *  when the BVH object is instantiated.
   * 
   * \warning The supplied boxes array must point to a buffer in a memory space
   *  that is compatible with the execution space. For example, when using
   *  CUDA_EXEC, boxes must be in unified memory or GPU memory. The code
   *  currently does not check for that.
   *
   * \pre boxes != nullptr
   * \pre numItems > 0
   */
  BVH(const FloatType* boxes,
      IndexType numItems,
      int allocatorID = axom::execution_space<ExecSpace>::allocatorID());

  /*!
   * \brief Destructor.
   */
  ~BVH();

  /*!
   * \brief Get the ID of the allocator used by the BVH.
   * \return allocatorID the ID of the allocator used by the BVH.
   */
  int getAllocatorID() const { return m_AllocatorID; };

  /*!
   * \brief Sets the scale factor for scaling the supplied bounding boxes.
   * \param [in] scale_factor the scale factor
   *
   * \note The default scale factor is set to 1.001
   */
  void setScaleFactor(FloatType scale_factor) { m_scaleFactor = scale_factor; };

  /*!
   * \brief Returns the scale factor used when constructing the BVH.
   * \return scale_factor the scale factor
   */
  FloatType getScaleFacor() const { return m_scaleFactor; };

  /*!
   * \brief Sets the tolerance used for querying the BVH.
   * \param [in] TOL the tolerance to use.
   *
   * \note Default tolerance set to floating_point_limits<FloatType>::epsilon()
   */
  void setTolerance(FloatType TOL) { m_Tolernace = TOL; };

  /*!
   * \brief Returns the tolerance value used for BVH queries.
   * \return TOL the tolerance
   */
  FloatType getTolerance() const { return m_Tolernace; };

  /*!
   * \brief Generates the BVH
   * \return status set to BVH_BUILD_OK on success.
   */
  int build();

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
  void getBounds(FloatType* min, FloatType* max) const;

  /*!
   * \brief Finds the candidate bins that contain each of the query points.
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
   *  pre-allocated by the caller before calling findPoints().
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
  void findPoints(IndexType* offsets,
                  IndexType* counts,
                  IndexType*& candidates,
                  IndexType numPts,
                  const FloatType* x,
                  const FloatType* y,
                  const FloatType* z = nullptr) const;

  /*!
   * \brief Finds the candidate bins that intersect the given rays.
   *
   * \param [out] offsets offset to the candidates array for each ray
   * \param [out] counts stores the number of candidates for each ray
   * \param [out] candidates array of candidate IDs for each ray
   * \param [in] numRays the total number of rays
   * \param [in] x0 array consisting the ray source point x-coordinates.
   * \param [in] nx array consisting the ray normal x-components.
   * \param [in] y0 array consisting the ray source point y-coordinates
   * \param [in] ny array consisting the ray normal y-components
   * \param [in] z0 array consisting the ray source point z-coorindates (in 3D)
   * \param [in] nz array consisting the ray normal z-components (in 3D)
   *
   * \note offsets and counts are arrays of size numRays that are pre-allocated
   *  by the caller, prior to the calling findRays().
   *
   * \note After the call to findRays(), the ith ray has:
   *  * counts[ i ] candidates
   *  * candidates stored in [ offsets[ i ], offsets[i]+counts[i] ]
   *
   * \pre offsets    != nullptr
   * \pre counts     != nullptr
   * \pre candidates == nullptr
   * \pre x0 != nullptr
   * \pre nx != nullptr
   * \pre y0 != nullptr
   * \pre ny != nullptr
   * \pre z0 != nullptr if dimension==3
   * \pre nz != nullptr if dimension==3
   */
  void findRays(IndexType* offsets,
                IndexType* counts,
                IndexType*& candidates,
                IndexType numRays,
                const FloatType* x0,
                const FloatType* nx,
                const FloatType* y0,
                const FloatType* ny,
                const FloatType* z0 = nullptr,
                const FloatType* nz = nullptr) const;

  /*!
   * \brief Finds the candidate bins that intersect the given bounding boxes.
   *
   * \param [out] offsets offset to the candidates array for each bounding box
   * \param [out] counts stores the number of candidates for each bounding box
   * \param [out] candidates array of candidate IDs for each bounding box
   * \param [in]  numBoxes the total number of bounding boxes
   * \param [in]  xmin array of x-coordinates of lower bounding box corner
   * \param [in]  xmax array of x-coordinates of upper bounding box corner
   * \param [in]  ymin array of y-coordinates of lower bounding box corner
   * \param [in]  ymax array of y-coordinates of upper bounding box corner
   * \param [in]  zmin array of z-coordinates of lower bounding box corner,
   *              may be nullptr if 2D
   * \param [in]  zmax array of z-coordinates of upper bounding box corner,
   *              may be nullptr if 2D
   *
   * \note offsets and counts are pointers to arrays of size numBoxes that are
   *  pre-allocated by the caller before calling findBoundingBoxes().
   *
   * \note After the call to findBoundingBoxes(), the ith bounding box has:
   *  * counts[ i ] candidates
   *  * candidates stored in [ offsets[ i ], offsets[i]+counts[i] ]
   *
   * \pre offsets    != nullptr
   * \pre counts     != nullptr
   * \pre candidates == nullptr
   * \pre xmin != nullptr
   * \pre xmax != nullptr
   * \pre ymin != nullptr
   * \pre ymax != nullptr
   * \pre zmin != nullptr if dimension==3
   * \pre zmax != nullptr if dimension==3
   */
  void findBoundingBoxes(IndexType* offsets,
                         IndexType* counts,
                         IndexType*& candidates,
                         IndexType numBoxes,
                         const FloatType* xmin,
                         const FloatType* xmax,
                         const FloatType* ymin,
                         const FloatType* ymax,
                         const FloatType* zmin = nullptr,
                         const FloatType* zmax = nullptr) const;

  /*!
   * \brief Writes the BVH to the specified VTK file for visualization.
   * \param [in] fileName the name of VTK file.
   * \note Primarily used for debugging.
   */
  void writeVtkFile(const std::string& fileName) const;

private:
  /// \name Private Members
  /// @{

  int m_AllocatorID;
  FloatType m_Tolernace;
  FloatType m_scaleFactor;
  IndexType m_numItems;
  const FloatType* m_boxes;
  internal::linear_bvh::BVHData<FloatType, NDIMS> m_bvh;

  static constexpr FloatType DEFAULT_SCALE_FACTOR = 1.001;
  /// @}

  DISABLE_COPY_AND_ASSIGNMENT(BVH);
  DISABLE_MOVE_AND_ASSIGNMENT(BVH);
};

}  // namespace spin
}  // namespace axom

#include "axom/spin/internal/linear_bvh/BVH_impl.hpp"

#endif  // AXOM_SPIN_BVH_H_
