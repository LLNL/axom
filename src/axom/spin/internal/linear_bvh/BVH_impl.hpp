// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_BVH_IMPL_HPP_
#define AXOM_SPIN_BVH_IMPL_HPP_

#include "axom/core/Types.hpp"              // fixed bitwidth types
#include "axom/core/execution/for_all.hpp"  // for generic for_all()
#include "axom/core/memory_management.hpp"  // for memory functions
#include "axom/core/numerics/floating_point_limits.hpp"  // floating_point_limits
#include "axom/core/utilities/AnnotationMacros.hpp"      // for annotations

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include "axom/primal/operators/intersect.hpp"

// slic includes
#include "axom/slic/interface/slic.hpp"  // for SLIC macros

// linear bvh includes
#include "axom/spin/internal/linear_bvh/QueryAccessor.hpp"

// RAJA includes
#include "RAJA/RAJA.hpp"

// C/C++ includes
#include <string>  // for std::string

namespace axom
{
namespace spin
{
template <typename FloatType>
using floating_point_limits = axom::numerics::floating_point_limits<FloatType>;

//------------------------------------------------------------------------------
//  BVH IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//  PRIVATE HELPER METHOD IMPLEMENTATION
//------------------------------------------------------------------------------
namespace lbvh = internal::linear_bvh;

/*!
 * \def BVH_PREDICATE
 *
 * \brief Macro that defines a traversal predicate functor.
 *
 * \param _predicateName the name of the predicate, e.g., `leftPredicate`
 * \param _p the primitive type
 * \param _bb the AABB for the current BVH bin
 *
 * \note The BVH_PREDICATE may be instantiated within or outside a kernel.
 *
 * \note This macro is intended to be used internally by the BVH implementation.
 */
#define BVH_PREDICATE(_predicateName, _p, _bb) \
  auto _predicateName = [=] AXOM_HOST_DEVICE(_p, _bb) -> bool

//------------------------------------------------------------------------------
//  PUBLIC API IMPLEMENTATION
//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace, typename FloatType, BVHType Impl>
BVH<NDIMS, ExecSpace, FloatType, Impl>::BVH(const FloatType* boxes,
                                            IndexType numItems,
                                            int allocatorID)
  : m_AllocatorID(allocatorID)
  , m_Tolernace(floating_point_limits<FloatType>::epsilon())
  , m_scaleFactor(DEFAULT_SCALE_FACTOR)
  , m_numItems(numItems)
  , m_boxes(boxes)
{ }

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace, typename FloatType, BVHType Impl>
int BVH<NDIMS, ExecSpace, FloatType, Impl>::build()
{
  AXOM_PERF_MARK_FUNCTION("BVH::build");

  using BoxType = primal::BoundingBox<FloatType, NDIMS>;
  using PointType = primal::Point<FloatType, NDIMS>;

  // STEP 1: Handle case when user supplied a single bounding box
  int numBoxes = m_numItems;
  BoxType* boxesptr = nullptr;
  if(m_numItems == 1)
  {
    numBoxes = 2;
    boxesptr = axom::allocate<BoxType>(numBoxes, m_AllocatorID);

    const BoxType* myboxes = reinterpret_cast<const BoxType*>(m_boxes);

    // copy first box and add a fake 2nd box
    for_all<ExecSpace>(
      2,
      AXOM_LAMBDA(IndexType i) {
        if(i == 0)
        {
          boxesptr[i] = myboxes[i];
        }
        else
        {
          BoxType empty_box;
          empty_box.addPoint(PointType(0.));
          boxesptr[i] = empty_box;
        }
      });

  }  // END if single item

  const BoxType* boxes_const;
  if(boxesptr)
  {
    boxes_const = boxesptr;
  }
  else
  {
    boxes_const = reinterpret_cast<const BoxType*>(m_boxes);
  }

  m_bvh.buildImpl(boxes_const, numBoxes, m_scaleFactor, m_AllocatorID);

  // STEP 5: deallocate boxesptr if user supplied a single box
  if(m_numItems == 1)
  {
    SLIC_ASSERT(boxesptr != nullptr);
    axom::deallocate(boxesptr);
  }

  return BVH_BUILD_OK;
}

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace, typename FloatType, BVHType Impl>
void BVH<NDIMS, ExecSpace, FloatType, Impl>::getBounds(FloatType* min,
                                                       FloatType* max) const
{
  SLIC_ASSERT(min != nullptr);
  SLIC_ASSERT(max != nullptr);
  primal::BoundingBox<FloatType, NDIMS> bounds = m_bvh.getBoundsImpl();
  for(int idim = 0; idim < NDIMS; idim++)
  {
    min[idim] = bounds.getMin()[idim];
    max[idim] = bounds.getMax()[idim];
  }
}

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace, typename FloatType, BVHType Impl>
void BVH<NDIMS, ExecSpace, FloatType, Impl>::findPoints(IndexType* offsets,
                                                        IndexType* counts,
                                                        IndexType*& candidates,
                                                        IndexType numPts,
                                                        const FloatType* x,
                                                        const FloatType* y,
                                                        const FloatType* z) const
{
  AXOM_PERF_MARK_FUNCTION("BVH::findPoints");

  SLIC_ASSERT(offsets != nullptr);
  SLIC_ASSERT(counts != nullptr);
  SLIC_ASSERT(candidates == nullptr);
  SLIC_ASSERT(x != nullptr);
  SLIC_ASSERT(y != nullptr);

  using PointType = primal::Point<FloatType, NDIMS>;
  using BoundingBoxType = primal::BoundingBox<FloatType, NDIMS>;
  using QueryAccessor = lbvh::QueryAccessor<NDIMS, FloatType>;

  // STEP 1: pack query points
  PointType* packed_points = axom::allocate<PointType>(numPts, m_AllocatorID);
  for_all<ExecSpace>(
    numPts,
    AXOM_LAMBDA(IndexType i) {
      QueryAccessor::getPoint(packed_points[i], i, x, y, z);
    });

  // STEP 2: define traversal predicates
  BVH_PREDICATE(predicate, const PointType& p, const BoundingBoxType& bb)
  {
    return bb.contains(p);
  };

  m_bvh.findCandidatesImpl(predicate,
                           offsets,
                           counts,
                           candidates,
                           numPts,
                           packed_points,
                           m_AllocatorID);

  axom::deallocate(packed_points);
}

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace, typename FloatType, BVHType Impl>
void BVH<NDIMS, ExecSpace, FloatType, Impl>::findRays(IndexType* offsets,
                                                      IndexType* counts,
                                                      IndexType*& candidates,
                                                      IndexType numRays,
                                                      const FloatType* x0,
                                                      const FloatType* nx,
                                                      const FloatType* y0,
                                                      const FloatType* ny,
                                                      const FloatType* z0,
                                                      const FloatType* nz) const
{
  AXOM_PERF_MARK_FUNCTION("BVH::findRays");

  SLIC_ASSERT(offsets != nullptr);
  SLIC_ASSERT(counts != nullptr);
  SLIC_ASSERT(candidates == nullptr);
  SLIC_ASSERT(x0 != nullptr);
  SLIC_ASSERT(nx != nullptr);
  SLIC_ASSERT(y0 != nullptr);
  SLIC_ASSERT(ny != nullptr);

  const FloatType TOL = m_Tolernace;

  using RayType = primal::Ray<FloatType, NDIMS>;
  using BoundingBoxType = primal::BoundingBox<FloatType, NDIMS>;
  using QueryAccessor = lbvh::QueryAccessor<NDIMS, FloatType>;

  // STEP 1: pack query rays
  RayType* packed_rays = axom::allocate<RayType>(numRays, m_AllocatorID);
  for_all<ExecSpace>(
    numRays,
    AXOM_LAMBDA(IndexType i) {
      typename RayType::PointType origin;
      typename RayType::VectorType direction;
      QueryAccessor::getPoint(origin, i, x0, y0, z0);
      QueryAccessor::getPoint(direction, i, nx, ny, nz);

      RayType ray {origin, direction};
      packed_rays[i] = ray;
    });

  // STEP 2: define traversal predicates
  BVH_PREDICATE(predicate, const RayType& r, const BoundingBoxType& bb)
  {
    primal::Point<FloatType, NDIMS> tmp;
    return primal::detail::intersect_ray(r, bb, tmp, TOL);
  };

  m_bvh.findCandidatesImpl(predicate,
                           offsets,
                           counts,
                           candidates,
                           numRays,
                           packed_rays,
                           m_AllocatorID);

  axom::deallocate(packed_rays);
}

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace, typename FloatType, BVHType Impl>
void BVH<NDIMS, ExecSpace, FloatType, Impl>::findBoundingBoxes(
  IndexType* offsets,
  IndexType* counts,
  IndexType*& candidates,
  IndexType numBoxes,
  const FloatType* xmin,
  const FloatType* xmax,
  const FloatType* ymin,
  const FloatType* ymax,
  const FloatType* zmin,
  const FloatType* zmax) const
{
  AXOM_PERF_MARK_FUNCTION("BVH::findBoundingBoxes");

  SLIC_ASSERT(offsets != nullptr);
  SLIC_ASSERT(counts != nullptr);
  SLIC_ASSERT(candidates == nullptr);
  SLIC_ASSERT(xmin != nullptr);
  SLIC_ASSERT(xmax != nullptr);
  SLIC_ASSERT(ymin != nullptr);
  SLIC_ASSERT(ymax != nullptr);

  using BoundingBoxType = primal::BoundingBox<FloatType, NDIMS>;
  using QueryAccessor = lbvh::QueryAccessor<NDIMS, FloatType>;

  // STEP 1: pack query boxes
  BoundingBoxType* packed_boxes =
    axom::allocate<BoundingBoxType>(numBoxes, m_AllocatorID);
  for_all<ExecSpace>(
    numBoxes,
    AXOM_LAMBDA(IndexType i) {
      QueryAccessor::getBoundingBox(packed_boxes[i],
                                    i,
                                    xmin,
                                    xmax,
                                    ymin,
                                    ymax,
                                    zmin,
                                    zmax);
    });

  // STEP 2: define traversal predicates
  BVH_PREDICATE(predicate, const BoundingBoxType& bb1, const BoundingBoxType& bb2)
  {
    return bb1.intersectsWith(bb2);
  };

  m_bvh.findCandidatesImpl(predicate,
                           offsets,
                           counts,
                           candidates,
                           numBoxes,
                           packed_boxes,
                           m_AllocatorID);

  axom::deallocate(packed_boxes);
}

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace, typename FloatType, BVHType Impl>
void BVH<NDIMS, ExecSpace, FloatType, Impl>::writeVtkFile(
  const std::string& fileName) const
{
  m_bvh.writeVtkFileImpl(fileName);
}

#undef BVH_PREDICATE

} /* namespace spin */
} /* namespace axom */

#endif /* AXOM_SPIN_BVH_IMPL_HPP_ */
