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

#if defined(AXOM_USE_RAJA)
  // RAJA includes
  #include "RAJA/RAJA.hpp"
#endif

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
  , m_boxes(reinterpret_cast<const BoxType*>(boxes))
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

    const BoxType* myboxes = m_boxes;

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
    boxes_const = m_boxes;
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
                                                        const PointType* pts) const
{
  AXOM_PERF_MARK_FUNCTION("BVH::findPoints");

  SLIC_ASSERT(offsets != nullptr);
  SLIC_ASSERT(counts != nullptr);
  SLIC_ASSERT(candidates == nullptr);
  SLIC_ASSERT(pts != nullptr);

  // Define traversal predicates
  BVH_PREDICATE(predicate, const PointType& p, const BoxType& bb)
  {
    return bb.contains(p);
  };

  m_bvh.findCandidatesImpl(predicate,
                           offsets,
                           counts,
                           candidates,
                           numPts,
                           pts,
                           m_AllocatorID);
}

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace, typename FloatType, BVHType Impl>
void BVH<NDIMS, ExecSpace, FloatType, Impl>::findRays(IndexType* offsets,
                                                      IndexType* counts,
                                                      IndexType*& candidates,
                                                      IndexType numRays,
                                                      const RayType* rays) const
{
  AXOM_PERF_MARK_FUNCTION("BVH::findRays");

  SLIC_ASSERT(offsets != nullptr);
  SLIC_ASSERT(counts != nullptr);
  SLIC_ASSERT(candidates == nullptr);
  SLIC_ASSERT(rays != nullptr);

  const FloatType TOL = m_tolerance;

  // Define traversal predicates
  BVH_PREDICATE(predicate, const RayType& r, const BoxType& bb)
  {
    primal::Point<FloatType, NDIMS> tmp;
    return primal::detail::intersect_ray(r, bb, tmp, TOL);
  };

  m_bvh.findCandidatesImpl(predicate,
                           offsets,
                           counts,
                           candidates,
                           numRays,
                           rays,
                           m_AllocatorID);
}

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace, typename FloatType, BVHType Impl>
void BVH<NDIMS, ExecSpace, FloatType, Impl>::findBoundingBoxes(
  IndexType* offsets,
  IndexType* counts,
  IndexType*& candidates,
  IndexType numBoxes,
  const BoxType* boxes) const
{
  AXOM_PERF_MARK_FUNCTION("BVH::findBoundingBoxes");

  SLIC_ASSERT(offsets != nullptr);
  SLIC_ASSERT(counts != nullptr);
  SLIC_ASSERT(candidates == nullptr);
  SLIC_ASSERT(boxes != nullptr);

  // STEP 2: define traversal predicates
  BVH_PREDICATE(predicate, const BoxType& bb1, const BoxType& bb2)
  {
    return bb1.intersectsWith(bb2);
  };

  m_bvh.findCandidatesImpl(predicate,
                           offsets,
                           counts,
                           candidates,
                           numBoxes,
                           boxes,
                           m_AllocatorID);
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
