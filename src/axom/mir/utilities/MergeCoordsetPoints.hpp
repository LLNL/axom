// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_MERGE_COORDSET_POINTS_HPP_
#define AXOM_MIR_MERGE_COORDSET_POINTS_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/mir/utilities/utilities.hpp"
#include "axom/mir/utilities/blueprint_utilities.hpp"
#include "axom/mir/utilities/CoordsetSlicer.hpp"

#include <conduit/conduit.hpp>

// Includes for appropriate round function.
#if defined(AXOM_DEVICE_CODE)
 #if defined(AXOM_USE_CUDA)
  #include <cuda_runtime.h>
 #endif
 #if defined(AXOM_USE_HIP)
  #include <hip/hip_runtime.h>
 #endif
#else
#include <cmath>
#endif

// RAJA
#if defined(AXOM_USE_RAJA)
  #include "RAJA/RAJA.hpp"
#endif

namespace axom
{
namespace mir
{
namespace utilities
{
namespace blueprint
{
namespace detail
{

/*!
 * \brief Round a number.
 */
template <typename value_type>
struct Rounder
{
  static AXOM_HOST_DEVICE inline value_type execute(value_type value)
  {
#if defined(AXOM_DEVICE_CODE)
    return round(value);
#else
    return std::round(value);
#endif
  }
};

/*!
 * \brief Round a float number.
 */
template <>
struct Rounder<float>
{
  static AXOM_HOST_DEVICE inline float execute(float value)
  {
#if defined(AXOM_DEVICE_CODE)
    return roundf(value);
#else
    return std::round(value);
#endif
  }
};

} // end namespace detail

/**
 * \brief Merge
 *
 * \tparam ExecSpace The execution space where the algorithm runs.
 * \tparam CoordsetView The coordset view type.
 */
template <typename ExecSpace, typename CoordsetView>
class MergeCoordsetPoints
{
public:
  using value_type = typename CoordsetView::value_type;

  /*!
   * \brief Constructor
   *
   * \param coordsetView The coordset view that wraps the coordset to be modified.
   */
  MergeCoordsetPoints(const CoordsetView &coordsetView) : m_coordsetView(coordsetView)
  {
  }

  /*!
   * \brief Merge the coordset points using a tolerance and pass out an array of the
   *        points that were selected from the original coordset since they can be
   *        used to slice nodal fields.
   *
   * \param[inout] n_coordset The Conduit node that contains the coordset. The coordset data
   *                          will be replaced.
   * \param n_options A Conduit node that contains the supported options. At present, there
   *                  is just a "tolerance" option that specifies how close points need to
   *                  be in order to be matched.
   * \param[out] selectedIds An array containing the ids of the nodes that are selected from
   *                         the old nodes.
   * \param[out] old2new An array with nnodes elements that contains the new node id for each
   *                     node in the old coordset. This can be useful in rewriting connectivity
   *                     or fields.
   *
   * \note The algorithm eliminates a little precision on the input coordinates to make a
   *       hashed name for the point. Points that are close enough should have the same
   *       hashed name. The names are passed through a Unique filter to select a set of
   *       unique names. Multiple points could result in the same hash name so one of the
   *       points is selected as the new definition for the unique point (the selection
   *       arises from sorting the hash names). This means that among a bunch of points
   *       that are close enough to be grouped, one of them will be picked and it might
   *       not be the best one "esthetically". The points are output in sorted hash name
   *       order too, which is a bit random.
   *
   * \return True if point merging happened; False if no point merging was needed.
   */
  bool execute(conduit::Node &n_coordset,
               const conduit::Node &n_options,
               axom::Array<axom::IndexType> &selectedIds,
               axom::Array<axom::IndexType> &old2new) const
  {
    namespace bputils = axom::mir::utilities::blueprint;
    AXOM_ANNOTATE_SCOPE("MergeCoordsetPoints");
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    // If the coordset is not explicit then there is nothing to do.
    if(n_coordset["type"].as_string() != "explicit")
    {
      SLIC_ERROR("The input coordset was not explicit.");
      return false;
    }

    // Get any options.
    const value_type DEFAULT_TOLERANCE = 1.e-10;
    value_type tolerance = DEFAULT_TOLERANCE;
    if(n_options.has_child("tolerance"))
    {
      tolerance = static_cast<value_type>(n_options["tolerance"].to_double());
    }

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("naming");
    using KeyType = std::uint64_t;

    const auto nnodes = m_coordsetView.numberOfNodes();

    axom::Array<KeyType> coordNames(nnodes, nnodes, allocatorID);
    auto coordNamesView = coordNames.view();
    const auto deviceCoordsetView = m_coordsetView;
    axom::for_all<ExecSpace>(nnodes, AXOM_LAMBDA(axom::IndexType index)
    {
      // Get the current point.
      const auto pt = deviceCoordsetView[index];

      // Truncate the point components using the tolerance so we can eliminate
      // precision beyond what we want with the tolerance. The idea is that points
      // that are close enough will hash to the same name.
      value_type truncated[CoordsetView::dimension()];
      for(int d = 0; d < CoordsetView::dimension(); d++)
      {
        value_type value = detail::Rounder<value_type>::execute(pt[d] / tolerance) * tolerance;
        if(value > -tolerance && value < tolerance)
        {
          value = value_type{0};
        }
        truncated[d] = value;
      }

      // Make a name for this point
      const void *tptr = static_cast<const void *>(truncated);
      coordNamesView[index] = axom::mir::utilities::hash_bytes(static_cast<const std::uint8_t *>(tptr), sizeof(value_type) * CoordsetView::dimension());
    });
    AXOM_ANNOTATE_END("naming");

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("unique");

    // Make faces unique.
    axom::Array<KeyType> uniqueNames;
    axom::mir::utilities::Unique<ExecSpace, KeyType>::execute(coordNamesView, uniqueNames, selectedIds);
    const auto uniqueNamesView = uniqueNames.view();
    const auto selectedIdsView = selectedIds.view();
    AXOM_ANNOTATE_END("unique");

    //--------------------------------------------------------------------------
    bool merged = false;
    old2new = axom::Array<axom::IndexType>(nnodes, nnodes, allocatorID);
    auto old2newView = old2new.view();
    if(selectedIds.size() < nnodes)
    {
      // There are fewer nodes in the selectedIds so we are able to combine nodes.

      AXOM_ANNOTATE_BEGIN("old2new");
      // Make a map of nodes in the old coordset to nodes in the new coordset. We
      // do it by looking up the old node name in the new coordset unique names.
      axom::for_all<ExecSpace>(nnodes, AXOM_LAMBDA(axom::IndexType index)
      {
        const auto newNodeId = axom::mir::utilities::bsearch(coordNamesView[index], uniqueNamesView);
        SLIC_ASSERT(newNodeId >= 0 && newNodeId < nnodes);
        old2newView[index] = newNodeId;
      });
      AXOM_ANNOTATE_END("old2new");

      //--------------------------------------------------------------------------
      // Use the selectedIds to slice the coordset to make a new coordset that
      // replaces the old one.
      bputils::CoordsetSlicer<ExecSpace, CoordsetView> css(m_coordsetView);
      bputils::SliceData slice;
      slice.m_indicesView = selectedIdsView;
      conduit::Node n_sliced;
      css.execute(slice, n_coordset, n_sliced);
      n_coordset.move(n_sliced);

      merged = true;
    }
    else
    {
      // It does not look like there is a need for point merging.
      AXOM_ANNOTATE_BEGIN("old2new");

      auto selectedIdsView = selectedIds.view();
      axom::for_all<ExecSpace>(nnodes, AXOM_LAMBDA(axom::IndexType index)
      {
        selectedIdsView[index] = index;
        old2newView[index] = index;
      });
    }

    return merged;
  }

  CoordsetView m_coordsetView;
};

}  // end namespace blueprint
}  // end namespace utilities
}  // end namespace mir
}  // end namespace axom

#endif
