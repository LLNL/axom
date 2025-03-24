// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_MERGE_COORDSET_POINTS_HPP_
#define AXOM_MIR_MERGE_COORDSET_POINTS_HPP_

#include "axom/core.hpp"
#include "axom/mir/utilities/utilities.hpp"
#include "axom/mir/utilities/blueprint_utilities.hpp"
#include "axom/mir/utilities/CoordsetSlicer.hpp"

#include <conduit/conduit.hpp>

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
/**
 * \brief Merge
 *
 * \tparam ExecSpace The execution space where the algorithm runs.
 * \tparam TopologyView The topology view type.
 */
template <typename ExecSpace, typename CoordsetView>
class MergeCoordsetPoints
{
public:

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
   *                     node in the old coordset. This can be useful in rewriting connectivity.
   */
  void execute(conduit::Node &n_coordset,
               const conduit::Node &n_options,
               axom::Array<axom::IndexType> &selectedIds,
               axom::Array<axom::IndexType> &old2new) const
  {
    namespace bputils = axom::mir::utilities::blueprint;
    AXOM_ANNOTATE_BEGIN("MergeCoordsetPoints");
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    // If the coordset is not explicit then there is nothing to do.
    if(n_coordset["type"].as_string() != "explicit")
    {
      return;
    }

    // Get any options.
    const double DEFAULT_TOLERANCE = 1.e-10;
    double tolerance = DEFAULT_TOLERANCE;
    if(n_options.has_child("tolerance"))
    {
      tolerance = n_options["tolerance"].to_double();
    }

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("naming");
    using value_type = typename CoordsetView::value_type;
    using KeyType = std::uint64_t;

    const auto nnodes = m_coordsetView.numberOfNodes();
std::cout << "Creating names for " << nnodes << " nodes." << std::endl;
    axom::Array<KeyType> coordNames(nnodes, nnodes, allocatorID);
    auto coordNamesView = coordNames.view();
    const auto deviceCoordsetView = m_coordsetView;
    axom::for_all<ExecSpace>(nnodes, AXOM_LAMBDA(axom::IndexType index)
    {
      // Get the current point.
      const auto pt = deviceCoordsetView[index];

      // Truncate the point components using the tolerance
      value_type truncated[CoordsetView::dimension()];
      for(int d = 0; d < CoordsetView::dimension(); d++)
      {
        value_type value = static_cast<value_type>(std::round(pt[d] / tolerance) * tolerance);
        if(value > -tolerance && value < tolerance)
        {
          value = value_type{0};
        }
        truncated[d] = value;
      }

      // Make a name for this point
      const void *tptr = static_cast<const void *>(truncated);
      coordNamesView[index] = axom::mir::utilities::hash_bytes(static_cast<const std::uint8_t *>(tptr), sizeof(value_type) * CoordsetView::dimension());
std::cout << "coordNames[" << index << "]=" << coordNamesView[index] << std::endl;
    });
    AXOM_ANNOTATE_END("naming");

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("unique");
    // The unique keys.
    axom::Array<KeyType> uniqueNames;

    // Make faces unique.
    axom::mir::utilities::Unique<ExecSpace, KeyType>::execute(coordNamesView, uniqueNames, selectedIds);
    const auto uniqueNamesView = uniqueNames.view();
    const auto selectedIdsView = selectedIds.view();
    AXOM_ANNOTATE_END("unique");
std::cout << "unique: selectedIds.size=" << selectedIds.size() << ", uniqueNames.size=" << uniqueNames.size() << std::endl;

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("old2new");
    old2new = axom::Array<axom::IndexType>(nnodes, nnodes, allocatorID);
std::cout << "creating old2new sized " << nnodes << " nodes." << std::endl;
    auto old2newView = old2new.view();
    axom::for_all<ExecSpace>(nnodes, AXOM_LAMBDA(axom::IndexType index)
    {
      const auto newNodeId = axom::mir::utilities::bsearch(coordNamesView[index], uniqueNamesView);
      SLIC_ASSERT(newNodeId >= 0 && newNodeId < nnodes);
      old2newView[index] = newNodeId;
std::cout << "old2new[" << index << "]=" << newNodeId << std::endl;
    });
    AXOM_ANNOTATE_END("old2new");
    
    //--------------------------------------------------------------------------
    // Use the selectedIds to slice the coordset to make a new coordset that
    // replaces the old one.
    bputils::CoordsetSlicer<ExecSpace, CoordsetView> css(m_coordsetView);
    bputils::SliceData slice;
    slice.m_indicesView = selectedIdsView;
    css.execute(slice, n_coordset, n_coordset);
std::cout << "n_coordset.size=" << n_coordset["values/x"].dtype().number_of_elements() << std::endl;
n_coordset.print();
  }

  CoordsetView m_coordsetView;
};

}  // end namespace blueprint
}  // end namespace utilities
}  // end namespace mir
}  // end namespace axom

#endif
