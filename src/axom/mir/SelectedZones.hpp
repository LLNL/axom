// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_SELECTED_ZONES_HPP_
#define AXOM_MIR_SELECTED_ZONES_HPP_

#include "axom/core.hpp"

#include <conduit/conduit.hpp>

namespace axom
{
namespace mir
{
/*!
 * \brief This class provides a kind of schema over options, as well
 *        as default values, and some utilities functions.
 */
template <typename ExecSpace>
class SelectedZones
{
public:
  /*!
   * \brief Constructor
   *
   * \param nzones The total number of zones in the associated topology.
   * \param options The node that contains the clipping options.
   */
  SelectedZones(axom::IndexType nzones, const conduit::Node &options)
    : m_selectedZones()
    , m_selectedZonesView()
  {
    buildSelectedZones(nzones, options);
  }

  /*!
   * \brief Return a view that contains the list of selected zone ids for the mesh.
   * \return A view that contains the list of selected zone ids for the mesh.
   */
  const axom::ArrayView<axom::IndexType> &view() const
  {
    return m_selectedZonesView;
  }

protected:
  /*!
   * \brief The options may contain a "selectedZones" member that is a list of zones
   *        that will be operated on. If such an array is present, copy and sort it.
   *        If the zone list is not present, make an array that selects every zone.
   *
   * \note selectedZones should contain local zone numbers, which in the case of
   *       strided-structured indexing are the [0..n) zone numbers that exist only
   *       within the selected window.
   */
  void buildSelectedZones(axom::IndexType nzones, const conduit::Node &options)
  {
    const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    if(options.has_child("selectedZones"))
    {
      // Store the zone list in m_selectedZones.
      int badValueCount = 0;
      views::IndexNode_to_ArrayView(options["selectedZones"], [&](auto zonesView) {
        using loop_policy =
          typename axom::execution_space<ExecSpace>::loop_policy;
        using reduce_policy =
          typename axom::execution_space<ExecSpace>::reduce_policy;

        // It probably does not make sense to request more zones than we have in the mesh.
        SLIC_ASSERT(zonesView.size() <= nzones);

        m_selectedZones = axom::Array<axom::IndexType>(zonesView.size(),
                                                       zonesView.size(),
                                                       allocatorID);
        auto szView = m_selectedZonesView = m_selectedZones.view();
        axom::for_all<ExecSpace>(
          szView.size(),
          AXOM_LAMBDA(axom::IndexType index) { szView[index] = zonesView[index]; });

        // Check that the selected zone values are in range.
        RAJA::ReduceSum<reduce_policy, int> errReduce(0);
        axom::for_all<ExecSpace>(
          szView.size(),
          AXOM_LAMBDA(axom::IndexType index) {
            const int err =
              (szView[index] < 0 || szView[index] >= nzones) ? 1 : 0;
            errReduce += err;
          });
        badValueCount = errReduce.get();

        // Make sure the selectedZones are sorted.
        RAJA::sort<loop_policy>(RAJA::make_span(szView.data(), szView.size()));
      });

      if(badValueCount > 0)
      {
        SLIC_ERROR("Out of range selectedZones values.");
      }
    }
    else
    {
      // Select all zones.
      m_selectedZones = axom::Array<axom::IndexType>(nzones, nzones, allocatorID);
      auto szView = m_selectedZonesView = m_selectedZones.view();
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType zoneIndex) { szView[zoneIndex] = zoneIndex; });
    }
  }

protected:
  axom::Array<axom::IndexType> m_selectedZones;  // Storage for a list of selected zone ids.
  axom::ArrayView<axom::IndexType> m_selectedZonesView;
};

}  // end namespace mir
}  // end namespace axom

#endif
