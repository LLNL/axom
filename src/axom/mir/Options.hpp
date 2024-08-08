// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_OPTIONS_HPP_
#define AXOM_MIR_OPTIONS_HPP_

#include "axom/core.hpp"

#include <conduit/conduit.hpp>

namespace axom
{
namespace mir
{

/**
 * \brief This class provides a kind of schema over options, as well
 *        as default values, and some utilities functions.
 */
template <typename ExecSpace>
class Options
{
public:
  /**
   * \brief Constructor
   *
   * \param nzones The total number of zones in the associated topology.
   * \param options The node that contains the clipping options.
   */
  Options(axom::IndexType nzones, const conduit::Node &options)
    : m_nzones(nzones)
    , m_options(options)
    , m_selectedZones()
  { }

  /**
   * \brief Return a view that contains the list of selected zone ids for the mesh.
   * \return A view that contains the list of selected zone ids for the mesh.
   *
   * \note The data for the view is generated if it has not yet been built.
   */
  axom::ArrayView<axom::IndexType> selectedZonesView()
  {
    if(m_selectedZones.size() == 0)
    {
      buildSelectedZones();
    }
    return m_selectedZones.view();
  }

  /**
   * \brief Invalidate the selected zones array (due to options changing) so we can rebuild it.
   */
  void invalidateSelectedZones() { m_selectedZones.clear(); }

  /**
   * \brief Return the name of the new topology to be created.
   * \param default_value The name to use if the option is not defined.
   * \return The name of the new topology to be created.
   */
  std::string topologyName(const std::string &default_value = std::string()) const
  {
    std::string name(default_value);
    if(m_options.has_child("topologyName"))
      name = m_options.fetch_existing("topologyName").as_string();
    return name;
  }

  /**
   * \brief Return the name of the new coordset to be created.
   * \param default_value The name to use if the option is not defined.
   * \return The name of the new coordset to be created.
   */
  std::string coordsetName(const std::string &default_value = std::string()) const
  {
    std::string name(default_value);
    if(m_options.has_child("coordsetName"))
      name = m_options.fetch_existing("coordsetName").as_string();
    return name;
  }

  /**
   * \brief Extract the names of the fields to process (and their output names) from the
   *        options or \a n_fields if the options do not contain fields.
   *
   * \param[out] f A map of the fields that will be processed, as well as their output name in the new fields.
   * \return True if the fields were present in the options. False otherwise.
   */
  bool fields(std::map<std::string, std::string> &f) const
  {
    bool retval = m_options.has_child("fields");
    f.clear();
    if(retval)
    {
      const conduit::Node &n_opt_fields = m_options.fetch_existing("fields");
      for(conduit::index_t i = 0; i < n_opt_fields.number_of_children(); i++)
      {
        if(n_opt_fields[i].dtype().is_string())
          f[n_opt_fields[i].name()] = n_opt_fields[i].as_string();
        else
          f[n_opt_fields[i].name()] = n_opt_fields[i].name();
      }
    }
    return retval;
  }

protected:
  /**
   * \brief The options may contain a "selectedZones" member that is a list of zones
   *        that will be operated on. If such an array is present, copy and sort it.
   *        If the zone list is not present, make an array that selects every zone.
   *
   * \note selectedZones should contain local zone numbers, which in the case of
   *       strided-structured indexing are the [0..n) zone numbers that exist only
   *       within the selected window.
   */
  void buildSelectedZones()
  {
    const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    if(m_options.has_child("selectedZones"))
    {
      // Store the zone list in m_selectedZones.
      int badValueCount = 0;
      views::IndexNode_to_ArrayView(m_options["selectedZones"], [&](auto zonesView) {
        using loop_policy =
          typename axom::execution_space<ExecSpace>::loop_policy;
        using reduce_policy =
          typename axom::execution_space<ExecSpace>::reduce_policy;

        // It probably does not make sense to request more zones than we have in the mesh.
        SLIC_ASSERT(zonesView.size() <= m_nzones);

        m_selectedZones = axom::Array<axom::IndexType>(zonesView.size(),
                                                       zonesView.size(),
                                                       allocatorID);
        auto szView = m_selectedZones.view();
        axom::for_all<ExecSpace>(
          szView.size(),
          AXOM_LAMBDA(auto index) { szView[index] = zonesView[index]; });

        // Check that the selected zone values are in range.
        const auto nzones = m_nzones;
        RAJA::ReduceSum<reduce_policy, int> errReduce(0);
        axom::for_all<ExecSpace>(
          szView.size(),
          AXOM_LAMBDA(auto index) {
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
      m_selectedZones =
        axom::Array<axom::IndexType>(m_nzones, m_nzones, allocatorID);
      auto szView = m_selectedZones.view();
      axom::for_all<ExecSpace>(
        m_nzones,
        AXOM_LAMBDA(auto zoneIndex) { szView[zoneIndex] = zoneIndex; });
    }
  }

protected:
  axom::IndexType m_nzones;  // The number of zones in the associated topology.
  const conduit::Node &m_options;  // A reference to the clipping options node.
  axom::Array<axom::IndexType> m_selectedZones;  // Storage for a list of selected zone ids.
};
}  // end namespace mir
}  // end namespace axom

#endif
