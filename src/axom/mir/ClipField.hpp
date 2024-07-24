// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_CLIP_FIELD_HPP_
#define AXOM_MIR_CLIP_FIELD_HPP_

#include "axom/core.hpp"
#include "axom/mir/clipping/ClipCases.h"
#include "axom/mir/clipping/ClipTableManager.hpp"
#include "axom/mir/views/Shapes.hpp"
#include "axom/mir/FieldBlender.hpp"
#include "axom/mir/CoordsetBlender.hpp"
#include "axom/mir/FieldSlicer.hpp"
#include "axom/mir/blueprint_utilities.hpp"
#include "axom/mir/utilities.hpp" // for cpp2conduit

#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint_mesh_utils.hpp>

#include <map>
#include <string>

namespace axom
{
namespace mir
{
namespace clipping
{
namespace details
{

#if 0
// NOTE: Longer term, we should hide more RAJA functionality behind axom wrappers
//       so we can write serial versions for when RAJA is not enabled.
namespace axom
{
namespace operators
{
template <typename DataType>
using plus = RAJA::operators::plus<DataType>;
} // end namespace operators
} // end namespace axom

template <typename ExecSpace, typename ArrayViewType, typename OperatorType>
void exclusive_scan(ArrayViewType &input, ArrayViewType &output, OperatorType &op)
{
  assert(input.size() == output.size());
  using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
  RAJA::exclusive_scan<loop_policy>(RAJA::make_span(input.data(), input.size()),
                                    RAJA::make_span(output.data(), output.size()),
                                    op);
}

#endif

std::map<std::string, int>
shapeMap_NameValue(const conduit::Node &n_shape_map)
{
  std::map<std::string, int> sm;
  for(conduit::index_t i = 0; i < n_shape_map.number_of_children(); i++)
  {
    sm[n_shape_map[i].name()] = n_shape_map[i].to_int();
  }
  return sm;
}

std::map<int, std::string>
shapeMap_ValueName(const conduit::Node &n_shape_map)
{
  std::map<int, std::string> sm;
  for(conduit::index_t i = 0; i < n_shape_map.number_of_children(); i++)
  {
    sm[n_shape_map[i].to_int()] = n_shape_map[i].name();
  }
  return sm;
}

/**
 * \brief Given an "ST_index" (e.g. ST_TET from clipping definitions), return an appropriate Shape::id() value.
 *
 * \param st_index The value we want to translate into a Shape::id() value.
 *
 * \return The Shape::id() value that matches the st_index, or 0 if there is no match.
 */
template <typename IntegerType>
AXOM_HOST_DEVICE
int ST_Index_to_ShapeID(IntegerType st_index)
{
  int shapeID = 0;
  switch(st_index)
  {
  case ST_LIN: shapeID = views::LineShape<int>::id(); break;
  case ST_TRI: shapeID = views::TriShape<int>::id(); break;
  case ST_QUA: shapeID = views::QuadShape<int>::id(); break;
  case ST_TET: shapeID = views::TetShape<int>::id(); break;
  case ST_PYR: shapeID = views::PyramidShape<int>::id(); break;
  case ST_WDG: shapeID = views::WedgeShape<int>::id(); break;
  case ST_HEX: shapeID = views::HexShape<int>::id(); break;
  }
  return shapeID;
}

/**
 * \brief Given a flag that includes bitwise-or'd shape ids, make a map that indicates which Conduit shapes are used.
 *
 * \param shapes This is a bitwise-or of various Shape::id() values.
 *
 * \return A map of Conduit shape name to Shape::id() value.
 */
std::map<std::string, int>
shapeMap_FromFlags(std::uint64_t shapes)
{
  std::map<std::string, int> sm;

  if((shapes & views::LineShape<int>::id()) > 0)
    sm["line"] = views::LineShape<int>::id();

  if((shapes & views::TriShape<int>::id()) > 0)
    sm["tri"] = views::TriShape<int>::id();

  if((shapes & views::QuadShape<int>::id()) > 0)
    sm["quad"] = views::QuadShape<int>::id();

  if((shapes & views::TetShape<int>::id()) > 0)
    sm["tet"] = views::TetShape<int>::id();

  if((shapes & views::PyramidShape<int>::id()) > 0)
    sm["pyramid"] = views::PyramidShape<int>::id();

  if((shapes & views::WedgeShape<int>::id()) > 0)
    sm["wedge"] = views::WedgeShape<int>::id();

  if((shapes & views::HexShape<int>::id()) > 0)
    sm["hex"] = views::HexShape<int>::id();

  return sm;
}

template <typename IndexT>
AXOM_HOST_DEVICE
IndexT getClipTableIndex(IndexT shapeId)
{
  IndexT index = 0;
  switch(shapeId)
  {
  case views::TriShape<IndexT>::id():     index = 0; break;
  case views::QuadShape<IndexT>::id():    index = 1; break;
  case views::TetShape<IndexT>::id():     index = 2; break;
  case views::PyramidShape<IndexT>::id(): index = 3; break;
  case views::WedgeShape<IndexT>::id():   index = 4; break;
  case views::HexShape<IndexT>::id():     index = 5; break;
  }
  return index;
}

AXOM_HOST_DEVICE
bool color0Selected(int selection)
{
  return axom::utilities::bitIsSet(selection, 0);
}

AXOM_HOST_DEVICE
bool color1Selected(int selection)
{
  return axom::utilities::bitIsSet(selection, 1);
}

AXOM_HOST_DEVICE
bool generatedPointIsSelected(unsigned char color, int selection)
{
  return color == NOCOLOR ||
         (color0Selected(selection) && color == COLOR0) ||
         (color1Selected(selection) && color == COLOR1);
}

AXOM_HOST_DEVICE
bool shapeIsSelected(unsigned char color, int selection)
{
  return (color0Selected(selection) && color == COLOR0) ||
         (color1Selected(selection) && color == COLOR1);
}

AXOM_HOST_DEVICE
float computeWeight(float d0, float d1, float clipValue)
{
  const float delta = d1 - d0;
  const float abs_delta = (delta < 0) ? -delta : delta;
  const float t = (abs_delta != 0.) ? ((clipValue - d0) / delta) : 0.;
  return t;
}


// TODO: Could we make ZoneType be a concept?
/**
 * \brief Use the ids for the provided zone and the values from the array view to determine a clip case.
 *
 * \tparam ZoneType A class that implements the Shape interface.
 * \tparam DataType The data type of the clip data.
 *
 * \param[in] zone The zone being clipped.
 * \param[in] view The view that can access the clipping distance field.
 * \param[in] clipValue The value used for clipping.
 *
 * \return The index of the clipping case.
 */
template <typename ZoneType, typename DataType>
AXOM_HOST_DEVICE
size_t clip_case(const ZoneType &zone, const axom::ArrayView<DataType> &view, DataType clipValue)
{
  size_t clipcase = 0;
  for(size_t i = 0; i < zone.numberOfNodes(); i++)
  {
    const auto id = zone.getId(i);
    const auto value = view[id] - clipValue;
    clipcase |= (value > 0) ? (1 << i) : 0;
  }
  return clipcase;
}

} // end namespace details

/**
 * \brief This class provides a kind of schema over the clipping options, as well
 *        as default values, and some utilities functions.
 */
template <typename ExecSpace>
class ClipOptions
{
public:
  ClipOptions(axom::IndexType nzones, const conduit::Node &options) : m_nzones(nzones), m_options(options), m_selectedZones()
  {
  }

  /**
   * \brief Return a view that contains the list of selected zone ids for the mesh.
   * \return A view that contains the list of selected zone ids for the mesh.
   */
  axom::ArrayView<axom::IndexType> selectedZonesView()
  {
    if(m_selectedZones.size() == 0)
      buildSelectedZones();
    return m_selectedZones.view();
  }

  /**
   * \brief Invalidate the selected zones array (due to options changing) so we can rebuild it.
   */
  void invalidateSelectedZones()
  {
    m_selectedZones.clear();
  }

  /**
   * \brief Return the name of the field used for clipping.
   * \return The name of the field used for clipping.
   */
  std::string clipField() const
  {
    return m_options.fetch_existing("clipField").as_string();
  }

  /**
   * \brief Return the clip value.
   * \return The clip value.
   */
  float clipValue() const
  {
    return m_options.has_child("clipValue") ? m_options.fetch_existing("clipValue").to_float() : 0.f;
  }

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
   * \brief Return the name of the new color field to be created.
   * \return The name of the new color field to be created.
   */
  std::string colorField() const
  {
    std::string name("color");
    if(m_options.has_child("colorField"))
      name = m_options.fetch_existing("colorField").as_string();
    return name;
  }

  /**
   * \brief Whether the "inside" of the clipping field is selected.
   * \return 1 of the inside clipping is selected, false otherwise.
   */
  bool inside() const
  {
    return m_options.has_path("inside") ? (m_options.fetch_existing("inside").to_int() > 0) : true;
  }

  /**
   * \brief Whether the "outside" of the clipping field is selected.
   * \return 1 of the outside clipping is selected, false otherwise.
   */
  bool outside() const
  {
    return m_options.has_path("outside") ? (m_options.fetch_existing("outside").to_int() > 0) : false;
  }

  /**
   * \brief Extract the names of the fields to process (and their output names) from the
   *        options or \a n_fields if the options do not contain fields.
   *
   * \param n_fields The Conduit node that contains mesh fields.
   *
   * \return A map of the fields that will be processed, as well as their output name in the new fields.
   */
  std::map<std::string, std::string> fields(const conduit::Node &n_fields) const
  {
    std::map<std::string, std::string> f;
    if(m_options.has_child("fields"))
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
    else
    {
      // No options were specified. Allow all fields with same topology as clipField.
      const conduit::Node &n_clipField = n_fields.fetch_existing(clipField());
      std::string topoName = n_clipField.fetch_existing("topology").as_string();
      for(conduit::index_t i = 0; i < n_fields.number_of_children(); i++)
      {
        if(topoName == n_fields[i].fetch_existing("topology").as_string())
          f[n_fields[i].name()] = n_fields[i].name();
      }
    }
    return f;
  }  

private:
  /**
   * \brief The options may contain a "selectedZones" member that is a list of zones
   *        that will be operated on. If such an array is present, copy and sort it.
   *        If the zone list is not present, make an array that selects every zone.
   */
  void buildSelectedZones()
  {
    const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    if(m_options.has_child("selectedZones"))
    {
      // Store the zone list in m_selectedZones.
      int badValueCount = 0;
      views::IndexNode_to_ArrayView(m_options["selectedZones"], [&](auto zonesView)
      {
        using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
        using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;
        using value_type = typename decltype(zonesView)::value_type;

        // It probably does not make sense to request more zones than we have in the mesh.
        SLIC_ASSERT(zonesView.size() <= m_nzones);

        m_selectedZones = axom::Array<axom::IndexType>(zonesView.size(), zonesView.size(), allocatorID);
        auto szView = m_selectedZones.view();
        axom::copy(szView.data(), zonesView.data(), zonesView.size() * sizeof(value_type));

        // Check that the selected zone values are in range.
        const auto nzones = m_nzones;
        RAJA::ReduceSum<reduce_policy, int> errReduce(0);
        axom::for_all<ExecSpace>(szView.size(), AXOM_LAMBDA(auto index)
        {
          const int err = (szView[index] < 0 || szView[index] >= nzones) ? 1 : 0;
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
      m_selectedZones = axom::Array<axom::IndexType>(m_nzones, m_nzones, allocatorID);
      auto szView = m_selectedZones.view();
      axom::for_all<ExecSpace>(m_nzones, AXOM_LAMBDA(auto zoneIndex)
      {
        szView[zoneIndex] = zoneIndex;
      });
    }
  }

private:
  axom::IndexType m_nzones;                     // The number of zones in the associated topology.
  const conduit::Node &m_options;               // A reference to the clipping options node.
  axom::Array<axom::IndexType> m_selectedZones; // Storage for a list of selected zone ids.
};




/**
 * \accelerated
 * \brief This class clips a topology using a field and puts the new topology into a new Conduit node.
 *
 * \tparam ExecSpace    The execution space where the compute-heavy kernels run.
 * \tparam TopologyView The topology view that can operate on the Blueprint topology.
 * \tparam CoordsetView The coordset view that can operate on the Blueprint coordset.
 */
template <typename ExecSpace, typename TopologyView, typename CoordsetView>
class ClipField
{
public:
  using BlendData = axom::mir::utilities::blueprint::BlendData;
  using SliceData = axom::mir::utilities::blueprint::SliceData;
  using ClipTableViews = axom::StackArray<axom::mir::clipping::TableView, 6>;

  /**
   * \brief Constructor
   *
   * \param topoView A topology view suitable for the supplied topology.
   * \param coordsetView A coordset view suitable for the supplied coordset.
   *
   */
  ClipField(const TopologyView &topoView, const CoordsetView &coordsetView) : m_topologyView(topoView), m_coordsetView(coordsetView), m_clipTables()
  {
  }

  /**
   * \brief Execute the clipping operation using the data stored in the specified \a clipField.
   *
   * \param[in] n_input The Conduit node that contains the topology, coordsets, and fields.
   * \param[in] n_options A Conduit node that contains clipping options.
   * \param[out] n_output A Conduit node that will hold the clipped output mesh. This should be a different node from \a n_input.
   *
   * \note The clipField field must currently be vertex-associated.
   */
  void execute(const conduit::Node &n_input,
               const conduit::Node &n_options,
               conduit::Node &n_output)
  {
     ClipOptions<ExecSpace> opts(0, n_options);
     const std::string clipFieldName = opts.clipField();

     const conduit::Node &n_fields = n_input.fetch_existing("fields");
     const conduit::Node &n_clipField = n_fields.fetch_existing(clipFieldName);
     const std::string &topoName = n_clipField["topology"].as_string();
     const conduit::Node &n_topo = n_input.fetch_existing("topologies/" + topoName);
     const std::string &coordsetName = n_topo["coordset"].as_string();
     const conduit::Node &n_coordset = n_input.fetch_existing("coordsets/" + coordsetName);
     
     execute(n_topo, n_coordset, n_fields,
             n_options,
             n_output["topologies/" + opts.topologyName(topoName)],
             n_output["coordsets/" + opts.coordsetName(coordsetName)],
             n_output["fields"]);
  }

  /**
   * \brief Execute the clipping operation using the data stored in the specified \a clipField.
   *
   * \param[in] n_topo The node that contains the input mesh topology.
   * \param[in] n_coordset The node that contains the input mesh coordset.
   * \param[in] n_fields The node that contains the input fields.
   * \param[in] n_options A Conduit node that contains clipping options.
   * \param[out] n_newTopo A node that will contain the new clipped topology.
   * \param[out] n_newCoordset A node that will contain the new coordset for the clipped topology.
   * \param[out] n_newFields A node that will contain the new fields for the clipped topology.
   *
   * \note The clipField field must currently be vertex-associated. Also, the output topology will be an unstructured topology with mixed shape types.
   */
  void execute(const conduit::Node &n_topo,
               const conduit::Node &n_coordset,
               const conduit::Node &n_fields,
               const conduit::Node &n_options,
               conduit::Node &n_newTopo,
               conduit::Node &n_newCoordset,
               conduit::Node &n_newFields)
  {
    using KeyType = typename BlendData::KeyType;
    using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
    using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;
    const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    using ConnectivityType = typename TopologyView::IndexType;
    constexpr auto connTypeID = axom::mir::utilities::blueprint::cpp2conduit<ConnectivityType>::id;

    const auto nzones = m_topologyView.numberOfZones();
    ClipOptions<ExecSpace> opts(nzones, n_options);

    // Get the clip field.
    std::string clipFieldName = opts.clipField();
    const conduit::Node &n_clip_field = n_fields.fetch_existing(opts.clipField());
    const conduit::Node &n_clip_field_values = n_clip_field["values"];
    SLIC_ASSERT(n_clip_field["association"].as_string() == "vertex");

    // Determine which parts of the shapes will be kept.
    int selection = 0;
    if(opts.inside())
      axom::utilities::setBit(selection, 0);
    if(opts.outside())
      axom::utilities::setBit(selection, 1);
    SLIC_ASSERT(selection > 0);

    // Load clip table data and make views.
    m_clipTables.load(m_topologyView.dimension());
    ClipTableViews clipTableViews;
    createClipTableViews(clipTableViews, m_topologyView.dimension());

    // ----------------------------------------------------------------------
    //
    // Stage 1: Iterate over elements and their respective clip cases to
    //          determine sizes of outputs.
    //
    // ----------------------------------------------------------------------
    RAJA::ReduceSum<reduce_policy, int> fragment_sum(0);
    RAJA::ReduceSum<reduce_policy, int> fragment_nids_sum(0);
    RAJA::ReduceSum<reduce_policy, int> blendGroups_sum(0);
    RAJA::ReduceSum<reduce_policy, int> blendGroupLen_sum(0);

    // Allocate some memory
    axom::Array<int> clipCases(nzones, nzones, allocatorID);      // The clip case for a zone.
    axom::Array<std::uint64_t> pointsUsed(nzones, nzones, allocatorID);   // Which points are used over all selected fragments in a zone
    axom::Array<int> blendGroups(nzones, nzones, allocatorID);    // Number of blend groups in a zone.
    axom::Array<int> blendGroupsLen(nzones, nzones, allocatorID); // Length of the blend groups in a zone.
    axom::Array<int> fragments(nzones, nzones, allocatorID);      // The number of fragments (child zones) produced for a zone.
    axom::Array<int> fragmentsSize(nzones, nzones, allocatorID);  // The connectivity size for all selected fragments in a zone.

    auto clipCasesView = clipCases.view();
    auto pointsUsedView = pointsUsed.view();
    auto blendGroupsView = blendGroups.view();
    auto blendGroupsLenView = blendGroupsLen.view();
    auto fragmentsView = fragments.view();
    auto fragmentsSizeView = fragmentsSize.view();

    views::Node_to_ArrayView(n_clip_field_values, [&](auto clipFieldView)
    {
      using clip_value_type = typename decltype(clipFieldView)::value_type;
      const auto clipValue = static_cast<clip_value_type>(opts.clipValue());

      m_topologyView.template for_selected_zones<ExecSpace>(opts.selectedZonesView(), AXOM_LAMBDA(auto zoneIndex, const auto &zone)
      {
        // Get the clip case for the current zone.
        const auto clipcase = details::clip_case(zone, clipFieldView, clipValue);
        clipCasesView[zoneIndex] = clipcase;

        // Iterate over the shapes in this clip case to determine the number of blend groups.
        const auto clipTableIndex = details::getClipTableIndex(zone.id());
        const auto &ctView = clipTableViews[clipTableIndex];

        int thisBlendGroups = 0;    // The number of blend groups produced in this case.
        int thisBlendGroupLen = 0;  // The total length of the blend groups.
        int thisFragments = 0;      // The number of zone fragments produced in this case.
        int thisFragmentsNumIds = 0;// The number of points used to make all the fragment zones.
        std::uint64_t ptused = 0;   // A bitset indicating which ST_XX nodes are used.

        auto it = ctView.begin(clipcase);
        const auto end = ctView.end(clipcase);
        for(; it != end; it++)
        {
          // Get the current shape in the clip case.
          const auto caseData = *it;

          if(caseData[0] == ST_PNT)
          {
            if(details::generatedPointIsSelected(caseData[2], selection))
            {
              const size_t nIds = caseData[3];
              for(size_t ni = 4; ni < nIds; ni++)
              {
                const auto pid = caseData[ni];

                // Increase the blend size to include this center point.
                if(pid <= P7)
                {
                   // corner point.
                   thisBlendGroupLen++;
                }
                else if(pid >= EA && pid <= EL)
                {
                  // edge points are derived from 2 corner points. If
                  // those appear here then we're probably creating a
                  // face point. We can store the 2 corner points in place
                  // of the edge point (along with some blending coeff).
                  thisBlendGroupLen += 2;
                }
              }

              // This center or face point counts as a blend group.
              thisBlendGroups++;

              // Mark the point used.
              axom::utilities::setBit(ptused, N0 + caseData[1]);
            }
          }
          else
          {
            if(details::shapeIsSelected(caseData[1], selection))
            {
              thisFragments++;
              const int nIdsThisFragment = caseData.size() - 2;
              thisFragmentsNumIds += nIdsThisFragment;

              // Mark the points in this fragment used.
              for(int i = 2; i < caseData.size(); i++)
              {
                axom::utilities::setBit(ptused, caseData[i]);
              }
            }
          }
        }

        // Save the flags for the points that were used in this zone
        pointsUsedView[zoneIndex] = ptused;

        // Count which points in the original cell are used.
        for(unsigned char pid = P0; pid <= P7; pid++)
        {
          const int incr = axom::utilities::bitIsSet(ptused, pid) ? 1 : 0;
          thisBlendGroupLen += incr; // {p0}
          thisBlendGroups += incr;
        }

        // Count edges that are used.
        for(unsigned char pid = EA; pid <= EL; pid++)
        {
          const int incr = axom::utilities::bitIsSet(ptused, pid) ? 1 : 0;
          thisBlendGroupLen += 2 * incr; // {p0 p1}
          thisBlendGroups += incr;
        }

        // Save the results.
        blendGroupsView[zoneIndex] = thisBlendGroups;
        blendGroupsLenView[zoneIndex] = thisBlendGroupLen;
        fragmentsView[zoneIndex] = thisFragments;
        fragmentsSizeView[zoneIndex] = thisFragmentsNumIds;

        // Sum up the sizes overall.
        fragment_sum += thisFragments;
        fragment_nids_sum += thisFragmentsNumIds;
        blendGroups_sum += thisBlendGroups;
        blendGroupLen_sum += thisBlendGroupLen;
      });
    });


// TODO: change int to topoView::IndexType

    // ----------------------------------------------------------------------
    //
    // Stage 2: Do some scans to fill out blendOffset and blendGroupOffsets,
    //          which is where we fill in the real data.
    //
    // blendOffset : Starting offset for blending data like blendIds, blendCoeff.
    // blendGroupOffset : Starting offset for blendNames, blendGroupSizes.
    // fragmentOffsets : Where an element's fragments begin in the output.
    // ----------------------------------------------------------------------
    axom::Array<int> blendOffset(nzones, nzones, allocatorID);
    axom::Array<int> blendGroupOffsets(nzones, nzones, allocatorID);
    axom::Array<int> fragmentOffsets(nzones, nzones, allocatorID);
    axom::Array<int> fragmentSizeOffsets(nzones, nzones, allocatorID);

    auto blendOffsetView = blendOffset.view();
    auto blendGroupOffsetsView = blendGroupOffsets.view();
    auto fragmentOffsetsView = fragmentOffsets.view();
    auto fragmentSizeOffsetsView = fragmentSizeOffsets.view();

    // Make offsets via scan.
    RAJA::exclusive_scan<loop_policy>(RAJA::make_span(blendGroupsLenView.data(), nzones),
                                      RAJA::make_span(blendOffsetView.data(), nzones),
                                      RAJA::operators::plus<int>{});

    RAJA::exclusive_scan<loop_policy>(RAJA::make_span(blendGroupsView.data(), nzones),
                                      RAJA::make_span(blendGroupOffsetsView.data(), nzones),
                                      RAJA::operators::plus<int>{});

    RAJA::exclusive_scan<loop_policy>(RAJA::make_span(fragmentsView.data(), nzones),
                                      RAJA::make_span(fragmentOffsetsView.data(), nzones),
                                      RAJA::operators::plus<int>{});

    RAJA::exclusive_scan<loop_policy>(RAJA::make_span(fragmentsSizeView.data(), nzones),
                                      RAJA::make_span(fragmentSizeOffsetsView.data(), nzones),
                                      RAJA::operators::plus<int>{});

    // ----------------------------------------------------------------------
    //
    // Stage 3: Iterate over the elements/cases again and fill in the blend
    //          groups that get produced: blendNames, blendGroupSizes,
    //          blendCoeff, blendIds. These are used to produce the new points.
    //
    //          NOTE: blendGroupStart is a scan of blendGroupSizes.
    //
    // ----------------------------------------------------------------------
    const auto blendGroupsSize = blendGroups_sum.get();
    const auto blendGroupLenSize = blendGroupLen_sum.get();

    axom::Array<KeyType> blendNames(blendGroupsSize, blendGroupsSize, allocatorID);
    axom::Array<IndexType> blendGroupSizes(blendGroupsSize, blendGroupsSize, allocatorID);
    axom::Array<IndexType> blendGroupStart(blendGroupsSize, blendGroupsSize, allocatorID);
    axom::Array<IndexType> blendIds(blendGroupLenSize, blendGroupLenSize, allocatorID);
    axom::Array<float> blendCoeff(blendGroupLenSize, blendGroupLenSize, allocatorID);

    auto blendNamesView = blendNames.view();
    auto blendGroupSizesView = blendGroupSizes.view();
    auto blendGroupStartView = blendGroupStart.view();
    auto blendIdsView = blendIds.view();
    auto blendCoeffView = blendCoeff.view();

    views::Node_to_ArrayView(n_clip_field_values, [&](auto clipFieldView)
    {
      const auto clipValue = opts.clipValue();
      m_topologyView.template for_selected_zones<ExecSpace>(opts.selectedZonesView(), AXOM_LAMBDA(auto zoneIndex, const auto &zone)
      {
        using ZoneType = typename std::remove_reference<decltype(zone)>::type;

        // Get the clip case for the current zone.
        const auto clipcase = clipCasesView[zoneIndex];

        // Iterate over the shapes in this clip case to determine the number of blend groups.
        const auto clipTableIndex = details::getClipTableIndex(zone.id());
        const auto &ctView = clipTableViews[clipTableIndex];

        const std::uint64_t ptused = pointsUsedView[zoneIndex];

        // Starting offset of where we store this element's blend groups.
        std::uint32_t bgStart = blendOffsetView[zoneIndex];
        std::uint32_t bgOffset = blendGroupOffsetsView[zoneIndex];

        auto it = ctView.begin(clipcase);
        const auto end = ctView.end(clipcase);
        for(; it != end; it++)
        {
          // Get the current shape in the clip case.
          const auto caseData = *it;

          if(caseData[0] == ST_PNT)
          {
            if(details::generatedPointIsSelected(caseData[2], selection))
            {
// TODO: put the ability to select on color back in... 0, 1, or both
              const int nIds = static_cast<int>(caseData[3]);
              const auto one_over_n = 1.f / static_cast<float>(nIds);
              const auto start = bgStart;

              for(int ni = 0; ni < nIds; ni++)
              {
                const auto ptid = caseData[4 + ni];

                // Add the point to the blend group.
                if(ptid <= P7)
                {
                  // corner point.
                  blendIdsView[bgStart] = zone.getId(ptid);
                  blendCoeffView[bgStart] = one_over_n;

                  bgStart++;
                }
                else if(ptid >= EA && ptid <= EL)
                {
                  // edge points are derived from 2 corner points. If
                  // those appear here then we're probably creating a
                  // face point. We can store the 2 corner points in place
                  // of the edge point (along with some blending coeff).
                  const auto edgeIndex = ptid - EA;
                  const auto edge = ZoneType::getEdge(edgeIndex);
                  const auto id0 = zone.getId(edge[0]);
                  const auto id1 = zone.getId(edge[1]);

                  // Figure out the blend for edge.
                  const float t = details::computeWeight(clipFieldView[id0], clipFieldView[id1], clipValue);
                
                  blendIdsView[bgStart]   = id0;
                  blendIdsView[bgStart+1] = id1;
                  blendCoeffView[bgStart] = one_over_n * (1. - t);
                  blendCoeffView[bgStart+1] = one_over_n * t;

                  bgStart += 2;
                }
              }

              // Store how many points make up this blend group. Note that the
              // size will not necessarily be equal to npts if edges were involved.
              std::uint32_t nblended = bgStart - blendGroupStartView[bgOffset];
              blendGroupSizesView[bgOffset] = nblended;

              // Store "name" of blend group.
              const auto blendName = axom::mir::utilities::make_name_n(blendIdsView.data() + start, nblended);
              blendNamesView[bgOffset++] = blendName;
            }
          }
        }

        // Add blend group for each original point that was used.
        for(unsigned char pid = P0; pid <= P7; pid++)
        {
          if(axom::utilities::bitIsSet(ptused, pid))
          {
            // Store blend group info
            blendIdsView[bgStart] = zone.getId(pid);
            blendCoeffView[bgStart] = 1.;

            // Store how many points make up this blend group.
            blendGroupSizesView[bgOffset] = 1;

            // Store where this blendGroup starts in the blendIds,blendCoeff.
            blendGroupStartView[bgOffset] = bgStart;

            // Store "name" of blend group.
            blendNamesView[bgOffset++] = axom::mir::utilities::make_name_1(zone.getId(pid));

            bgStart++;
          }
        }

        // Add blend group for each edge point that was used.
        for(unsigned char pid = EA; pid <= EL; pid++)
        {
          if(axom::utilities::bitIsSet(ptused, pid))
          {
            const auto edgeIndex = pid - EA;
            const auto edge = ZoneType::getEdge(edgeIndex);
            const auto id0 = zone.getId(edge[0]);
            const auto id1 = zone.getId(edge[1]);

            // Figure out the blend for edge.
            const float t = details::computeWeight(clipFieldView[id0], clipFieldView[id1], clipValue);

            // Store blend group info                
            blendIdsView[bgStart]   = id0;
            blendIdsView[bgStart+1] = id1;
            blendCoeffView[bgStart] = (1. - t);
            blendCoeffView[bgStart+1] = t;

            // Store how many points make up this blend group.
            blendGroupSizesView[bgOffset] = 2;

            // Store where this blendGroup starts in the blendIds,blendCoeff.
            blendGroupStartView[bgOffset] = bgStart;

            // Store "name" of blend group.
            blendNamesView[bgOffset++] = axom::mir::utilities::make_name_2(id0, id1);

            bgStart += 2;
          }
        }
      });
    });

    // ----------------------------------------------------------------------
    //
    // Stage 4 - Make the blend groups unique based on their blendName.
    //
    // ----------------------------------------------------------------------
    // At this point, we have created the blend group data. We can now use the
    // blendNames to make unique blend groups. uNames contains a sorted list of
    // the unique blend group names while uIndices is their original index in
    // blendNames/blendGroupOffsets/blendGroupSizes.
    axom::Array<KeyType> uNames;
    axom::Array<axom::IndexType> uIndices;
    axom::mir::utilities::unique<ExecSpace, KeyType>(blendNames, uNames, uIndices);

    auto uNamesView = uNames.view();

    // Bundle up blend data.
    axom::mir::utilities::blueprint::BlendData blend;
    blend.m_blendNamesView = uNames.view();
    blend.m_uniqueIndicesView = uIndices.view();
    blend.m_blendGroupSizesView = blendGroupSizes.view();
    blend.m_blendGroupStartView = blendGroupStart.view();
    blend.m_blendIdsView = blendIds.view();
    blend.m_blendCoeffView = blendCoeff.view();

    // ----------------------------------------------------------------------
    //
    // Stage 5 - Make new connectivity.
    //
    // ----------------------------------------------------------------------
    const auto finalNumZones = fragment_sum.get();
    const auto finalConnSize = fragment_nids_sum.get();

    n_newTopo.reset();
    n_newTopo["type"] = "unstructured";
    n_newTopo["coordset"] = n_newCoordset.name();

    // Get the ID of a Conduit allocator that will allocate through Axom with device allocator allocatorID.
    utilities::blueprint::ConduitAllocateThroughAxom c2a(allocatorID);
    const int conduitAllocatorID = c2a.getConduitAllocatorID();

    // Allocate connectivity.
    conduit::Node &n_conn = n_newTopo["elements/connectivity"];
    n_conn.set_allocator(conduitAllocatorID);
    n_conn.set(conduit::DataType(connTypeID, finalConnSize));
    auto connView = axom::ArrayView<ConnectivityType>(static_cast<ConnectivityType *>(n_conn.data_ptr()), finalConnSize);

    // Allocate shapes.
    conduit::Node &n_shapes = n_newTopo["elements/shapes"];
    n_shapes.set_allocator(conduitAllocatorID);
    n_shapes.set(conduit::DataType(connTypeID, finalNumZones));
    auto shapesView = axom::ArrayView<ConnectivityType>(static_cast<ConnectivityType *>(n_shapes.data_ptr()), finalNumZones);

    // Allocate sizes.
    conduit::Node &n_sizes = n_newTopo["elements/sizes"];
    n_sizes.set_allocator(conduitAllocatorID);
    n_sizes.set(conduit::DataType(connTypeID, finalNumZones));
    auto sizesView = axom::ArrayView<ConnectivityType>(static_cast<ConnectivityType *>(n_sizes.data_ptr()), finalNumZones);

    // Allocate offsets.
    conduit::Node &n_offsets = n_newTopo["elements/offsets"];
    n_offsets.set_allocator(conduitAllocatorID);
    n_offsets.set(conduit::DataType(connTypeID, finalNumZones));
    auto offsetsView = axom::ArrayView<ConnectivityType>(static_cast<ConnectivityType *>(n_offsets.data_ptr()), finalNumZones);

    // Allocate a color variable to keep track of the "color" of the fragments.
    conduit::Node &n_color = n_newFields[opts.colorField()];
    n_color["topology"] = opts.topologyName(n_newTopo.name());
    n_color["association"] = "element";
    conduit::Node &n_color_values = n_color["values"];
    n_color_values.set_allocator(conduitAllocatorID);
    n_color_values.set(conduit::DataType::int32(finalNumZones));
    auto colorView = axom::ArrayView<int>(static_cast<ConnectivityType *>(n_color_values.data_ptr()), finalNumZones);

    // Here we fill in the new connectivity, sizes, shapes.
    // We get the node ids from the unique blend names, de-duplicating points when making the new connectivity.
    RAJA::ReduceBitOr<reduce_policy, std::uint64_t> shapesUsed_reduce(0);
    m_topologyView.template for_selected_zones<ExecSpace>(opts.selectedZonesView(), AXOM_LAMBDA(auto zoneIndex, const auto &zone)
    {
      // If there are no fragments, return from lambda.
      if(fragmentsView[zoneIndex] == 0)
        return;

      // Seek to the start of the blend groups for this zone.
      std::uint32_t bgStart = blendGroupOffsetsView[zoneIndex];

      // Go through the points in the order they would have been added as blend
      // groups, get their blendName, and then overall index of that blendName
      // in uNames, the unique list of new dof names. That will be their index
      // in the final points.
      const std::uint64_t ptused = pointsUsedView[zoneIndex];
      ConnectivityType point_2_new[N3 + 1];
      for(unsigned char pid = N0; pid <= N3; pid++)
      {
        if(axom::utilities::bitIsSet(ptused, pid))
        {
          const auto name = blendNamesView[bgStart++];
          point_2_new[pid] = axom::mir::utilities::bsearch(name, uNamesView);
        }
      }
      for(unsigned char pid = P0; pid <= P7; pid++)
      {
        if(axom::utilities::bitIsSet(ptused, pid))
        {
          const auto name = blendNamesView[bgStart++];
          point_2_new[pid] = axom::mir::utilities::bsearch(name, uNamesView);
        }
      }
      for(unsigned char pid = EA; pid <= EL; pid++)
      {
        if(axom::utilities::bitIsSet(ptused, pid))
        {
          const auto name = blendNamesView[bgStart++];
          point_2_new[pid] = axom::mir::utilities::bsearch(name, uNamesView);
        }
      }

      // This is where the output fragment connectivity start for this zone
      int outputIndex = fragmentSizeOffsetsView[zoneIndex];
      // This is where the output fragment sizes/shapes start for this zone.
      int sizeIndex = fragmentOffsetsView[zoneIndex];
      std::uint64_t shapesUsed = 0;

      // Iterate over the selected fragments and emit connectivity for them.
      const auto clipcase = clipCasesView[zoneIndex];
      const auto clipTableIndex = details::getClipTableIndex(zone.id());
      const auto ctView = clipTableViews[clipTableIndex];
      auto it = ctView.begin(clipcase);
      const auto end = ctView.end(clipcase);
      for(; it != end; it++)
      {
        // Get the current shape in the clip case.
        const auto caseData = *it;
        const auto fragmentShape = static_cast<ConnectivityType>(caseData[0]);

        if(fragmentShape != ST_PNT)
        {
          if(details::shapeIsSelected(caseData[1], selection))
          {
            // Output the nodes used in this zone.
            for(int i = 2; i < caseData.size(); i++)
              connView[outputIndex++] = point_2_new[caseData[i]];

            const auto nIdsInFragment = caseData.size() - 2;
            const auto shapeID = details::ST_Index_to_ShapeID(fragmentShape);
            sizesView[sizeIndex] = nIdsInFragment;
            shapesView[sizeIndex] = shapeID;
            colorView[sizeIndex] = caseData[1] - COLOR0;

            sizeIndex++;

            // Record which shape type was used. (ids are powers of 2)
            shapesUsed |= shapeID;
          }
        }
      }

      shapesUsed_reduce |= shapesUsed;
    });

    // If inside and outside are not selected, remove the color field since we should not need it.
    if(!(opts.inside() && opts.outside()))
    {
      n_newFields.remove(opts.colorField());
    }

    // Make offsets
    RAJA::exclusive_scan<loop_policy>(RAJA::make_span(sizesView.data(), finalNumZones),
                                      RAJA::make_span(offsetsView.data(), finalNumZones),
                                      RAJA::operators::plus<ConnectivityType>{});

    // Add shape information to the connectivity.
    const auto shapesUsed = shapesUsed_reduce.get();
    const auto shapeMap = details::shapeMap_FromFlags(shapesUsed);
    if(axom::utilities::countBits(shapesUsed) > 1)
    {
      n_newTopo["elements/shape"] = "mixed";
      conduit::Node &n_shape_map = n_newTopo["elements/shape_map"];
      for(auto it = shapeMap.cbegin(); it != shapeMap.cend(); it++)
         n_shape_map[it->first] = it->second;
    }
    else
    {
      n_shapes.reset();
      n_newTopo["elements"].remove("shapes");
      n_newTopo["elements/shape"] = shapeMap.begin()->first;
    }

    //-----------------------------------------------------------------------
    // STAGE 6 - Make new coordset.
    //-----------------------------------------------------------------------   
    makeCoordset(blend, n_coordset, n_newCoordset);

    //-----------------------------------------------------------------------
    // STAGE 7 - Make new fields.
    //-----------------------------------------------------------------------
    axom::mir::utilities::blueprint::SliceData slice;
    axom::Array<IndexType> sliceIndices(finalNumZones, finalNumZones, allocatorID);
    auto sliceIndicesView = sliceIndices.view();
    axom::for_all<ExecSpace>(nzones, AXOM_LAMBDA(auto zoneIndex)
    {
      const auto start = fragmentOffsetsView[zoneIndex];
      for(int i = 0; i < fragmentsView[zoneIndex]; i++)
        sliceIndicesView[start + i] = zoneIndex;
    });
    slice.m_indicesView = sliceIndicesView;
    makeFields(blend, slice, opts.topologyName(n_topo.name()), opts.fields(n_fields), n_fields, n_newFields);

    //-----------------------------------------------------------------------
    // STAGE 8 - make originalElements (this will later be optional)
    //-----------------------------------------------------------------------
    if(n_fields.has_child("originalElements"))
    {
      // originalElements already exists. We need to map it forward.
      const conduit::Node &n_orig = n_fields["originalElements"];
      const conduit::Node &n_orig_values = n_orig["values"];
      views::IndexNode_to_ArrayView(n_orig_values, [&](auto origValuesView)
      {
        using value_type = typename decltype(origValuesView)::value_type;
        conduit::Node &n_origElem = n_newFields["originalElements"];       
        n_origElem["association"] = "element";
        n_origElem["topology"] = opts.topologyName(n_newTopo.name());
        conduit::Node &n_values = n_origElem["values"];
        n_values.set_allocator(conduitAllocatorID);
        n_values.set(conduit::DataType(n_orig_values.dtype().id(), finalNumZones));
        auto valuesView = axom::ArrayView<value_type>(static_cast<value_type *>(n_values.data_ptr()), finalNumZones);
        axom::for_all<ExecSpace>(nzones, AXOM_LAMBDA(auto zoneIndex)
        {
          int sizeIndex = fragmentOffsetsView[zoneIndex];
          int nFragments = fragmentsView[zoneIndex];
          for(int i = 0; i < nFragments; i++)
            valuesView[sizeIndex + i] = origValuesView[zoneIndex];
        });
      });
    }
    else
    {
      // Make a new node and populate originalElement.
      conduit::Node &n_orig = n_newFields["originalElements"];
      n_orig["association"] = "element";
      n_orig["topology"] = opts.topologyName(n_newTopo.name());
      conduit::Node &n_values = n_orig["values"];
      n_values.set_allocator(conduitAllocatorID);
      n_values.set(conduit::DataType(connTypeID, finalNumZones));
      auto valuesView = axom::ArrayView<ConnectivityType>(static_cast<ConnectivityType *>(n_values.data_ptr()), finalNumZones);
      axom::for_all<ExecSpace>(nzones, AXOM_LAMBDA(auto zoneIndex)
      {
        int sizeIndex = fragmentOffsetsView[zoneIndex];
        int nFragments = fragmentsView[zoneIndex];
        for(int i = 0; i < nFragments; i++)
          valuesView[sizeIndex + i] = zoneIndex;
      });
    }

//-----------------------------------------------------------------------------------
  } // end of execute

private:
  /**
   * \brief Create views for the clip tables of various shapes.
   *
   * \param[out] views The views array that will contain the table views.
   * \param dimension The dimension the topology (so we can load a subset of tables)
   */
  void createClipTableViews(ClipTableViews &views, int dimension)
  {
    if(dimension == -1 || dimension == 2)
    {
      views[details::getClipTableIndex(axom::mir::views::TriShape<IndexType>::id())] = m_clipTables[ST_TRI].view();
      views[details::getClipTableIndex(axom::mir::views::QuadShape<IndexType>::id())] = m_clipTables[ST_QUA].view();
    }
    if(dimension == -1 || dimension == 3)
    {
      views[details::getClipTableIndex(axom::mir::views::TetShape<IndexType>::id())] = m_clipTables[ST_TET].view();
      views[details::getClipTableIndex(axom::mir::views::PyramidShape<IndexType>::id())] = m_clipTables[ST_PYR].view();
      views[details::getClipTableIndex(axom::mir::views::WedgeShape<IndexType>::id())] = m_clipTables[ST_WDG].view();
      views[details::getClipTableIndex(axom::mir::views::HexShape<IndexType>::id())] = m_clipTables[ST_HEX].view();
    }
  }

  /**
   * \brief Make the new coordset using the blend data and the input coordset/coordsetview.
   *
   * \param blend The BlendData that we need to construct the new coordset.
   * \param n_coordset The input coordset, which is passed for metadata.
   * \param[out] n_newCoordset The new coordset.
   */
  void makeCoordset(const BlendData &blend, const conduit::Node &n_coordset, conduit::Node &n_newCoordset) const
  {
    axom::mir::utilities::blueprint::CoordsetBlender<ExecSpace, CoordsetView> cb;
    n_newCoordset.reset();
    cb.execute(blend, m_coordsetView, n_coordset, n_newCoordset);
  }

  /**
   * \brief Make new fields for the output topology.
   *
   * \param blend The BlendData that we need to construct the new vertex fields.
   * \param slice The SliceData we need to construct new element fields.
   * \param topologyName The name of the new field's topology.
   * \param fieldMap A map containing the names of the fields that we'll operate on.
   * \param n_fields The source fields.
   * \param[out] n_out_fields The node that will contain the new fields.
   */
  void makeFields(const BlendData &blend,
                  const SliceData &slice,
                  const std::string &topologyName,
                  const std::map<std::string, std::string> &fieldMap,
                  const conduit::Node &n_fields,
                  conduit::Node &n_out_fields) const
  {
    for(auto it = fieldMap.begin(); it != fieldMap.end(); it++)
    {
      const conduit::Node &n_field = n_fields.fetch_existing(it->first);
      const std::string association = n_field["association"].as_string();
      if(association == "element")
      {
        axom::mir::utilities::blueprint::FieldSlicer<ExecSpace> s;
        s.execute(slice, n_field, n_out_fields[it->second]);
        n_out_fields[it->second]["topology"] = topologyName;
      }
      else if(association == "vertex")
      {
        axom::mir::utilities::blueprint::FieldBlender<ExecSpace> b;
        b.execute(blend, n_field, n_out_fields[it->second]);
        n_out_fields[it->second]["topology"] = topologyName;
      }
    }
  }

private:
  TopologyView m_topologyView;
  CoordsetView m_coordsetView;
  axom::mir::clipping::ClipTableManager<ExecSpace> m_clipTables;
};

} // end namespace clipping
} // end namespace mir
} // end namespace axom

#endif
