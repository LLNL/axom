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
#include "axom/mir/utilities.hpp"
#include "axom/mir/views/view_traits.hpp"
#include "axom/mir/ClipOptions.hpp"
#include "axom/mir/BlendGroupBuilder.hpp"
#include "axom/mir/utilities.hpp"
#include "axom/mir/SelectedZones.hpp"
#include "axom/slic.hpp"

#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint_mesh_utils.hpp>
#include <conduit/conduit_relay_io.hpp>

#include <map>
#include <string>

//#define AXOM_DEBUG_CLIP_FIELD
#define AXOM_CLIP_FILTER_DEGENERATES

namespace axom
{
namespace mir
{
namespace clipping
{
namespace details
{
/**
 * \brief Given an "ST_index" (e.g. ST_TET from clipping definitions), return an appropriate ShapeID value.
 *
 * \param st_index The value we want to translate into a ShapeID value.
 *
 * \return The ShapeID value that matches the st_index, or 0 if there is no match.
 */
template <typename IntegerType>
inline AXOM_HOST_DEVICE IntegerType ST_Index_to_ShapeID(IntegerType st_index)
{
  IntegerType shapeID = 0;
  switch(st_index)
  {
  case ST_LIN:
    shapeID = views::Line_ShapeID;
    break;
  case ST_TRI:
    shapeID = views::Tri_ShapeID;
    break;
  case ST_QUA:
    shapeID = views::Quad_ShapeID;
    break;
  case ST_TET:
    shapeID = views::Tet_ShapeID;
    break;
  case ST_PYR:
    shapeID = views::Pyramid_ShapeID;
    break;
  case ST_WDG:
    shapeID = views::Wedge_ShapeID;
    break;
  case ST_HEX:
    shapeID = views::Hex_ShapeID;
    break;
  }
  return shapeID;
}

/**
 * \brief Returns a clip table index for the input shapeId.
 * \param shapeId A shapeID (e.g. Tet_ShapeID)
 * \return The clip table index for the shape.
 */
AXOM_HOST_DEVICE
inline int getClipTableIndex(int shapeId)
{
  int index = 0;
  switch(shapeId)
  {
  case views::Tri_ShapeID:
    index = 0;
    break;
  case views::Quad_ShapeID:
    index = 1;
    break;
  case views::Tet_ShapeID:
    index = 2;
    break;
  case views::Pyramid_ShapeID:
    index = 3;
    break;
  case views::Wedge_ShapeID:
    index = 4;
    break;
  case views::Hex_ShapeID:
    index = 5;
    break;
  }
  return index;
}

AXOM_HOST_DEVICE
inline bool color0Selected(int selection)
{
  return axom::utilities::bitIsSet(selection, 0);
}

AXOM_HOST_DEVICE
inline bool color1Selected(int selection)
{
  return axom::utilities::bitIsSet(selection, 1);
}

AXOM_HOST_DEVICE
inline bool generatedPointIsSelected(unsigned char color, int selection)
{
  return color == NOCOLOR || (color0Selected(selection) && color == COLOR0) ||
    (color1Selected(selection) && color == COLOR1);
}

AXOM_HOST_DEVICE
inline bool shapeIsSelected(unsigned char color, int selection)
{
  return (color0Selected(selection) && color == COLOR0) ||
    (color1Selected(selection) && color == COLOR1);
}

template <typename IdType, int MAXSIZE>
inline AXOM_HOST_DEVICE int unique_count(const IdType *values, int n)
{
  IdType v[MAXSIZE];
  // Start off with one unique element
  int nv = 1;
  v[0] = values[0];
  // Scan the rest
  for(int j = 1; j < n; j++)
  {
    int fi = -1;
    for(int i = 0; i < nv; i++)
    {
      if(values[j] == v[i])
      {
        fi = i;
        break;
      }
    }
    if(fi == -1)
    {
      v[nv++] = values[j];
    }
  }
  return nv;
}

#if defined(AXOM_DEBUG_CLIP_FIELD)
template <typename ViewType>
void printHost(const std::string &name, const ViewType &deviceView)
{
  using value_type = typename ViewType::value_type;
  int nn = deviceView.size();
  value_type *host = new value_type[nn];
  axom::copy(host, deviceView.data(), sizeof(value_type) * nn);
  std::cout << name << "[" << nn << "] = {";
  for(int ii = 0; ii < nn; ii++)
  {
    std::cout << ", " << host[ii];
  }
  std::cout << "}" << std::endl;
  delete[] host;
}
#endif

}  // end namespace details

//------------------------------------------------------------------------------
/**
 * \brief This class helps ClipField determine intersection cases and weights
 *        using a field designated by the options.
 */
template <typename ExecSpace, typename ConnectivityType>
class FieldIntersector
{
public:
  using ClipFieldType = float;
  using ConnectivityView = axom::ArrayView<ConnectivityType>;

  /**
   * \brief This is a view class for FieldIntersector that can be used in device code.
   */
  struct View
  {
    /**
     * \brief Given a zone index and the node ids that comprise the zone, return
     *        the appropriate clip case, taking into account the clip field and
     *        clip value.
     *
     * \param zoneIndex The zone index.
     * \param nodeIds A view containing node ids for the zone.
     */
    AXOM_HOST_DEVICE
    axom::IndexType determineClipCase(axom::IndexType AXOM_UNUSED_PARAM(zoneIndex),
                                      const ConnectivityView &nodeIds) const
    {
      axom::IndexType clipcase = 0;
      for(IndexType i = 0; i < nodeIds.size(); i++)
      {
        const auto id = nodeIds[i];
        const auto value = m_clipFieldView[id] - m_clipValue;
        clipcase |= (value > 0) ? (1 << i) : 0;
      }
      return clipcase;
    }

    /**
     * \brief Compute the weight of a clip value along an edge (id0, id1) using the clip field and value.
     *
     * \param id0 The mesh node at the start of the edge.
     * \param id1 The mesh node at the end of the edge.
     *
     * \return A parametric position t [0,1] where we locate \a clipValues in [d0,d1].
     */
    AXOM_HOST_DEVICE
    ClipFieldType computeWeight(axom::IndexType AXOM_UNUSED_PARAM(zoneIndex),
                                ConnectivityType id0,
                                ConnectivityType id1) const
    {
      const ClipFieldType d0 = m_clipFieldView[id0];
      const ClipFieldType d1 = m_clipFieldView[id1];
      constexpr ClipFieldType tiny = 1.e-09;
      return axom::utilities::clampVal(axom::utilities::abs(m_clipValue - d0) /
                                         (axom::utilities::abs(d1 - d0) + tiny),
                                       ClipFieldType(0),
                                       ClipFieldType(1));
    }

    axom::ArrayView<ClipFieldType> m_clipFieldView {};
    ClipFieldType m_clipValue {};
  };

  /**
   * \brief Initialize the object from options.
   * \param n_options The node that contains the options.
   * \param n_fields The node that contains fields.
   */
  void initialize(const conduit::Node &n_options, const conduit::Node &n_fields)
  {
    namespace bputils = axom::mir::utilities::blueprint;
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    // Get the clip field and value.
    ClipOptions opts(n_options);
    std::string clipFieldName = opts.clipField();
    m_view.m_clipValue = opts.clipValue();

    // Make sure the clipField is the right data type and store access to it in the view.
    const conduit::Node &n_clip_field = n_fields.fetch_existing(opts.clipField());
    const conduit::Node &n_clip_field_values = n_clip_field["values"];
    SLIC_ASSERT(n_clip_field["association"].as_string() == "vertex");
    if(n_clip_field_values.dtype().id() == bputils::cpp2conduit<ClipFieldType>::id)
    {
      // Make a view.
      m_view.m_clipFieldView =
        bputils::make_array_view<ClipFieldType>(n_clip_field_values);
    }
    else
    {
      // Convert to ClipFieldType.
      const IndexType n =
        static_cast<IndexType>(n_clip_field_values.dtype().number_of_elements());
      m_clipFieldData = axom::Array<ClipFieldType>(n, n, allocatorID);
      auto clipFieldView = m_view.m_clipFieldView = m_clipFieldData.view();
      views::Node_to_ArrayView(n_clip_field_values, [&](auto clipFieldViewSrc) {
        axom::for_all<ExecSpace>(
          n,
          AXOM_LAMBDA(auto index) {
            clipFieldView[index] =
              static_cast<ClipFieldType>(clipFieldViewSrc[index]);
          });
      });
    }
  }

  /**
   * \brief Determine the name of the topology on which to operate.
   * \param n_input The input mesh node.
   * \param n_options The clipping options.
   * \return The name of the toplogy on which to operate.
   */
  std::string getTopologyName(const conduit::Node &n_input,
                              const conduit::Node &n_options) const
  {
    // Get the clipField's topo name.
    ClipOptions opts(n_options);
    const conduit::Node &n_fields = n_input.fetch_existing("fields");
    const conduit::Node &n_clipField = n_fields.fetch_existing(opts.clipField());
    return n_clipField["topology"].as_string();
  }

  /**
   * \brief Return a new instance of the view.
   * \return A new instance of the view.
   */
  View view() const { return m_view; }

private:
  axom::Array<ClipFieldType> m_clipFieldData {};
  View m_view {};
};

//------------------------------------------------------------------------------
/**
 * \accelerated
 * \brief This class clips a topology using a field and puts the new topology into a new Conduit node.
 *
 * \tparam ExecSpace    The execution space where the compute-heavy kernels run.
 * \tparam TopologyView The topology view that can operate on the Blueprint topology.
 * \tparam CoordsetView The coordset view that can operate on the Blueprint coordset.
 * \tparam IntersectPolicy The intersector policy that can helps with cases and weights.
 * \tparam NamingPolicy The policy for making names from arrays of ids.
 */
template <typename ExecSpace,
          typename TopologyView,
          typename CoordsetView,
          typename IntersectPolicy = axom::mir::clipping::
            FieldIntersector<ExecSpace, typename TopologyView::ConnectivityType>,
          typename NamingPolicy = axom::mir::utilities::HashNaming>
class ClipField
{
public:
  using BlendData = axom::mir::utilities::blueprint::BlendData;
  using SliceData = axom::mir::utilities::blueprint::SliceData;
  using ClipTableViews = axom::StackArray<axom::mir::clipping::TableView, 6>;
  using Intersector = IntersectPolicy;

  using BitSet = std::uint32_t;
  using KeyType = typename NamingPolicy::KeyType;
  using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
  using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;
  using ConnectivityType = typename TopologyView::ConnectivityType;
  using BlendGroupBuilderType = BlendGroupBuilder<ExecSpace, NamingPolicy>;
  using ClipFieldType = float;

  /**
   * \brief Constructor
   *
   * \param topoView A topology view suitable for the supplied topology.
   * \param coordsetView A coordset view suitable for the supplied coordset.
   *
   */
  ClipField(const TopologyView &topoView,
            const CoordsetView &coordsetView,
            const Intersector &intersector = Intersector())
    : m_topologyView(topoView)
    , m_coordsetView(coordsetView)
    , m_intersector(intersector)
    , m_clipTables()
  { }

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
    // Get the topo/coordset names in the input.
    ClipOptions opts(n_options);
    const std::string topoName =
      m_intersector.getTopologyName(n_input, n_options);
    const conduit::Node &n_topo =
      n_input.fetch_existing("topologies/" + topoName);
    const std::string coordsetName = n_topo["coordset"].as_string();
    const conduit::Node &n_coordset =
      n_input.fetch_existing("coordsets/" + coordsetName);
    const conduit::Node &n_fields = n_input.fetch_existing("fields");

    execute(n_topo,
            n_coordset,
            n_fields,
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
    namespace bputils = axom::mir::utilities::blueprint;
    const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    AXOM_ANNOTATE_SCOPE("ClipField");

    // Make the selected zones and get the size.
    ClipOptions opts(n_options);
    axom::mir::SelectedZones<ExecSpace> selectedZones(
      m_topologyView.numberOfZones(),
      n_options);
    const auto nzones = selectedZones.view().size();

    // Give the intersector a chance to further initialize.
    {
      AXOM_ANNOTATE_SCOPE("Initialize intersector");
      m_intersector.initialize(n_options, n_fields);
    }

    // Load clip table data and make views.
    m_clipTables.load(m_topologyView.dimension());
    ClipTableViews clipTableViews;
    createClipTableViews(clipTableViews, m_topologyView.dimension());

    // Allocate some memory and store views in ZoneData, FragmentData.
    AXOM_ANNOTATE_BEGIN("allocation");
    axom::Array<int> clipCases(nzones,
                               nzones,
                               allocatorID);  // The clip case for a zone.
    axom::Array<BitSet> pointsUsed(
      nzones,
      nzones,
      allocatorID);  // Which points are used over all selected fragments in a zone
    ZoneData zoneData;
    zoneData.m_clipCasesView = clipCases.view();
    zoneData.m_pointsUsedView = pointsUsed.view();

    axom::Array<IndexType> fragments(
      nzones,
      nzones,
      allocatorID);  // The number of fragments (child zones) produced for a zone.
    axom::Array<IndexType> fragmentsSize(
      nzones,
      nzones,
      allocatorID);  // The connectivity size for all selected fragments in a zone.
    axom::Array<IndexType> fragmentOffsets(nzones, nzones, allocatorID);
    axom::Array<IndexType> fragmentSizeOffsets(nzones, nzones, allocatorID);

    FragmentData fragmentData;
    fragmentData.m_fragmentsView = fragments.view();
    fragmentData.m_fragmentsSizeView = fragmentsSize.view();
    fragmentData.m_fragmentOffsetsView = fragmentOffsets.view();
    fragmentData.m_fragmentSizeOffsetsView = fragmentSizeOffsets.view();

    axom::Array<IndexType> blendGroups(
      nzones,
      nzones,
      allocatorID);  // Number of blend groups in a zone.
    axom::Array<IndexType> blendGroupsLen(
      nzones,
      nzones,
      allocatorID);  // Length of the blend groups in a zone.
    axom::Array<IndexType> blendOffset(
      nzones,
      nzones,
      allocatorID);  // Start of zone's blend group indices
    axom::Array<IndexType> blendGroupOffsets(
      nzones,
      nzones,
      allocatorID);  // Start of zone's blend group offsets in definitions.
    AXOM_ANNOTATE_END("allocation");

    // Make an object to help manage building the blend groups.
    BlendGroupBuilderType builder;
    builder.setBlendGroupSizes(blendGroups.view(), blendGroupsLen.view());

    // Compute sizes and offsets
    computeSizes(clipTableViews, builder, zoneData, fragmentData, opts, selectedZones);
    computeFragmentSizes(fragmentData, selectedZones);
    computeFragmentOffsets(fragmentData);

    IndexType blendGroupsSize = 0, blendGroupLenSize = 0;
    builder.computeBlendGroupSizes(blendGroupsSize, blendGroupLenSize);
    builder.setBlendGroupOffsets(blendOffset.view(), blendGroupOffsets.view());
    builder.computeBlendGroupOffsets();

    // Allocate memory for blend groups.
    AXOM_ANNOTATE_BEGIN("allocation2");
    axom::Array<KeyType> blendNames(blendGroupsSize, blendGroupsSize, allocatorID);
    axom::Array<IndexType> blendGroupSizes(blendGroupsSize,
                                           blendGroupsSize,
                                           allocatorID);
    axom::Array<IndexType> blendGroupStart(blendGroupsSize,
                                           blendGroupsSize,
                                           allocatorID);
    axom::Array<IndexType> blendIds(blendGroupLenSize,
                                    blendGroupLenSize,
                                    allocatorID);
    axom::Array<float> blendCoeff(blendGroupLenSize,
                                  blendGroupLenSize,
                                  allocatorID);

    // Make the blend groups.
    builder.setBlendViews(blendNames.view(),
                          blendGroupSizes.view(),
                          blendGroupStart.view(),
                          blendIds.view(),
                          blendCoeff.view());
    AXOM_ANNOTATE_END("allocation2");
    makeBlendGroups(clipTableViews, builder, zoneData, opts, selectedZones);

    // Make the blend groups unique
    axom::Array<KeyType> uNames;
    axom::Array<axom::IndexType> uIndices;
    {
      AXOM_ANNOTATE_SCOPE("unique");
      axom::mir::utilities::unique<ExecSpace, KeyType>(builder.blendNames(),
                                                       uNames,
                                                       uIndices);
      builder.setUniqueNames(uNames.view(), uIndices.view());
    }
    bputils::BlendData blend = builder.makeBlendData();

    // Make the clipped mesh
    makeTopology(clipTableViews,
                 builder,
                 zoneData,
                 fragmentData,
                 opts,
                 selectedZones,
                 n_newTopo,
                 n_newCoordset,
                 n_newFields);
    makeCoordset(blend, n_coordset, n_newCoordset);

    // Get the fields that we want to process.
    std::map<std::string, std::string> fieldsToProcess;
    int numElementFields = 0;
    if(opts.fields(fieldsToProcess))
    {
      // Fields were present in the options. Count the element fields.
      for(auto it = fieldsToProcess.begin(); it != fieldsToProcess.end(); it++)
      {
        const conduit::Node &n_field = n_fields.fetch_existing(it->first);
        if(n_field.fetch_existing("topology").as_string() == n_topo.name())
        {
          numElementFields +=
            (n_field.fetch_existing("association").as_string() == "element") ? 1
                                                                             : 0;
        }
      }
    }
    else
    {
      // Fields were not present in the options. Select all fields that have the same topology as n_topo.
      for(conduit::index_t i = 0; i < n_fields.number_of_children(); i++)
      {
        const conduit::Node &n_field = n_fields[i];
        if(n_field.fetch_existing("topology").as_string() == n_topo.name())
        {
          numElementFields +=
            (n_field.fetch_existing("association").as_string() == "element") ? 1
                                                                             : 0;

          fieldsToProcess[n_field.name()] = n_field.name();
        }
      }
    }

    // Make slice indices if we have element fields.
    bputils::SliceData slice;
    axom::Array<IndexType> sliceIndices;
    if(numElementFields > 0)
    {
      AXOM_ANNOTATE_SCOPE("sliceIndices");
      sliceIndices = axom::Array<IndexType>(fragmentData.m_finalNumZones,
                                            fragmentData.m_finalNumZones,
                                            allocatorID);
      auto sliceIndicesView = sliceIndices.view();

      // Fill in sliceIndicesView.
      const auto selectedZonesView = selectedZones.view();
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(auto index) {
          const auto zoneIndex = selectedZonesView[index];
          const auto start = fragmentData.m_fragmentOffsetsView[index];
          for(int i = 0; i < fragmentData.m_fragmentsView[index]; i++)
            sliceIndicesView[start + i] = zoneIndex;
        });
      slice.m_indicesView = sliceIndicesView;
    }

    makeFields(blend,
               slice,
               opts.topologyName(n_topo.name()),
               fieldsToProcess,
               n_fields,
               n_newFields);
    makeOriginalElements(fragmentData,
                         opts,
                         selectedZones,
                         n_fields,
                         n_newTopo,
                         n_newFields);
  }

private:
  /**
   * \brief Contains data that describes the number and size of zone fragments in the output.
   */
  struct FragmentData
  {
    IndexType m_finalNumZones {0};
    IndexType m_finalConnSize {0};
    axom::ArrayView<IndexType> m_fragmentsView {};
    axom::ArrayView<IndexType> m_fragmentsSizeView {};
    axom::ArrayView<IndexType> m_fragmentOffsetsView {};
    axom::ArrayView<IndexType> m_fragmentSizeOffsetsView {};
  };

  /**
   * \brief Contains some per-zone data that we want to hold onto between methods.
   */
  struct ZoneData
  {
    axom::ArrayView<int> m_clipCasesView {};
    axom::ArrayView<BitSet> m_pointsUsedView {};
  };

  /**
   * \brief Make a bitset that indicates the parts of the selection that are selected.
   */
  int getSelection(const ClipOptions &opts) const
  {
    int selection = 0;
    if(opts.inside()) axom::utilities::setBitOn(selection, 0);
    if(opts.outside()) axom::utilities::setBitOn(selection, 1);
    SLIC_ASSERT(selection > 0);
    return selection;
  }

  /**
   * \brief Create views for the clip tables of various shapes.
   *
   * \param[out] views The views array that will contain the table views.
   * \param dimension The dimension the topology (so we can load a subset of tables)
   */
  void createClipTableViews(ClipTableViews &views, int dimension)
  {
    AXOM_ANNOTATE_SCOPE("createClipTableViews");
    if(dimension == -1 || dimension == 2)
    {
      views[details::getClipTableIndex(views::Tri_ShapeID)] =
        m_clipTables[ST_TRI].view();
      views[details::getClipTableIndex(views::Quad_ShapeID)] =
        m_clipTables[ST_QUA].view();
    }
    if(dimension == -1 || dimension == 3)
    {
      views[details::getClipTableIndex(views::Tet_ShapeID)] =
        m_clipTables[ST_TET].view();
      views[details::getClipTableIndex(views::Pyramid_ShapeID)] =
        m_clipTables[ST_PYR].view();
      views[details::getClipTableIndex(views::Wedge_ShapeID)] =
        m_clipTables[ST_WDG].view();
      views[details::getClipTableIndex(views::Hex_ShapeID)] =
        m_clipTables[ST_HEX].view();
    }
  }

  /**
   * \brief Iterate over zones and their respective fragments to determine sizes
   *        for fragments and blend groups.
   *
   * \param[in] clipTableViews An object that holds views of the clipping table data.
   * \param[in] builder This object holds views to blend group data and helps with building/access.
   * \param[in] zoneData This object holds views to per-zone data.
   * \param[in] fragmentData This object holds views to per-fragment data.
   * \param[inout] opts Clipping options.
   *
   * \note Objects that we need to capture into kernels are passed by value (they only contain views anyway). Data can be modified through the views.
   */
  void computeSizes(ClipTableViews clipTableViews,
                    BlendGroupBuilderType builder,
                    ZoneData zoneData,
                    FragmentData fragmentData,
                    const ClipOptions &opts,
                    const SelectedZones<ExecSpace> &selectedZones) const
  {
    AXOM_ANNOTATE_SCOPE("computeSizes");
    const auto selection = getSelection(opts);

    auto blendGroupsView = builder.state().m_blendGroupsView;
    auto blendGroupsLenView = builder.state().m_blendGroupsLenView;

    const auto deviceIntersector = m_intersector.view();
    m_topologyView.template for_selected_zones<ExecSpace>(
      selectedZones.view(),
      AXOM_LAMBDA(auto szIndex, auto zoneIndex, const auto &zone) {
        // Get the clip case for the current zone.
        const auto clipcase =
          deviceIntersector.determineClipCase(zoneIndex, zone.getIds());
        zoneData.m_clipCasesView[szIndex] = clipcase;

        // Iterate over the shapes in this clip case to determine the number of blend groups.
        const auto clipTableIndex = details::getClipTableIndex(zone.id());
        const auto &ctView = clipTableViews[clipTableIndex];

        int thisBlendGroups =
          0;  // The number of blend groups produced in this case.
        int thisBlendGroupLen = 0;  // The total length of the blend groups.
        int thisFragments =
          0;  // The number of zone fragments produced in this case.
        int thisFragmentsNumIds =
          0;  // The number of points used to make all the fragment zones.
        BitSet ptused = 0;  // A bitset indicating which ST_XX nodes are used.

        auto it = ctView.begin(clipcase);
        const auto end = ctView.end(clipcase);
        for(; it != end; it++)
        {
          // Get the current shape in the clip case.
          const auto fragment = *it;

          if(fragment[0] == ST_PNT)
          {
            if(details::generatedPointIsSelected(fragment[2], selection))
            {
              const int nIds = static_cast<int>(fragment[3]);

              for(int ni = 0; ni < nIds; ni++)
              {
                const auto pid = fragment[4 + ni];

                // Increase the blend size to include this center point.
                if(pid <= P7)
                {
                  // corner point
                  thisBlendGroupLen++;
                }
                else if(pid >= EA && pid <= EL)
                {
                  // edge point
                  thisBlendGroupLen += 2;
                }
              }

              // This center or face point counts as a blend group.
              thisBlendGroups++;

              // Mark the point used.
              axom::utilities::setBitOn(ptused, N0 + fragment[1]);
            }
          }
          else
          {
            if(details::shapeIsSelected(fragment[1], selection))
            {
              thisFragments++;
              const int nIdsThisFragment = fragment.size() - 2;
              thisFragmentsNumIds += nIdsThisFragment;

              // Mark the points this fragment used.
              for(int i = 2; i < fragment.size(); i++)
              {
                axom::utilities::setBitOn(ptused, fragment[i]);
              }
            }
          }
        }

        // Save the flags for the points that were used in this zone
        zoneData.m_pointsUsedView[szIndex] = ptused;

        // Count which points in the original cell are used.
        for(IndexType pid = P0; pid <= P7; pid++)
        {
          const int incr = axom::utilities::bitIsSet(ptused, pid) ? 1 : 0;

          thisBlendGroupLen += incr;  // {p0}
          thisBlendGroups += incr;
        }

        // Count edges that are used.
        for(IndexType pid = EA; pid <= EL; pid++)
        {
          const int incr = axom::utilities::bitIsSet(ptused, pid) ? 1 : 0;

          thisBlendGroupLen += 2 * incr;  // {p0 p1}
          thisBlendGroups += incr;
        }

        // Save the results.
        fragmentData.m_fragmentsView[szIndex] = thisFragments;
        fragmentData.m_fragmentsSizeView[szIndex] = thisFragmentsNumIds;

        // Set blend group sizes for this zone.
        blendGroupsView[szIndex] = thisBlendGroups;
        blendGroupsLenView[szIndex] = thisBlendGroupLen;
      });

#if defined(AXOM_DEBUG_CLIP_FIELD)
    std::cout
      << "------------------------ computeSizes ------------------------"
      << std::endl;
    details::printHost("fragmentData.m_fragmentsView",
                       fragmentData.m_fragmentsView);
    details::printHost("fragmentData.m_fragmentsSizeView",
                       fragmentData.m_fragmentsSizeView);
    details::printHost("blendGroupsView", blendGroupsView);
    details::printHost("blendGroupsLenView", blendGroupsLenView);
    details::printHost("zoneData.m_pointsUsedView", zoneData.m_pointsUsedView);
    details::printHost("zoneData.m_clipCasesView", zoneData.m_clipCasesView);
    std::cout
      << "--------------------------------------------------------------"
      << std::endl;
#endif
  }

  /**
   * \brief Compute the total number of fragments and their size.
   *
   * \param[inout] fragmentData The object that contains data about the zone fragments.
   */
  void computeFragmentSizes(FragmentData &fragmentData,
                            const SelectedZones<ExecSpace> &selectedZones) const
  {
    AXOM_ANNOTATE_SCOPE("computeFragmentSizes");
    const auto nzones = selectedZones.view().size();

    // Sum the number of fragments.
    RAJA::ReduceSum<reduce_policy, IndexType> fragment_sum(0);
    const auto fragmentsView = fragmentData.m_fragmentsView;
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(auto szIndex) { fragment_sum += fragmentsView[szIndex]; });
    fragmentData.m_finalNumZones = fragment_sum.get();

    // Sum the fragment connectivity sizes.
    RAJA::ReduceSum<reduce_policy, IndexType> fragment_nids_sum(0);
    const auto fragmentsSizeView = fragmentData.m_fragmentsSizeView;
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(auto szIndex) {
        fragment_nids_sum += fragmentsSizeView[szIndex];
      });
    fragmentData.m_finalConnSize = fragment_nids_sum.get();
  }

  /**
   * \brief Compute fragment offsets.
   *
   * \param[inout] fragmentData The object that contains data about the zone fragments.
   */
  void computeFragmentOffsets(FragmentData &fragmentData) const
  {
    AXOM_ANNOTATE_SCOPE("computeFragmentOffsets");
    axom::exclusive_scan<ExecSpace>(fragmentData.m_fragmentsView,
                                    fragmentData.m_fragmentOffsetsView);
    axom::exclusive_scan<ExecSpace>(fragmentData.m_fragmentsSizeView,
                                    fragmentData.m_fragmentSizeOffsetsView);

#if defined(AXOM_DEBUG_CLIP_FIELD)
    std::cout << "------------------------ computeFragmentOffsets "
                 "------------------------"
              << std::endl;
    details::printHost("fragmentData.m_fragmentOffsetsView",
                       fragmentData.m_fragmentOffsetsView);
    details::printHost("fragmentData.m_fragmentSizeOffsetsView",
                       fragmentData.m_fragmentSizeOffsetsView);
    std::cout << "-------------------------------------------------------------"
                 "-----------"
              << std::endl;
#endif
  }

  /**
   * \brief Fill in the data for the blend group views.
   *
   * \param[in] clipTableViews An object that holds views of the clipping table data.
   * \param[in] builder This object holds views to blend group data and helps with building/access.
   * \param[in] zoneData This object holds views to per-zone data.
   * \param[inout] opts Clipping options.
   *
   * \note Objects that we need to capture into kernels are passed by value (they only contain views anyway). Data can be modified through the views.
   */
  void makeBlendGroups(ClipTableViews clipTableViews,
                       BlendGroupBuilderType builder,
                       ZoneData zoneData,
                       const ClipOptions &opts,
                       const SelectedZones<ExecSpace> &selectedZones) const
  {
    AXOM_ANNOTATE_SCOPE("makeBlendGroups");
    const auto selection = getSelection(opts);

    const auto deviceIntersector = m_intersector.view();
    m_topologyView.template for_selected_zones<ExecSpace>(
      selectedZones.view(),
      AXOM_LAMBDA(auto szIndex, auto zoneIndex, const auto &zone) {
        // Get the clip case for the current zone.
        const auto clipcase = zoneData.m_clipCasesView[szIndex];

        // Iterate over the shapes in this clip case to determine the number of blend groups.
        const auto clipTableIndex = details::getClipTableIndex(zone.id());
        const auto &ctView = clipTableViews[clipTableIndex];

        // These are the points used in this zone's fragments.
        const BitSet ptused = zoneData.m_pointsUsedView[szIndex];

        // Get the blend groups for this zone.
        auto groups = builder.blendGroupsForZone(szIndex);

        auto it = ctView.begin(clipcase);
        const auto end = ctView.end(clipcase);
        for(; it != end; it++)
        {
          // Get the current shape in the clip case.
          const auto fragment = *it;

          if(fragment[0] == ST_PNT)
          {
            if(details::generatedPointIsSelected(fragment[2], selection))
            {
              const int nIds = static_cast<int>(fragment[3]);
              const auto one_over_n = 1.f / static_cast<float>(nIds);

              groups.beginGroup();
              for(int ni = 0; ni < nIds; ni++)
              {
                const auto ptid = fragment[4 + ni];

                // Add the point to the blend group.
                if(ptid <= P7)
                {
                  // corner point.
                  groups.add(zone.getId(ptid), one_over_n);
                }
                else if(ptid >= EA && ptid <= EL)
                {
                  // edge point.
                  const auto edgeIndex = ptid - EA;
                  const auto edge = zone.getEdge(edgeIndex);
                  const auto id0 = zone.getId(edge[0]);
                  const auto id1 = zone.getId(edge[1]);

                  // Figure out the blend for edge.
                  const auto t =
                    deviceIntersector.computeWeight(zoneIndex, id0, id1);

                  groups.add(id0, one_over_n * (1.f - t));
                  groups.add(id1, one_over_n * t);
                }
              }
              groups.endGroup();
            }
          }
        }

        // Add blend group for each original point that was used.
        for(IndexType pid = P0; pid <= P7; pid++)
        {
          if(axom::utilities::bitIsSet(ptused, pid))
          {
            groups.beginGroup();
            groups.add(zone.getId(pid), 1.f);
            groups.endGroup();
          }
        }

        // Add blend group for each edge point that was used.
        for(IndexType pid = EA; pid <= EL; pid++)
        {
          if(axom::utilities::bitIsSet(ptused, pid))
          {
            const auto edgeIndex = pid - EA;
            const auto edge = zone.getEdge(edgeIndex);
            const auto id0 = zone.getId(edge[0]);
            const auto id1 = zone.getId(edge[1]);

            // Figure out the blend for edge.
            const auto t = deviceIntersector.computeWeight(zoneIndex, id0, id1);

            // Close to the endpoints, just count the edge blend group
            // as an endpoint to ensure better blend group matching later.
            constexpr decltype(t) LOWER = 1.e-4;
            constexpr decltype(t) UPPER = 1. - LOWER;
            groups.beginGroup();
            if(t < LOWER)
            {
              groups.add(id0, 1.f);
            }
            else if(t > UPPER)
            {
              groups.add(id1, 1.f);
            }
            else
            {
              groups.add(id0, 1.f - t);
              groups.add(id1, t);
            }
            groups.endGroup();
          }
        }
      });
  }

  /**
   * \brief Make the clipped mesh topology.
   *
   * \param[in] clipTableViews An object that holds views of the clipping table data.
   * \param[in] builder This object holds views to blend group data and helps with building/access.
   * \param[in] zoneData This object holds views to per-zone data.
   * \param[in] fragmentData This object holds views to per-fragment data.
   * \param[in] opts Clipping options.
   * \param[in] selectedZones The selected zones.
   * \param[out] n_newTopo The node that will contain the new topology.
   * \param[out] n_newCoordset The node that will contain the new coordset.
   * \param[out] n_newFields The node that will contain the new fields.
   *
   * \note Objects that we need to capture into kernels are passed by value (they only contain views anyway). Data can be modified through the views.
   */
  void makeTopology(ClipTableViews clipTableViews,
                    BlendGroupBuilderType builder,
                    ZoneData zoneData,
                    FragmentData fragmentData,
                    const ClipOptions &opts,
                    const SelectedZones<ExecSpace> &selectedZones,
                    conduit::Node &n_newTopo,
                    conduit::Node &n_newCoordset,
                    conduit::Node &n_newFields) const
  {
    AXOM_ANNOTATE_SCOPE("makeTopology");
    namespace bputils = axom::mir::utilities::blueprint;
    constexpr auto connTypeID = bputils::cpp2conduit<ConnectivityType>::id;
    const auto selection = getSelection(opts);

    AXOM_ANNOTATE_BEGIN("allocation");
    n_newTopo.reset();
    n_newTopo["type"] = "unstructured";
    n_newTopo["coordset"] = n_newCoordset.name();

    // Get the ID of a Conduit allocator that will allocate through Axom with device allocator allocatorID.
    utilities::blueprint::ConduitAllocateThroughAxom<ExecSpace> c2a;
    const int conduitAllocatorID = c2a.getConduitAllocatorID();

    // Allocate connectivity.
    conduit::Node &n_conn = n_newTopo["elements/connectivity"];
    n_conn.set_allocator(conduitAllocatorID);
    n_conn.set(conduit::DataType(connTypeID, fragmentData.m_finalConnSize));
    auto connView = bputils::make_array_view<ConnectivityType>(n_conn);

    // Allocate shapes.
    conduit::Node &n_shapes = n_newTopo["elements/shapes"];
    n_shapes.set_allocator(conduitAllocatorID);
    n_shapes.set(conduit::DataType(connTypeID, fragmentData.m_finalNumZones));
    auto shapesView = bputils::make_array_view<ConnectivityType>(n_shapes);

    // Allocate sizes.
    conduit::Node &n_sizes = n_newTopo["elements/sizes"];
    n_sizes.set_allocator(conduitAllocatorID);
    n_sizes.set(conduit::DataType(connTypeID, fragmentData.m_finalNumZones));
    auto sizesView = bputils::make_array_view<ConnectivityType>(n_sizes);

    // Allocate offsets.
    conduit::Node &n_offsets = n_newTopo["elements/offsets"];
    n_offsets.set_allocator(conduitAllocatorID);
    n_offsets.set(conduit::DataType(connTypeID, fragmentData.m_finalNumZones));
    auto offsetsView = bputils::make_array_view<ConnectivityType>(n_offsets);

    // Allocate a color variable to keep track of the "color" of the fragments.
    conduit::Node &n_color = n_newFields[opts.colorField()];
    n_color["topology"] = opts.topologyName(n_newTopo.name());
    n_color["association"] = "element";
    conduit::Node &n_color_values = n_color["values"];
    n_color_values.set_allocator(conduitAllocatorID);
    n_color_values.set(conduit::DataType::int32(fragmentData.m_finalNumZones));
    auto colorView = bputils::make_array_view<int>(n_color_values);

    // Fill in connectivity values in case we leave empty slots later.
    axom::for_all<ExecSpace>(
      connView.size(),
      AXOM_LAMBDA(auto index) { connView[index] = 0; });

#if defined(AXOM_DEBUG_CLIP_FIELD)
    // Initialize the values beforehand. For debugging.
    axom::for_all<ExecSpace>(
      shapesView.size(),
      AXOM_LAMBDA(auto index) {
        shapesView[index] = -2;
        sizesView[index] = -3;
        offsetsView[index] = -4;
        colorView[index] = -5;
      });
#endif
    AXOM_ANNOTATE_END("allocation");

    // Here we fill in the new connectivity, sizes, shapes.
    // We get the node ids from the unique blend names, de-duplicating points when making the new connectivity.
    //
    // NOTE: During development, I ran into problems with this kernel not executing
    //       due to point_2_new being too large. The solution was to reduce the values
    //       for EA-EL, N0-N3 to shrink the array to the point where it can fit in
    //       memory available to the thread.
    //
#if defined(AXOM_CLIP_FILTER_DEGENERATES)
    RAJA::ReduceBitOr<reduce_policy, BitSet> degenerates_reduce(0);
#endif
    {
      AXOM_ANNOTATE_SCOPE("build");
      m_topologyView.template for_selected_zones<ExecSpace>(
        selectedZones.view(),
        AXOM_LAMBDA(auto szIndex,
                    auto AXOM_UNUSED_PARAM(zoneIndex),
                    const auto &zone) {
          // If there are no fragments, return from lambda.
          if(fragmentData.m_fragmentsView[szIndex] == 0) return;

          // Seek to the start of the blend groups for this zone.
          auto groups = builder.blendGroupsForZone(szIndex);

          // Go through the points in the order they would have been added as blend
          // groups, get their blendName, and then overall index of that blendName
          // in uNames, the unique list of new dof names. That will be their index
          // in the final points.
          const BitSet ptused = zoneData.m_pointsUsedView[szIndex];
          ConnectivityType point_2_new[N3 + 1];
          for(BitSet pid = N0; pid <= N3; pid++)
          {
            if(axom::utilities::bitIsSet(ptused, pid))
            {
              point_2_new[pid] = groups.uniqueBlendGroupIndex();
              groups++;
            }
          }
          for(BitSet pid = P0; pid <= P7; pid++)
          {
            if(axom::utilities::bitIsSet(ptused, pid))
            {
              point_2_new[pid] = groups.uniqueBlendGroupIndex();
              groups++;
            }
          }
          for(BitSet pid = EA; pid <= EL; pid++)
          {
            if(axom::utilities::bitIsSet(ptused, pid))
            {
              point_2_new[pid] = groups.uniqueBlendGroupIndex();
              groups++;
            }
          }

          // This is where the output fragment connectivity start for this zone
          int outputIndex = fragmentData.m_fragmentSizeOffsetsView[szIndex];
          // This is where the output fragment sizes/shapes start for this zone.
          int sizeIndex = fragmentData.m_fragmentOffsetsView[szIndex];
#if defined(AXOM_CLIP_FILTER_DEGENERATES)
          bool degenerates = false;
#endif
          // Iterate over the selected fragments and emit connectivity for them.
          const auto clipcase = zoneData.m_clipCasesView[szIndex];
          const auto clipTableIndex = details::getClipTableIndex(zone.id());
          const auto ctView = clipTableViews[clipTableIndex];
          auto it = ctView.begin(clipcase);
          const auto end = ctView.end(clipcase);
          for(; it != end; it++)
          {
            // Get the current shape in the clip case.
            const auto fragment = *it;
            const auto fragmentShape = fragment[0];

            if(fragmentShape != ST_PNT)
            {
              if(details::shapeIsSelected(fragment[1], selection))
              {
                // Output the nodes used in this zone.
                const int fragmentSize = fragment.size();
#if defined(AXOM_CLIP_FILTER_DEGENERATES)
                int connStart = outputIndex;
#endif
                offsetsView[sizeIndex] = outputIndex;
                for(int i = 2; i < fragmentSize; i++)
                  connView[outputIndex++] = point_2_new[fragment[i]];

                const auto nIdsThisFragment = fragmentSize - 2;
#if defined(AXOM_CLIP_FILTER_DEGENERATES)
                // Check for degenerate
                int nUniqueIds = details::unique_count<ConnectivityType, 8>(
                  connView.data() + connStart,
                  nIdsThisFragment);
                bool thisFragmentDegenerate = nUniqueIds < (nIdsThisFragment - 1);
                degenerates |= thisFragmentDegenerate;

                // Rewind the outputIndex so we don't emit it in the connectivity.
                if(thisFragmentDegenerate)
                {
                  //std::cout << "degenerate " << szIndex << " {";
                  //for(int i = 0; i < nIdsThisFragment; i++)
                  //{
                  //   std::cout << connView[connStart + i] << ", ";
                  //}
                  //std::cout << std::endl;

                  outputIndex = connStart;

                  // There is one less fragment than we're expecting in the output.
                  fragmentData.m_fragmentsView[szIndex] -= 1;
                }
                // Mark empty size.
                sizesView[sizeIndex] =
                  thisFragmentDegenerate ? 0 : nIdsThisFragment;
#else
                sizesView[sizeIndex] = nIdsThisFragment;
#endif
                shapesView[sizeIndex] =
                  details::ST_Index_to_ShapeID(fragmentShape);
                colorView[sizeIndex] = fragment[1] - COLOR0;
                sizeIndex++;
              }
            }
          }

#if defined(AXOM_CLIP_FILTER_DEGENERATES)
          // Reduce overall whether there are degenerates.
          degenerates_reduce |= degenerates;
#endif
        });

#if defined(AXOM_DEBUG_CLIP_FIELD)
      std::cout
        << "------------------------ makeTopology ------------------------"
        << std::endl;
      std::cout << "degenerates_reduce=" << degenerates_reduce.get() << std::endl;
      //      details::printHost("selectedZones", selectedZones.view());
      details::printHost("m_fragmentsView", fragmentData.m_fragmentsView);
      //      details::printHost("zoneData.m_clipCasesView", zoneData.m_clipCasesView);
      //      details::printHost("zoneData.m_pointsUsedView", zoneData.m_pointsUsedView);
      details::printHost("conn", connView);
      details::printHost("sizes", sizesView);
      details::printHost("offsets", offsetsView);
      details::printHost("shapes", shapesView);
      details::printHost("color", colorView);
      std::cout
        << "--------------------------------------------------------------"
        << std::endl;
#endif
    }

#if defined(AXOM_CLIP_FILTER_DEGENERATES)
    // We get into this block when degenerate zones were detected where
    // all of their nodes are the same. We need to filter those out.
    if(degenerates_reduce.get())
    {
      AXOM_ANNOTATE_SCOPE("degenerates");

      // There were degenerates so the expected number of fragments per zone (m_fragmentsView)
      // was adjusted down. That means redoing the offsets. These need to be up
      // to date to handle zonal fields later.
      axom::exclusive_scan<ExecSpace>(fragmentData.m_fragmentsView,
                                      fragmentData.m_fragmentOffsetsView);

      // Use sizesView to make a mask that has 1's where size > 0.
      axom::IndexType nz = fragmentData.m_finalNumZones;
      axom::Array<int> mask(nz,
                            nz,
                            axom::execution_space<ExecSpace>::allocatorID());
      axom::Array<int> maskOffsets(
        nz,
        nz,
        axom::execution_space<ExecSpace>::allocatorID());
      auto maskView = mask.view();
      auto maskOffsetsView = maskOffsets.view();
      RAJA::ReduceSum<reduce_policy, axom::IndexType> mask_reduce(0);
      axom::for_all<ExecSpace>(
        nz,
        AXOM_LAMBDA(auto index) {
          const int ival = (sizesView[index] > 0) ? 1 : 0;
          maskView[index] = ival;
          mask_reduce += ival;
        });
      const axom::IndexType filteredZoneCount = mask_reduce.get();

      // Make offsets
      axom::exclusive_scan<ExecSpace>(maskView, maskOffsetsView);

      // Replace data in the input Conduit node with a denser version using the mask.
      auto filter =
        [&](conduit::Node &n_src, const auto srcView, axom::IndexType newSize) {
          using value_type = typename decltype(srcView)::value_type;
          conduit::Node n_values;
          n_values.set_allocator(conduitAllocatorID);
          n_values.set(
            conduit::DataType(bputils::cpp2conduit<value_type>::id, newSize));
          auto valuesView = bputils::make_array_view<value_type>(n_values);
          const auto nValues = maskView.size();
          axom::for_all<ExecSpace>(
            nValues,
            AXOM_LAMBDA(auto index) {
              if(maskView[index] > 0)
              {
                const auto destIndex = maskOffsetsView[index];
                valuesView[destIndex] = srcView[index];
              }
            });

          n_src.swap(n_values);
          return bputils::make_array_view<value_type>(n_src);
        };

      // Filter sizes, shapes, color using the mask
      sizesView = filter(n_sizes, sizesView, filteredZoneCount);
      offsetsView = filter(n_offsets, offsetsView, filteredZoneCount);
      shapesView = filter(n_shapes, shapesView, filteredZoneCount);
      colorView = filter(n_color_values, colorView, filteredZoneCount);

      // Record the filtered size.
      fragmentData.m_finalNumZones = filteredZoneCount;
    }
#endif

    // Figure out which shapes were used.
    BitSet shapesUsed {};
    {
      AXOM_ANNOTATE_SCOPE("shapesUsed");
      RAJA::ReduceBitOr<reduce_policy, BitSet> shapesUsed_reduce(0);
      axom::for_all<ExecSpace>(
        shapesView.size(),
        AXOM_LAMBDA(auto index) {
          BitSet shapeBit {};
          axom::utilities::setBitOn(shapeBit, shapesView[index]);
          shapesUsed_reduce |= shapeBit;
        });
      shapesUsed = shapesUsed_reduce.get();
    }

#if defined(AXOM_DEBUG_CLIP_FIELD)
    std::cout
      << "------------------------ makeTopology ------------------------"
      << std::endl;
    details::printHost("selectedZones", selectedZones.view());
    details::printHost("m_fragmentsView", fragmentData.m_fragmentsView);
    details::printHost("zoneData.m_clipCasesView", zoneData.m_clipCasesView);
    details::printHost("zoneData.m_pointsUsedView", zoneData.m_pointsUsedView);
    details::printHost("conn", connView);
    details::printHost("sizes", sizesView);
    details::printHost("offsets", offsetsView);
    details::printHost("shapes", shapesView);
    details::printHost("color", colorView);
    std::cout
      << "--------------------------------------------------------------"
      << std::endl;
#endif

    // If inside and outside are not selected, remove the color field since we should not need it.
    if(!(opts.inside() && opts.outside()))
    {
      n_newFields.remove(opts.colorField());
    }

#if 1
    // Handle some quad->tri degeneracies
    if(axom::utilities::bitIsSet(shapesUsed, views::Quad_ShapeID))
    {
      AXOM_ANNOTATE_SCOPE("quadtri");
      const axom::IndexType numOutputZones = shapesView.size();
      RAJA::ReduceBitOr<reduce_policy, BitSet> shapesUsed_reduce(0);
      axom::for_all<ExecSpace>(
        numOutputZones,
        AXOM_LAMBDA(auto index) {
          if(shapesView[index] == views::Quad_ShapeID)
          {
            const auto offset = offsetsView[index];
            ConnectivityType pts[4];
            int npts = 0;
            for(int current = 0; current < 4; current++)
            {
              int next = (current + 1) % 4;
              ConnectivityType curNode = connView[offset + current];
              ConnectivityType nextNode = connView[offset + next];
              if(curNode != nextNode) pts[npts++] = curNode;
            }

            if(npts == 3)
            {
              shapesView[index] = views::Tri_ShapeID;
              sizesView[index] = 3;
              connView[offset] = pts[0];
              connView[offset + 1] = pts[1];
              connView[offset + 2] = pts[2];
              connView[offset + 3] =
                pts[2];  // Repeat the last point (it won't be used though).
            }
          }

          BitSet shapeBit {};
          axom::utilities::setBitOn(shapeBit, shapesView[index]);
          shapesUsed_reduce |= shapeBit;
        });
      // We redid shapesUsed reduction in case triangles appeared.
      shapesUsed = shapesUsed_reduce.get();
    }
#endif

    // Add shape information to the connectivity.
    SLIC_ASSERT_MSG(shapesUsed != 0, "No shapes were produced!");
    const auto shapeMap = shapeMap_FromFlags(shapesUsed);
    SLIC_ASSERT_MSG(shapeMap.empty() == false, "The shape map is empty!");
    if(axom::utilities::countBits(shapesUsed) > 1)
    {
      n_newTopo["elements/shape"] = "mixed";
      conduit::Node &n_shape_map = n_newTopo["elements/shape_map"];
      for(auto it = shapeMap.cbegin(); it != shapeMap.cend(); it++)
        n_shape_map[it->first] = it->second;
    }
    else
    {
      n_newTopo["elements"].remove("shapes");
      n_newTopo["elements/shape"] = shapeMap.begin()->first;
    }
  }

  /**
   * \brief Make the new coordset using the blend data and the input coordset/coordsetview.
   *
   * \param blend The BlendData that we need to construct the new coordset.
   * \param n_coordset The input coordset, which is passed for metadata.
   * \param[out] n_newCoordset The new coordset.
   */
  void makeCoordset(const BlendData &blend,
                    const conduit::Node &n_coordset,
                    conduit::Node &n_newCoordset) const
  {
    AXOM_ANNOTATE_SCOPE("makeCoordset");
    axom::mir::utilities::blueprint::
      CoordsetBlender<ExecSpace, CoordsetView, axom::mir::utilities::blueprint::SelectSubsetPolicy>
        cb;
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
    AXOM_ANNOTATE_SCOPE("makeFields");
    for(auto it = fieldMap.begin(); it != fieldMap.end(); it++)
    {
      const conduit::Node &n_field = n_fields.fetch_existing(it->first);
      const std::string association = n_field["association"].as_string();
      if(association == "element")
      {
        bool handled = false;
        // Conditionally support strided-structured.
        if constexpr(axom::mir::views::view_traits<
                       TopologyView>::supports_strided_structured())
        {
          if(n_field.has_path("offsets") && n_field.has_path("strides"))
          {
            using Indexing = typename TopologyView::IndexingPolicy;
            using IndexingPolicy =
              axom::mir::utilities::blueprint::SSElementFieldIndexing<Indexing>;
            IndexingPolicy indexing;
            indexing.m_indexing = m_topologyView.indexing();
            indexing.update(n_field);

            axom::mir::utilities::blueprint::FieldSlicer<ExecSpace, IndexingPolicy> s(
              indexing);
            s.execute(slice, n_field, n_out_fields[it->second]);
            handled = true;
          }
        }
        if(!handled)
        {
          axom::mir::utilities::blueprint::FieldSlicer<ExecSpace> s;
          s.execute(slice, n_field, n_out_fields[it->second]);
        }

        n_out_fields[it->second]["topology"] = topologyName;
      }
      else if(association == "vertex")
      {
        bool handled = false;

        // Node indices in the blend groups are global indices. This means that, provided the
        // field's offsets/strides match the topology's (and why would they not?), we can skip
        // strided structured support for now. Enable this code if the field offsets/strides
        // do not match the topo's offsets/strides.

        // Conditionally support strided-structured.
        if constexpr(axom::mir::views::view_traits<
                       TopologyView>::supports_strided_structured())
        {
          if(n_field.has_path("offsets") && n_field.has_path("strides"))
          {
            // Make node indexing that the field blender can use.
            using Indexing = typename TopologyView::IndexingPolicy;
            using IndexingPolicy =
              axom::mir::utilities::blueprint::SSVertexFieldIndexing<Indexing>;
            IndexingPolicy indexing;
            indexing.m_topoIndexing = m_topologyView.indexing().expand();
            indexing.m_fieldIndexing = m_topologyView.indexing().expand();
            indexing.update(n_field);

            // If the topo and field offsets/strides are different then we need to go through
            // SSVertexFieldIndexing. Otherwise, we can let the normal case further below
            // handle the field.
            if(indexing.m_topoIndexing.m_offsets !=
                 indexing.m_fieldIndexing.m_offsets ||
               indexing.m_topoIndexing.m_strides !=
                 indexing.m_fieldIndexing.m_strides)
            {
              // Blend the field.
              axom::mir::utilities::blueprint::FieldBlender<
                ExecSpace,
                axom::mir::utilities::blueprint::SelectSubsetPolicy,
                IndexingPolicy>
                b(indexing);
              b.execute(blend, n_field, n_out_fields[it->second]);
              handled = true;
            }
          }
        }
        if(!handled)
        {
          // Blend the field.
          axom::mir::utilities::blueprint::FieldBlender<
            ExecSpace,
            axom::mir::utilities::blueprint::SelectSubsetPolicy>
            b;
          b.execute(blend, n_field, n_out_fields[it->second]);
        }

        n_out_fields[it->second]["topology"] = topologyName;
      }
    }
  }

  /**
   * \brief Make an originalElements field so we can know each output zone's original zone number in the input mesh.
   *
   * \param[in] fragmentData This object holds views to per-fragment data.
   * \param[in] opts Clipping options.
   * \param[in] selectedZones The selected zones.
   * \param[in] n_fields The node that contains the input mesh's fields.
   * \param[out] n_newTopo The node that will contain the new topology.
   * \param[out] n_newFields The node that will contain the new fields.
   *
   * \note Objects that we need to capture into kernels are passed by value (they only contain views anyway). Data can be modified through the views.
   */
  void makeOriginalElements(FragmentData fragmentData,
                            const ClipOptions &opts,
                            const SelectedZones<ExecSpace> &selectedZones,
                            const conduit::Node &n_fields,
                            conduit::Node &n_newTopo,
                            conduit::Node &n_newFields) const
  {
    AXOM_ANNOTATE_SCOPE("makeOriginalElements");
    namespace bputils = axom::mir::utilities::blueprint;
    constexpr auto connTypeID = bputils::cpp2conduit<ConnectivityType>::id;

    utilities::blueprint::ConduitAllocateThroughAxom<ExecSpace> c2a;
    const int conduitAllocatorID = c2a.getConduitAllocatorID();

    const auto selectedZonesView = selectedZones.view();
    const auto nzones = selectedZonesView.size();
    const std::string originalElements(opts.originalElementsField());

    if(n_fields.has_child(originalElements))
    {
      // originalElements already exists. We need to map it forward.
      const conduit::Node &n_orig = n_fields[originalElements];
      const conduit::Node &n_orig_values = n_orig["values"];
      views::IndexNode_to_ArrayView(n_orig_values, [&](auto origValuesView) {
        using value_type = typename decltype(origValuesView)::value_type;
        conduit::Node &n_origElem = n_newFields[originalElements];
        n_origElem["association"] = "element";
        n_origElem["topology"] = opts.topologyName(n_newTopo.name());
        conduit::Node &n_values = n_origElem["values"];
        n_values.set_allocator(conduitAllocatorID);
        n_values.set(conduit::DataType(n_orig_values.dtype().id(),
                                       fragmentData.m_finalNumZones));
        auto valuesView = bputils::make_array_view<value_type>(n_values);
        axom::for_all<ExecSpace>(
          nzones,
          AXOM_LAMBDA(auto index) {
            const int sizeIndex = fragmentData.m_fragmentOffsetsView[index];
            const int nFragments = fragmentData.m_fragmentsView[index];
            const auto zoneIndex = selectedZonesView[index];
            for(int i = 0; i < nFragments; i++)
              valuesView[sizeIndex + i] = origValuesView[zoneIndex];
          });
      });
    }
    else
    {
      // Make a new node and populate originalElement.
      conduit::Node &n_orig = n_newFields[originalElements];
      n_orig["association"] = "element";
      n_orig["topology"] = opts.topologyName(n_newTopo.name());
      conduit::Node &n_values = n_orig["values"];
      n_values.set_allocator(conduitAllocatorID);
      n_values.set(conduit::DataType(connTypeID, fragmentData.m_finalNumZones));
      auto valuesView = bputils::make_array_view<ConnectivityType>(n_values);
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(auto index) {
          const int sizeIndex = fragmentData.m_fragmentOffsetsView[index];
          const int nFragments = fragmentData.m_fragmentsView[index];
          const auto zoneIndex = selectedZonesView[index];
          for(int i = 0; i < nFragments; i++)
            valuesView[sizeIndex + i] = zoneIndex;
        });
    }
  }

  /**
   * \brief Given a flag that includes bitwise-or'd shape ids, make a map that indicates which Conduit shapes are used.
   *
   * \param shapes This is a bitwise-or of various (1 << ShapeID) values.
   *
   * \return A map of Conduit shape name to ShapeID value.
   */
  std::map<std::string, int> shapeMap_FromFlags(std::uint64_t shapes) const
  {
    std::map<std::string, int> sm;

    if(axom::utilities::bitIsSet(shapes, views::Line_ShapeID))
      sm["line"] = views::Line_ShapeID;

    if(axom::utilities::bitIsSet(shapes, views::Tri_ShapeID))
      sm["tri"] = views::Tri_ShapeID;

    if(axom::utilities::bitIsSet(shapes, views::Quad_ShapeID))
      sm["quad"] = views::Quad_ShapeID;

    if(axom::utilities::bitIsSet(shapes, views::Polygon_ShapeID))
      sm["polygon"] = views::Polygon_ShapeID;

    if(axom::utilities::bitIsSet(shapes, views::Tet_ShapeID))
      sm["tet"] = views::Tet_ShapeID;

    if(axom::utilities::bitIsSet(shapes, views::Pyramid_ShapeID))
      sm["pyramid"] = views::Pyramid_ShapeID;

    if(axom::utilities::bitIsSet(shapes, views::Wedge_ShapeID))
      sm["wedge"] = views::Wedge_ShapeID;

    if(axom::utilities::bitIsSet(shapes, views::Hex_ShapeID))
      sm["hex"] = views::Hex_ShapeID;

    if(axom::utilities::bitIsSet(shapes, views::Polyhedron_ShapeID))
      sm["polyhedron"] = views::Polyhedron_ShapeID;

    return sm;
  }

private:
  TopologyView m_topologyView;
  CoordsetView m_coordsetView;
  Intersector m_intersector;
  axom::mir::clipping::ClipTableManager<ExecSpace> m_clipTables;
};

}  // end namespace clipping
}  // end namespace mir
}  // end namespace axom

#endif
