// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_EQUIZ_ALGORITHM_HPP_
#define AXOM_MIR_EQUIZ_ALGORITHM_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/mir.hpp"

// Include these directly for now.
#include "axom/mir/views/MaterialView.hpp"
#include "axom/mir/MIRAlgorithm.hpp"
#include "axom/mir/RecenterField.hpp"
#include "axom/mir/NodeToZoneRelationBuilder.hpp"

#include <conduit/conduit.hpp>
#include <conduit/conduit_relay_io_blueprint.hpp>

#include <string>

#if defined(AXOM_USE_RAJA)
  #include <RAJA/RAJA.hpp>
#endif

namespace axom
{
namespace mir
{
using MaterialID = int;
using MaterialIDArray = axom::Array<MaterialID>;
using MaterialIDView = axom::ArrayView<MaterialID>;
using MaterialVF = float;
using MaterialVFArray = axom::Array<MaterialVF>;
using MaterialVFView = axom::ArrayView<MaterialVF>;

constexpr static int NULL_MATERIAL = -1;
constexpr static MaterialVF NULL_MATERIAL_VF = 0.f;

#if 1
AXOM_HOST_DEVICE
inline bool zoneMatters(int zoneIndex)
{
   const int zones[] = {26,
33,34,35,
36,37,38,39,40,41,
49,50,
59,60,
69,
78,79};
  for(size_t i = 0; i < sizeof(zones)/sizeof(int); i++)
  {
    if(zones[i] == zoneIndex)
      return true;
  }
  return false;
}

template <typename IntType>
struct RelationView
{
  axom::ArrayView<IntType> m_zonesView {};
  axom::ArrayView<IntType> m_sizesView {};
  axom::ArrayView<IntType> m_offsetsView {};
};

#if 1
/**
 * \brief This class is an intersection policy compatible with ClipField. It
 *        helps determine clip cases and weights using material-aware logic.
 */
template <typename ConnectivityT, int MAXMATERIALS = 10>
class MaterialIntersector
{
public:
  using ConnectivityType = ConnectivityT;
  using ConnectivityView = axom::ArrayView<ConnectivityType>;

  /**
   * \brief This is a view class for MatsetIntersector that can be used in device code.
   */
  struct View
  {
    static constexpr int INVALID_INDEX = -1;

    /**
     * \brief Determine the clipping case, taking into account the zone's material
     *        and the current material being added.
     *
     * \param zoneIndex The zone index in the zoneMaterialView field.
     * \param nodeIds A view containing node ids for the current zone.
     *
     * \return The clip case number for the zone.
     */
    AXOM_HOST_DEVICE
    axom::IndexType determineClipCase(axom::IndexType zoneIndex, const ConnectivityView &nodeIdsView) const
    {
      // Determine the matvf view index for the material that owns the zone.
      int backgroundIndex = INVALID_INDEX;
      int zoneMatID = m_zoneMatNumberView[zoneIndex];
      if(zoneMatID != NULL_MATERIAL)
        backgroundIndex = matNumberToIndex(zoneMatID);
      // Determine the matvf view index for the current material.
      int currentIndex = matNumberToIndex(m_currentMaterial);

      axom::IndexType clipcase = 0;
      const auto n = nodeIdsView.size();
      for(IndexType i = 0; i < n; i++)
      {
        const auto nid = nodeIdsView[i];
        MaterialVF vf1 = (backgroundIndex != INVALID_INDEX) ? m_matvfViews[backgroundIndex][nid] : 0.;
        MaterialVF vf2 = (currentIndex != INVALID_INDEX) ? m_matvfViews[currentIndex][nid] : 0.;

        clipcase |= (vf2 > vf1) ? (1 << i) : 0;
      }
      return clipcase;
    }

    /**
     * \brief Compute the weight of a clip value along an edge (id0, id1) using the clip field and value.
     *
     * \param id0 The mesh node at the start of the edge.
     * \param id1 The mesh node at the end of the edge.
     */
    AXOM_HOST_DEVICE
    float computeWeight(axom::IndexType zoneIndex, ConnectivityType id0, ConnectivityType id1) const
    {
      // Determine the matvf view index for the material that owns the zone.
      int backgroundIndex = INVALID_INDEX;
      int zoneMatID = m_zoneMatNumberView[zoneIndex];
      if(zoneMatID != NULL_MATERIAL)
        backgroundIndex = matNumberToIndex(zoneMatID);
      // Determine the matvf view index for the current material.
      int currentIndex = matNumberToIndex(m_currentMaterial);

      // Get the volume fractions for mat1, mat2 at the edge endpoints id0, id1.
      MaterialVF vf1[2], vf2[2];
      vf1[0] = (backgroundIndex != INVALID_INDEX) ? m_matvfViews[backgroundIndex][id0] : 0.;
      vf1[1] = (backgroundIndex != INVALID_INDEX) ? m_matvfViews[backgroundIndex][id1] : 0.;
      vf2[0] = (currentIndex != INVALID_INDEX) ? m_matvfViews[currentIndex][id0] : 0.;
      vf2[1] = (currentIndex != INVALID_INDEX) ? m_matvfViews[currentIndex][id1] : 0.;

      float numerator = vf2[0] - vf1[0];
      float denominator = -vf1[0] + vf1[1] + vf2[0] - vf2[1];

      float t = 0.f;
      if(denominator != 0.f)
      {
        t = numerator / denominator;
      }
      t = axom::utilities::clampVal(t, 0.f, 1.f);

      return t;
    }

    AXOM_HOST_DEVICE
    int matNumberToIndex(int matNumber) const
    {
      auto index = axom::mir::utilities::bsearch(matNumber, m_matNumbersView);
      return (index != -1) ? m_matIndicesView[index] : INVALID_INDEX;
    }

    /// Helper initialization methods for the host.

    void addMaterial(const MaterialVFView &matvf)
    {
      m_matvfViews.push_back(matvf);
    }

    void setMaterialNumbers(const axom::ArrayView<int> &matNumbersView)
    {
      m_matNumbersView = matNumbersView;
    }

    void setMaterialIndices(const axom::ArrayView<int> &matIndicesView)
    {
      m_matIndicesView = matIndicesView;
    }

    void setZoneMaterialID(const axom::ArrayView<int> &zoneMatsView)
    {
      m_zoneMatNumberView = zoneMatsView;
    }

    void setCurrentMaterial(int matNumber)
    {
      m_currentMaterial = matNumber;
    }

    axom::StaticArray<MaterialVFView, MAXMATERIALS> m_matvfViews {};
    axom::ArrayView<int> m_matNumbersView {};
    axom::ArrayView<int> m_matIndicesView {};
    axom::ArrayView<int> m_zoneMatNumberView {}; //!< Contains the current material number that owns each zone.
    int m_currentMaterial {};                    //!< The current material.
  };

  /**
   * \brief Initialize the object from options.
   * \param n_options The node that contains the options.
   * \param n_fields The node that contains fields.
   */
  void initialize(const conduit::Node &AXOM_UNUSED_PARAM(n_options), const conduit::Node &AXOM_UNUSED_PARAM(n_fields))
  {
  }

  /**
   * \brief Determine the name of the topology on which to operate.
   * \param n_input The input mesh node.
   * \param n_options The clipping options.
   * \return The name of the toplogy on which to operate.
   */
  std::string getTopologyName(const conduit::Node &AXOM_UNUSED_PARAM(n_input), const conduit::Node &n_options) const
  {
    return n_options["topology"].as_string();
  }

  void addMaterial(const MaterialVFView &matvf)
  {
    m_view.addMaterial(matvf);
  }

  void setMaterialNumbers(const axom::ArrayView<int> &matNumbers)
  {
    m_view.setMaterialNumbers(matNumbers);
  }

  void setMaterialIndices(const axom::ArrayView<int> &matIndices)
  {
    m_view.setMaterialIndices(matIndices);
  }

  void setZoneMaterialID(const axom::ArrayView<int> &zoneMatsView)
  {
    m_view.setZoneMaterialID(zoneMatsView);
  }

  void setCurrentMaterial(int matNumber)
  {
    m_view.setCurrentMaterial(matNumber);
  }

  /**
   * \brief Return a new instance of the view.
   * \return A new instance of the view.
   * \note Call this after all values are set.
   */
  View view() const { return m_view; }

private:
  View m_view {};
};
#else
/**
 * \brief This class is an intersection policy compatible with ClipField. It
 *        helps determine clip cases and weights using material-aware logic.
 */
template <typename ConnectivityT>
class MaterialIntersector
{
public:
  using ConnectivityType = ConnectivityT;
  using ConnectivityView = axom::ArrayView<ConnectivityType>;

  MaterialIntersector(const MaterialVFView &accumVFView, const MaterialVFView &currentVFView) : 
    m_view(accumVFView, currentVFView)
  { }

  /**
   * \brief This is a view class for MatsetIntersector that can be used in device code.
   */
  struct View
  {
    /// Constructor
    View() = default;

    /// Constructor
    View(const MaterialVFView &accumVFView, const MaterialVFView &currentVFView) : m_accumVFView(accumVFView), m_currentVFView(currentVFView)
    { }

    /**
     * \brief Determine the clipping case, taking into account the zone's material
     *        and the current material being added.
     *
     * \param zoneIndex The zone index in the zoneMaterialView field.
     * \param nodeIds A view containing node ids for the current zone.
     *
     * \return The clip case number for the zone.
     */
    AXOM_HOST_DEVICE
    axom::IndexType determineClipCase(axom::IndexType AXOM_UNUSED_PARAM(zoneIndex), const ConnectivityView &nodeIdsView) const
    {
      axom::IndexType clipcase = 0;
      const auto n = nodeIdsView.size();
      for(IndexType i = 0; i < n; i++)
      {
        const auto nid = nodeIdsView[i];
        MaterialVF vf1 = m_accumVFView[nid];
        MaterialVF vf2 = m_currentVFView[nid];

        clipcase |= (vf2 > vf1) ? (1 << i) : 0;
      }
      return clipcase;
    }

    /**
     * \brief Compute the weight of a clip value along an edge (id0, id1) using the clip field and value.
     *
     * \param id0 The mesh node at the start of the edge.
     * \param id1 The mesh node at the end of the edge.
     */
    AXOM_HOST_DEVICE
    float computeWeight(axom::IndexType AXOM_UNUSED_PARAM(zoneIndex), ConnectivityType id0, ConnectivityType id1) const
    {
      // Get the volume fractions for mat1, mat2 at the edge endpoints id0, id1.
      MaterialVF vf1[2], vf2[2];
      vf1[0] = m_accumVFView[id0];
      vf1[1] = m_accumVFView[id1];
      vf2[0] = m_currentVFView[id0];
      vf2[1] = m_currentVFView[id1];

      float numerator = vf2[0] - vf1[0];
      float denominator = -vf1[0] + vf1[1] + vf2[0] - vf2[1];

      float t = 0.f;
      if(denominator != 0.f)
      {
        t = numerator / denominator;
      }
      t = axom::utilities::clampVal(t, 0.f, 1.f);

      return t;
    }

    MaterialVFView m_accumVFView {};
    MaterialVFView m_currentVFView {};
  };

  /**
   * \brief Initialize the object from options.
   * \param n_options The node that contains the options.
   * \param n_fields The node that contains fields.
   */
  void initialize(const conduit::Node &AXOM_UNUSED_PARAM(n_options), const conduit::Node &AXOM_UNUSED_PARAM(n_fields))
  {
  }

  /**
   * \brief Determine the name of the topology on which to operate.
   * \param n_input The input mesh node.
   * \param n_options The clipping options.
   * \return The name of the toplogy on which to operate.
   */
  std::string getTopologyName(const conduit::Node &AXOM_UNUSED_PARAM(n_input), const conduit::Node &n_options) const
  {
    return n_options["topology"].as_string();
  }

  /**
   * \brief Return a new instance of the view.
   * \return A new instance of the view.
   */
  View view() const { return m_view; }

private:
  View m_view {};
};
#endif
#endif

/**
 * \accelerated
 * \brief Implements Meredith's Equi-Z algorithm on the GPU using Blueprint inputs/outputs.
 */
template <typename ExecSpace, typename TopologyView, typename CoordsetView, typename MatsetView>
class EquiZAlgorithm : public axom::mir::MIRAlgorithm
{
public:
  /**
   * \brief Constructor
   *
   * \param topoView The topology view to use for the input data.
   * \param coordsetView The coordset view to use for the input data.
   * \param matsetView The matset view to use for the input data.
   */
  EquiZAlgorithm(const TopologyView &topoView,
                 const CoordsetView &coordsetView,
                 const MatsetView &matsetView)
    : m_topologyView(topoView)
    , m_coordsetView(coordsetView)
    , m_matsetView(matsetView)
  { }

  /// Destructor
  virtual ~EquiZAlgorithm() = default;

protected:
  void printNode(const conduit::Node &n) const
  {
    conduit::Node options;
    options["num_children_threshold"] = 10000;
    options["num_elements_threshold"] = 10000;
    n.to_summary_string_stream(std::cout, options);
  }

  /**
   * \brief Perform material interface reconstruction on a single domain.
   *
   * \param[in] n_topo The Conduit node containing the topology that will be used for MIR.
   * \param[in] n_coordset The Conduit node containing the coordset.
   * \param[in] n_fields The Conduit node containing the fields.
   * \param[in] n_matset The Conduit node containing the matset.
   * \param[in] n_options The Conduit node containing the options that help govern MIR execution.
   *
   * \param[out] n_newTopo A node that will contain the new clipped topology.
   * \param[out] n_newCoordset A node that will contain the new coordset for the clipped topology.
   * \param[out] n_newFields A node that will contain the new fields for the clipped topology.
   * \param[out] n_newMatset A Conduit node that will contain the new matset.
   * 
   */
  virtual void executeDomain(const conduit::Node &n_topo,
                             const conduit::Node &n_coordset,
                             const conduit::Node &n_fields,
                             const conduit::Node &n_matset,
                             const conduit::Node &n_options,
                             conduit::Node &n_newTopo,
                             conduit::Node &n_newCoordset,
                             conduit::Node &n_newFields,
                             conduit::Node &n_newMatset) override
  {
    namespace bputils = axom::mir::utilities::blueprint;

    // Copy the options.
    conduit::Node n_options_copy;
    bputils::copy<ExecSpace>(n_options_copy, n_options);
    n_options_copy["topology"] = n_topo.name();

    // Make some nodes that will contain the inputs to subsequent iterations.
    // Store them under a single node so the nodes will have names.
    conduit::Node n_Input;
    conduit::Node &n_InputTopo = n_Input[n_topo.path()];
    conduit::Node &n_InputCoordset = n_Input[n_coordset.path()];
    conduit::Node &n_InputFields = n_Input[n_fields.path()];

    // Get the materials from the matset and determine which of them are clean/mixed.
    axom::mir::views::MaterialInformation allMats, cleanMats, mixedMats;
    classifyMaterials(n_matset, allMats, cleanMats, mixedMats);

    //--------------------------------------------------------------------------
    //
    // Make node-centered VF fields and add various working fields.
    //
    //--------------------------------------------------------------------------
    n_InputFields.reset();
    for(conduit::index_t i = 0; i < n_fields.number_of_children(); i++)
    {
      n_InputFields[n_fields[i].name()].set_external(n_fields[i]);
    }
    makeNodeCenteredVFs(n_topo, n_coordset, n_InputFields, mixedMats);
    makeWorkingFields(n_topo, n_InputFields, cleanMats, mixedMats);

    // Iterate over mixed materials. Note the first material was already accounted for so start at 1.
    for(size_t i = 1; i < mixedMats.size(); i++)
    {
      if(i == 1)
      {
#if 1
        // Print input mesh.
        std::cout << "-------------------------------------------------------------------" << std::endl;
        conduit::Node mesh;
        mesh[n_topo.path()].set_external(n_topo);
        mesh[n_coordset.path()].set_external(n_coordset);
        mesh[n_InputFields.path()].set_external(n_InputFields);
        mesh[n_matset.path()].set_external(n_matset);
        printNode(mesh);
        std::cout << "-------------------------------------------------------------------" << std::endl;
#endif

        // The first time through, we can use the supplied views.
        iteration<TopologyView, CoordsetView>(m_topologyView,
                                              m_coordsetView,

                                              allMats,
                                              mixedMats[i],

                                              n_topo,
                                              n_coordset,
                                              n_InputFields,

                                              n_options_copy,

                                              n_newTopo,
                                              n_newCoordset,
                                              n_newFields);

        // In later iterations, we do not want to pass selectedZones through
        // since they are only valid on the current input topology. Also, if they
        // were passed then the new topology only has those selected zones.
        if(n_options_copy.has_child("selectedZones"))
        {
          n_options_copy.remove("selectedZones");
        }
      }
      else
      {
        // Clear the inputs from the last iteration.
        n_InputTopo.reset();
        n_InputCoordset.reset();
        n_InputFields.reset();

        // Move the outputs of the last iteration to the inputs of this iteration.       
        n_InputTopo.move(n_newTopo);
        n_InputCoordset.move(n_newCoordset);
        n_InputFields.move(n_newFields);

        // The data are now an unstructured view, probably a mixed shape view.
        // Dispatch to an appropriate topo view
        // clangformat-off
        views::dispatch_explicit_coordset(n_InputCoordset, [&](auto coordsetView) {
          using ICoordsetView = decltype(coordsetView);
          using ConnectivityType = typename TopologyView::ConnectivityType;
          views::typed_dispatch_unstructured_topology<ConnectivityType, views::view_traits<TopologyView>::selected_shapes()>(
            n_InputTopo,
            [&](const auto &AXOM_UNUSED_PARAM(shape), auto topologyView) {
              using ITopologyView = decltype(topologyView);

              // Do the next iteration.
              iteration<ITopologyView, ICoordsetView>(
                topologyView,
                coordsetView,

                allMats,
                mixedMats[i],

                n_InputTopo,
                n_InputCoordset,
                n_InputFields,

                n_options_copy,

                n_newTopo,
                n_newCoordset,
                n_newFields);
            });
        });
        // clangformat-on
      }
    }

    // TODO: we need to build the new matset.

  }

  void classifyMaterials(const conduit::Node &n_matset,
                         axom::mir::views::MaterialInformation &allMats,
                         axom::mir::views::MaterialInformation &cleanMats,
                         axom::mir::views::MaterialInformation &mixedMats) const
  {
    cleanMats.clear();
    mixedMats.clear();                      
    allMats = axom::mir::views::materials(n_matset);

    // TODO: actually determine which materials are clean/mixed. It's probably
    //       best to ask the matsetView since it takes some work to determine
    //       this.

    mixedMats = allMats;
  }

  std::string zonalFieldName(int matId) const
  {
    std::stringstream ss;
    ss << "__equiz_zonal_volume_fraction_" << matId;
    return ss.str();
  }

  std::string nodalFieldName(int matId) const
  {
    std::stringstream ss;
    ss << "__equiz_nodal_volume_fraction_" << matId;
    return ss.str();
  }

  std::string accumulatedVFFieldName() const  { return "__equiz_accumulated_VF"; }

  axom::mir::views::MaterialInformation mixedMaterials(const axom::mir::views::MaterialInformation &matInfo) const
  {
    // TODO: figure out which mats are mixed and which are clean.
    return matInfo;
  }

  void makeNodeCenteredVFs(const conduit::Node &n_topo,
                           const conduit::Node &n_coordset,
                           conduit::Node &n_fields,
                           const axom::mir::views::MaterialInformation &mixedMats) const
  {
    namespace bputils = axom::mir::utilities::blueprint;
    // Make a node to zone relation.
    conduit::Node relation;
    {
      bputils::NodeToZoneRelationBuilder<ExecSpace> rb;
      rb.execute(n_topo, n_coordset, relation);
      printNode(relation);
      std::cout.flush();
    }

    // Get the ID of a Conduit allocator that will allocate through Axom with device allocator allocatorID.
    bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;

    // Make nodal VFs for each mixed material.
    const auto nzones = m_topologyView.numberOfZones();
    const auto nnodes = m_coordsetView.numberOfNodes();
    for(const auto &mat : mixedMats)
    {
      const int matNumber = mat.number;
      const std::string zonalName = zonalFieldName(matNumber);
      conduit::Node &n_zonalField = n_fields[zonalName];
      n_zonalField["topology"] = n_topo.name();
      n_zonalField["association"] = "element";
      n_zonalField["values"].set_allocator(c2a.getConduitAllocatorID());
      n_zonalField["values"].set(conduit::DataType(bputils::cpp2conduit<MaterialVF>::id, nzones));
      auto zonalFieldView = bputils::make_array_view<MaterialVF>(n_zonalField["values"]);

      // Fill the zonal field from the matset.
      MatsetView deviceMatsetView(m_matsetView);
      axom::for_all<ExecSpace>(
        m_topologyView.numberOfZones(),
        AXOM_LAMBDA(auto zoneIndex) {
          typename MatsetView::FloatType vf {};
          deviceMatsetView.zoneContainsMaterial(zoneIndex, matNumber, vf);
          zonalFieldView[zoneIndex] = static_cast<MaterialVF>(vf);
      });

      // Make a nodal field for the current material by recentering.
      const std::string nodalName = nodalFieldName(matNumber);
      conduit::Node &n_nodalField = n_fields[nodalName];
      n_nodalField["topology"] = n_topo.name();
      n_nodalField["association"] = "vertex";
      n_nodalField["values"].set_allocator(c2a.getConduitAllocatorID());
      n_nodalField["values"].set(conduit::DataType(bputils::cpp2conduit<MaterialVF>::id, nnodes));
      bputils::RecenterField<ExecSpace> z2n;
      z2n.execute(n_zonalField, relation, n_nodalField);

#ifndef DEBUG_EQUIZ
      n_fields.remove(zonalName);
#endif
    }
  }

  void makeWorkingFields(const conduit::Node &n_topo,
                         conduit::Node &n_fields,
                         const axom::mir::views::MaterialInformation &cleanMats,
                         const axom::mir::views::MaterialInformation &mixedMats) const
  {
    namespace bputils = axom::mir::utilities::blueprint;

    // Get the ID of a Conduit allocator that will allocate through Axom with device allocator allocatorID.
    bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;

    const auto nzones = m_topologyView.numberOfZones();
    const auto nnodes = m_coordsetView.numberOfNodes();

    // Make the accumulation field.
    conduit::Node &n_accumField = n_fields[accumulatedVFFieldName()];
    n_accumField["topology"] = n_topo.name();
    n_accumField["association"] = "vertex";
    n_accumField["values"].set_allocator(c2a.getConduitAllocatorID());
    n_accumField["values"].set(conduit::DataType(bputils::cpp2conduit<MaterialVF>::id, nnodes));
    auto accumFieldView = bputils::make_array_view<MaterialVF>(n_accumField["values"]);
    if(!mixedMats.empty())
    {
      // Copy the first material's data into the accumulated vfs.
      conduit::Node &n_firstMat = n_fields.fetch_existing(nodalFieldName(mixedMats[0].number) + "/values");
      axom::copy(accumFieldView.data(), n_firstMat.data_ptr(), sizeof(MaterialVF) * nnodes);
    }
    else
    {
      axom::for_all<ExecSpace>(nnodes, AXOM_LAMBDA(auto nodeIndex)
      {
        accumFieldView[nodeIndex] = 0;
      });
    }

    // Make the zonal id and vf fields too.
    conduit::Node &n_zonalIDField = n_fields["__equiz_zonalMaterialID"];
    n_zonalIDField["topology"] = n_topo.name();
    n_zonalIDField["association"] = "element";
    n_zonalIDField["values"].set_allocator(c2a.getConduitAllocatorID());
    n_zonalIDField["values"].set(conduit::DataType(bputils::cpp2conduit<MaterialID>::id, nzones));

    conduit::Node &n_zonalVFField = n_fields["__equiz_zonalMaterialVF"];
    n_zonalVFField["topology"] = n_topo.name();
    n_zonalVFField["association"] = "element";
    n_zonalVFField["values"].set_allocator(c2a.getConduitAllocatorID());
    n_zonalVFField["values"].set(conduit::DataType(bputils::cpp2conduit<MaterialVF>::id, nzones));

    auto zonalIDFieldView = bputils::make_array_view<MaterialID>(n_zonalIDField["values"]);
    auto zonalVFFieldView = bputils::make_array_view<MaterialVF>(n_zonalVFField["values"]);
    // Fill all zones with empty.
    axom::for_all<ExecSpace>(nzones, AXOM_LAMBDA(auto nodeIndex)
    {
      zonalIDFieldView[nodeIndex] = NULL_MATERIAL;
      zonalVFFieldView[nodeIndex] = NULL_MATERIAL_VF;
    });

    // Fill in the clean zones.
    using FloatType = typename MatsetView::FloatType;
    MatsetView deviceMatsetView(m_matsetView);
    for(const auto &mat : cleanMats)
    {
      const int matNumber = mat.number;
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(auto zoneIndex) {
          FloatType vf {};
          if(deviceMatsetView.zoneContainsMaterial(zoneIndex, matNumber, vf))
          {
            zonalIDFieldView[zoneIndex] = matNumber;
            zonalVFFieldView[zoneIndex] = static_cast<MaterialVF>(vf);
          }
        });
    }
    // Fill in the mixed zones for the first mixed material.
    if(!mixedMats.empty())
    {
      const int matNumber = mixedMats[0].number;
std::cout << "************ START\n";
#if 1
      const std::string matFieldName = nodalFieldName(matNumber);
      auto matVFView = bputils::make_array_view<MaterialVF>(n_fields.fetch_existing(matFieldName + "/values"));

// NOTE - I think this will cause too many zones to get subdivided.

      // Fill in any zone that has nodes where the nodal matVF is greater than zero.
      m_topologyView.template for_all_zones<ExecSpace>(AXOM_LAMBDA(auto zoneIndex, const auto &zone)
      {
        constexpr MaterialVF NO_MATERIAL = 1.e-6;
        const axom::IndexType n = zone.getIds().size();
        MaterialVF matvfSum {};
        for(axom::IndexType i = 0; i < n; i++)
        {
          matvfSum += matVFView[zone.getId(i)];
        }
        if(matvfSum > NO_MATERIAL)
        {
          zonalIDFieldView[zoneIndex] = matNumber;
          zonalVFFieldView[zoneIndex] = matvfSum / static_cast<MaterialVF>(n);
        }
      });
#else
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(auto zoneIndex) {
          FloatType vf {};
          if(deviceMatsetView.zoneContainsMaterial(zoneIndex, matNumber, vf))
          {
            zonalIDFieldView[zoneIndex] = matNumber;
            zonalVFFieldView[zoneIndex] = static_cast<MaterialVF>(vf);
std::cout << "   storing VF " << zonalVFFieldView[zoneIndex] << "\n";
          }
      });
#endif
std::cout << "************ END\n";
    }
  }

  /**
   * \brief Perform one round of material clipping.
   *
   * \tparam ITopologyView The topology view type for the intermediate topology.
   * \tparam ICoordsetView The topology view type for the intermediate coordset.
   */
  template <typename ITopologyView, typename ICoordsetView>
  void iteration(const ITopologyView &topoView,
                 const ICoordsetView &coordsetView,

                 const axom::mir::views::MaterialInformation &allMats,
                 const axom::mir::views::Material &currentMat,

                 const conduit::Node &n_topo,
                 const conduit::Node &n_coordset,
                 conduit::Node &n_fields,

                 const conduit::Node &n_options,

                 conduit::Node &n_newTopo,
                 conduit::Node &n_newCoordset,
                 conduit::Node &n_newFields)
  {
    namespace bputils = axom::mir::utilities::blueprint;
    namespace bpmeshutils = conduit::blueprint::mesh::utils;
    std::cout << "------------------------ start of iteration "
                 "--------------------------------\n";

    const std::string colorField("__equiz__colors");

#if 1
    //--------------------------------------------------------------------------
    //
    // Save the iteration inputs.
    //
    //--------------------------------------------------------------------------
    conduit::Node n_mesh_input;
    n_mesh_input[n_topo.path()].set_external(n_topo);
    n_mesh_input[n_coordset.path()].set_external(n_coordset);
    n_mesh_input[n_fields.path()].set_external(n_fields);

    // save
    std::stringstream ss1;
    ss1 << "debug_equiz_input." << currentMat.number;
    conduit::relay::io::blueprint::save_mesh(n_mesh_input, ss1.str(), "hdf5");
#endif

    //--------------------------------------------------------------------------
    //
    // Invoke clipping with a material intersector.
    //
    //--------------------------------------------------------------------------
    const std::string matFieldName = nodalFieldName(currentMat.number);

    // Get the accumulated VFs and the VFs for the current material.
    auto accumVFView = bputils::make_array_view<MaterialVF>(n_fields.fetch_existing(accumulatedVFFieldName() + "/values"));
    auto currentMatVFView = bputils::make_array_view<MaterialVF>(n_fields.fetch_existing(matFieldName + "/values"));

    // Make a material intersector.
    using ConnectivityType = typename ITopologyView::ConnectivityType;
    using IntersectorType = MaterialIntersector<ConnectivityType, MatsetView::MaxMaterials>;
#if 1
    IntersectorType intersector;
    // Populate intersector, including making a number:index map
    axom::Array<int> matNumber, matIndex;
    const int nmats = static_cast<int>(allMats.size());
    for(int index = 0; index < nmats; index++)
    {
      // Add a matvf view to the intersector.
      const std::string matFieldName = nodalFieldName(allMats[index].number);
      auto matVFView = bputils::make_array_view<MaterialVF>(n_fields.fetch_existing(matFieldName + "/values"));
      intersector.addMaterial(matVFView);

      matNumber.push_back(allMats[index].number);
      matIndex.push_back(index);
    }
    // Sort indices by matNumber.
    std::sort(matIndex.begin(), matIndex.end(), [&](auto idx1, auto idx2)
    {
      return matNumber[idx1] < matNumber[idx2];
    });
    std::sort(matNumber.begin(), matNumber.end());

    // Store the number:index map into the intersector.
    int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    axom::Array<int> matNumberDevice(nmats, nmats, allocatorID), matIndexDevice(nmats, nmats, allocatorID);
    axom::copy(matNumberDevice.data(), matNumber.data(), sizeof(int) * nmats);
    axom::copy(matIndexDevice.data(), matIndex.data(), sizeof(int) * nmats);
    intersector.setMaterialNumbers(matNumberDevice.view());
    intersector.setMaterialIndices(matIndexDevice.view());

    intersector.setZoneMaterialID(bputils::make_array_view<MaterialID>(n_fields.fetch_existing("__equiz_zonalMaterialID/values")));
    intersector.setCurrentMaterial(currentMat.number);
#else
    IntersectorType intersector(accumVFView, currentMatVFView);
#endif

    // Make clip options.
    conduit::Node options;
    options["inside"] = 1;
    options["outside"] = 1;
    options["colorField"] = colorField;
    if(n_options.has_child("selectedZones"))
    {
      // Pass selectedZones along in the clip options, if present.
      options["selectedZones"].set_external(
        n_options.fetch_existing("selectedZones"));
    }
    options["topology"] = n_options["topology"];

    // Clip the topology using the matset.
    {
      using ClipperType = axom::mir::clipping::ClipField<ExecSpace, ITopologyView, ICoordsetView, IntersectorType>;
      ClipperType clipper(topoView, coordsetView, intersector);
      clipper.execute(n_topo,
                      n_coordset,
                      n_fields,
                      options,
                      n_newTopo,
                      n_newCoordset,
                      n_newFields);
    }

    // In the output fields, add the current material VFs to the accumulated fields.
    // This raises the floor for the next iteration.
    accumVFView = bputils::make_array_view<MaterialVF>(n_newFields.fetch_existing(accumulatedVFFieldName() + "/values"));
    currentMatVFView = bputils::make_array_view<MaterialVF>(n_newFields.fetch_existing(matFieldName + "/values"));
    axom::for_all<ExecSpace>(accumVFView.size(), AXOM_LAMBDA(auto nodeIndex)
    {
      accumVFView[nodeIndex] += currentMatVFView[nodeIndex];
    });

    //--------------------------------------------------------------------------
    //
    // Updated zoneMaterialID and zoneMaterialVF based on color.
    //
    //--------------------------------------------------------------------------
    const auto originalElementsView = bputils::make_array_view<ConnectivityType>(n_newFields.fetch_existing("originalElements/values"));
    const auto colorView = bputils::make_array_view<ConnectivityType>(n_newFields.fetch_existing(colorField + "/values"));
    const auto nzonesNew = originalElementsView.size();

    // Allocate new material ID/VF arrays.
    const int currentMatNumber = currentMat.number;
    conduit::Node &n_zonalMaterialID = n_newFields.fetch_existing("__equiz_zonalMaterialID/values");
    conduit::Node &n_zonalMaterialVF = n_newFields.fetch_existing("__equiz_zonalMaterialVF/values");
    auto zonalMaterialID = bputils::make_array_view<MaterialID>(n_zonalMaterialID);
    auto zonalMaterialVF = bputils::make_array_view<MaterialVF>(n_zonalMaterialVF);
    axom::for_all<ExecSpace>(nzonesNew, AXOM_LAMBDA(auto zoneIndex)
    {
      // Color the part we want with the current material.
      if(colorView[zoneIndex] == 1)
      {
        zonalMaterialID[zoneIndex] = currentMatNumber;

        // NOTE: I'm not 100% sure whether I need to track the VF.

        const auto origIndex = originalElementsView[zoneIndex];
        typename MatsetView::FloatType vf;
        m_matsetView.zoneContainsMaterial(origIndex, currentMatNumber, vf);
        zonalMaterialVF[zoneIndex] = static_cast<MaterialVF>(vf);
      }
    });

#if 1
    //--------------------------------------------------------------------------
    //
    // Save the clip results.
    //
    //--------------------------------------------------------------------------
    conduit::Node mesh;
    mesh[n_newTopo.path()].set_external(n_newTopo);
    mesh[n_newCoordset.path()].set_external(n_newCoordset);
    mesh[n_newFields.path()].set_external(n_newFields);

    // print
    printNode(mesh);
    std::cout.flush();

    // save
    std::stringstream ss;
    ss << "debug_equiz_output." << currentMat.number;
    conduit::relay::io::blueprint::save_mesh(mesh, ss.str(), "hdf5");
#endif

    // We can remove the current material's nodal from the new fields
//    if(n_newFields.has_child(matFieldName))
//    {
//      n_newFields.remove(matFieldName);
//    }
    // We do not want the color field to survive into the next iteration.
    n_newFields.remove(colorField);

    std::cout << "------------------------ end of iteration "
                 "--------------------------------\n";
  }

private:
  TopologyView m_topologyView;
  CoordsetView m_coordsetView;
  MatsetView   m_matsetView;
};

}  // end namespace mir
}  // end namespace axom

#endif
