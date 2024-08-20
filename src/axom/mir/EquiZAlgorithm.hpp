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

#include <algorithm>
#include <string>

#if defined(AXOM_USE_RAJA)
  #include <RAJA/RAJA.hpp>
#endif

// Uncomment to save inputs and outputs.
// #define AXOM_DEBUG_EQUIZ

// This enables a tweak to the algorithm that tries to skip the first iteration
// by incorporating the first material's ids into the zonalMaterialID field. It
// could be faster but it might not be as robust.
// #define AXOM_EQUIZ_SKIP_FIRST_ITERATION

#if defined(AXOM_DEBUG_EQUIZ)
  #include <conduit/conduit_relay_io_blueprint.hpp>
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
constexpr static MaterialVF NULL_MATERIAL_VF = -1.f;

#if 0
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
#endif

namespace detail
{
/**
 * \brief This class is an intersection policy compatible with ClipField. It
 *        helps determine clip cases and weights using material-aware logic.
 *
 * \tparam ConnectivityT The type of index we'd see in the associated mesh's
 *                       connectivity. We template on it so we can pass the
 *                       array views of connectivity (node lists) to methods here.
 * \tparam MAXMATERIALS The max number of materials to handle.
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
    axom::IndexType determineClipCase(axom::IndexType zoneIndex,
                                      const ConnectivityView &nodeIdsView) const
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
        // clangformat-off
        MaterialVF vf1 = (backgroundIndex != INVALID_INDEX) ? m_matvfViews[backgroundIndex][nid] : NULL_MATERIAL_VF;
        MaterialVF vf2 = (currentIndex != INVALID_INDEX) ? m_matvfViews[currentIndex][nid] : 0;
        // clangformat-on

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
    float computeWeight(axom::IndexType zoneIndex,
                        ConnectivityType id0,
                        ConnectivityType id1) const
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
      // clangformat-off
      vf1[0] = (backgroundIndex != INVALID_INDEX) ? m_matvfViews[backgroundIndex][id0] : NULL_MATERIAL_VF;
      vf1[1] = (backgroundIndex != INVALID_INDEX) ? m_matvfViews[backgroundIndex][id1] : NULL_MATERIAL_VF;
      vf2[0] = (currentIndex != INVALID_INDEX) ? m_matvfViews[currentIndex][id0] : 0;
      vf2[1] = (currentIndex != INVALID_INDEX) ? m_matvfViews[currentIndex][id1] : 0;
      // clangformat-on

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

    /**
     * \brief Return the volume fraction array index in m_matIndicesView for the
     *        given material number \a matNumber.
     *
     * \param matNumber A material number that occurs in the matset material ids.
     *
     * \return The m_matNumbersView index on success; INVALID_INDEX on failure.
     */
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

    void setCurrentMaterial(int matNumber) { m_currentMaterial = matNumber; }

    axom::StaticArray<MaterialVFView, MAXMATERIALS>
      m_matvfViews {};  //!< Array of volume fraction views
    axom::ArrayView<int> m_matNumbersView {};  //!< Sorted array of material numbers.
    axom::ArrayView<int> m_matIndicesView {};  //!< Array of indices into m_matvfViews for the material numbers.
    axom::ArrayView<int> m_zoneMatNumberView {};  //!< Contains the current material number that owns each zone.
    int m_currentMaterial {};  //!< The current material.
  };

  /**
   * \brief Initialize the object from options.
   * \param n_options The node that contains the options.
   * \param n_fields The node that contains fields.
   */
  void initialize(const conduit::Node &AXOM_UNUSED_PARAM(n_options),
                  const conduit::Node &AXOM_UNUSED_PARAM(n_fields))
  { }

  /**
   * \brief Determine the name of the topology on which to operate.
   * \param n_input The input mesh node.
   * \param n_options The clipping options.
   * \return The name of the toplogy on which to operate.
   */
  std::string getTopologyName(const conduit::Node &AXOM_UNUSED_PARAM(n_input),
                              const conduit::Node &n_options) const
  {
    return n_options["topology"].as_string();
  }

  /// Set various attributes.

  void addMaterial(const MaterialVFView &matvf) { m_view.addMaterial(matvf); }

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

}  // end namespace detail

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
#if defined(AXOM_DEBUG_EQUIZ)
  void printNode(const conduit::Node &n) const
  {
    conduit::Node options;
    options["num_children_threshold"] = 10000;
    options["num_elements_threshold"] = 10000;
    n.to_summary_string_stream(std::cout, options);
  }
#endif

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
    AXOM_ANNOTATE_SCOPE("EquizAlgorithm");

    // Copy the options.
    conduit::Node n_options_copy;
    bputils::copy<ExecSpace>(n_options_copy, n_options);
    n_options_copy["topology"] = n_topo.name();

#if defined(AXOM_DEBUG_EQUIZ)
    //--------------------------------------------------------------------------
    //
    // Save the MIR input.
    //
    //--------------------------------------------------------------------------
    conduit::Node n_tmpInput;
    n_tmpInput[n_topo.path()].set_external(n_topo);
    n_tmpInput[n_coordset.path()].set_external(n_coordset);
    n_tmpInput[n_fields.path()].set_external(n_fields);
    n_tmpInput[n_matset.path()].set_external(n_matset);
    conduit::relay::io::blueprint::save_mesh(n_tmpInput,
                                             "debug_equiz_input",
                                             "hdf5");
#endif

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

    //--------------------------------------------------------------------------
    //
    // Iterate over mixed materials.
    //
    //--------------------------------------------------------------------------
#if defined(AXOM_EQUIZ_SKIP_FIRST_ITERATION)
    constexpr int first = 1;
#else
    constexpr int first = 0;
#endif
    for(size_t i = first; i < mixedMats.size(); i++)
    {
      if(i == first)
      {
        // The first time through, we can use the supplied views.
        iteration<TopologyView, CoordsetView>(i,
                                              m_topologyView,
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
        // Dispatch to an appropriate topo view.
        // clangformat-off
        views::dispatch_explicit_coordset(n_InputCoordset, [&](auto coordsetView) {
          using ICoordsetView = decltype(coordsetView);
          using ConnectivityType = typename TopologyView::ConnectivityType;
          // Dispatch to an appropriate topo view, taking into account the connectivity
          // type and the possible shapes that would be supported for the input topology.
          views::typed_dispatch_unstructured_topology<
            ConnectivityType,
            views::view_traits<TopologyView>::selected_shapes()>(
            n_InputTopo,
            [&](const auto &AXOM_UNUSED_PARAM(shape), auto topologyView) {
              using ITopologyView = decltype(topologyView);

              // Do the next iteration.
              iteration<ITopologyView, ICoordsetView>(i,
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

    // Build the new matset.
    buildNewMatset(n_matset, n_newFields, n_newMatset);

    // Cleanup.
    {
      AXOM_ANNOTATE_SCOPE("cleanup");
      for(const auto &mat : allMats)
      {
        const std::string nodalMatName(nodalFieldName(mat.number));
        if(n_newFields.has_child(nodalMatName))
        {
          n_newFields.remove(nodalMatName);
        }
#if defined(AXOM_DEBUG_EQUIZ)
        const std::string zonalMatName(zonalFieldName(mat.number));
        if(n_newFields.has_child(zonalMatName))
        {
          n_newFields.remove(zonalMatName);
        }
#endif
      }
      n_newFields.remove(zonalMaterialIDName());
    }

#if defined(AXOM_DEBUG_EQUIZ)
    //--------------------------------------------------------------------------
    //
    // Save the MIR output.
    //
    //--------------------------------------------------------------------------
    conduit::Node n_output;
    n_output[n_newTopo.path()].set_external(n_newTopo);
    n_output[n_newCoordset.path()].set_external(n_newCoordset);
    n_output[n_newFields.path()].set_external(n_newFields);
    n_output[n_newMatset.path()].set_external(n_newMatset);
    conduit::relay::io::blueprint::save_mesh(n_output,
                                             "debug_equiz_output",
                                             "hdf5");
    //printNode(n_output);
#endif
  }

  /**
   * \brief Examine the materials and determine which are clean/mixed.
   *
   * \param n_matset A Conduit node containing the matset.
   * \param[out] allMats A vector of all of the materials.
   * \param[out] cleanMats A vector of the clean materials.
   * \param[out] mixedMats A vector of the mixed materials.
   */
  void classifyMaterials(const conduit::Node &n_matset,
                         axom::mir::views::MaterialInformation &allMats,
                         axom::mir::views::MaterialInformation &cleanMats,
                         axom::mir::views::MaterialInformation &mixedMats) const
  {
    AXOM_ANNOTATE_SCOPE("classifyMaterials");

    cleanMats.clear();
    mixedMats.clear();
    allMats = axom::mir::views::materials(n_matset);

    // TODO: actually determine which materials are clean/mixed. It's probably
    //       best to ask the matsetView since it takes some work to determine
    //       this.

    mixedMats = allMats;
  }

  /**
   * \brief Return the name of the zonal material field for a given matId.
   * \return The name of the zonal material field.
   */
  std::string zonalFieldName(int matId) const
  {
    std::stringstream ss;
    ss << "__equiz_zonal_volume_fraction_" << matId;
    return ss.str();
  }

  /**
   * \brief Return the name of the nodal material field for a given matId.
   * \return The name of the nodal material field.
   */
  std::string nodalFieldName(int matId) const
  {
    std::stringstream ss;
    ss << "__equiz_nodal_volume_fraction_" << matId;
    return ss.str();
  }

  /**
   * \brief Return the name of the zonal material id field.
   * \return The name of the zonal material id field.
   */
  std::string zonalMaterialIDName() const { return "__equiz_zonalMaterialID"; }

  /**
   * \brief Makes node-cenetered volume fractions for the materials in the matset
   *        and attaches them as fields.
   *
   * \param n_topo A Conduit node containing the input topology.
   * \param n_coordset A Conduit node containin the input coordset.
   * \param[inout] A Conduit node where the new fields will be added.
   * \param mixedMats A vector of mixed materials.
   */
  void makeNodeCenteredVFs(const conduit::Node &n_topo,
                           const conduit::Node &n_coordset,
                           conduit::Node &n_fields,
                           const axom::mir::views::MaterialInformation &mixedMats) const
  {
    AXOM_ANNOTATE_SCOPE("makeNodeCenteredVFs");

    namespace bputils = axom::mir::utilities::blueprint;
    // Make a node to zone relation so we know for each node, which zones it touches.
    conduit::Node relation;
    {
      bputils::NodeToZoneRelationBuilder<ExecSpace> rb;
      rb.execute(n_topo, n_coordset, relation);
      //printNode(relation);
      //std::cout.flush();
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
      n_zonalField["values"].set(
        conduit::DataType(bputils::cpp2conduit<MaterialVF>::id, nzones));
      auto zonalFieldView =
        bputils::make_array_view<MaterialVF>(n_zonalField["values"]);

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
      n_nodalField["values"].set(
        conduit::DataType(bputils::cpp2conduit<MaterialVF>::id, nnodes));
      bputils::RecenterField<ExecSpace> z2n;
      z2n.execute(n_zonalField, relation, n_nodalField);

#if !defined(AXOM_DEBUG_EQUIZ)
      // Remove the zonal field that we don't normally need (unless we're debugging).
      n_fields.remove(zonalName);
#endif
    }
  }

  /**
   * \brief Set up the "working fields", mainly a zonalMaterialID that includes
   *        the contributions from the clean materials and the first mixed material.
   *
   * \param n_topo A Conduit node containing the input topology pre-MIR.
   * \param n_fields A Conduit node containing the fields pre-MIR.
   * \param cleanMats A vector of clean materials.
   * \param mixedMats A vector of mixed materials.
   */
  void makeWorkingFields(const conduit::Node &n_topo,
                         conduit::Node &n_fields,
                         const axom::mir::views::MaterialInformation &cleanMats,
#if defined(AXOM_EQUIZ_SKIP_FIRST_ITERATION)
                         const axom::mir::views::MaterialInformation &mixedMats
#else
                         const axom::mir::views::MaterialInformation &AXOM_UNUSED_PARAM(mixedMats)
#endif
                        ) const
  {
    namespace bputils = axom::mir::utilities::blueprint;
    AXOM_ANNOTATE_SCOPE("makeWorkingFields");

    // Get the ID of a Conduit allocator that will allocate through Axom with device allocator allocatorID.
    bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;

    const auto nzones = m_topologyView.numberOfZones();

    // Make the zonal id field.
    conduit::Node &n_zonalIDField = n_fields[zonalMaterialIDName()];
    n_zonalIDField["topology"] = n_topo.name();
    n_zonalIDField["association"] = "element";
    n_zonalIDField["values"].set_allocator(c2a.getConduitAllocatorID());
    n_zonalIDField["values"].set(
      conduit::DataType(bputils::cpp2conduit<MaterialID>::id, nzones));
    auto zonalIDFieldView =
      bputils::make_array_view<MaterialID>(n_zonalIDField["values"]);

    // Fill all zones with NULL_MATERIAL.
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(auto nodeIndex) {
        zonalIDFieldView[nodeIndex] = NULL_MATERIAL;
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
          }
        });
    }
#if defined(AXOM_EQUIZ_SKIP_FIRST_ITERATION)
    // Fill in the mixed zones for the first mixed material.
    if(!mixedMats.empty())
    {
      const int matNumber = mixedMats[0].number;
      const std::string matFieldName = nodalFieldName(matNumber);
      auto matVFView = bputils::make_array_view<MaterialVF>(
        n_fields.fetch_existing(matFieldName + "/values"));

      // Fill in any zone that has nodes where the nodal matVF is greater than zero.
      // This fuzzes it out to more zones so we get better blending with the next
      // material we try to overlay.
      m_topologyView.template for_all_zones<ExecSpace>(
        AXOM_LAMBDA(auto zoneIndex, const auto &zone) {
          constexpr MaterialVF VOLUME_FRACTION_CUTOFF = 1.e-6;
          MaterialVF matvfSum {};
          for(const auto nid : zone.getIds())
          {
            matvfSum += matVFView[nid];
          }
          // Overwrite the existing material.
          if(matvfSum > VOLUME_FRACTION_CUTOFF)
          {
            zonalIDFieldView[zoneIndex] = matNumber;
          }
        });
    }
#endif
  }

  /**
   * \brief Perform one iteration of material clipping.
   *
   * \tparam ITopologyView The topology view type for the intermediate topology.
   * \tparam ICoordsetView The topology view type for the intermediate coordset.
   *
   * \param iter The iteration number.
   * \param topoView The topology view for the intermediate input topology.
   * \param coordsetView The coordset view for the intermediate input coordset.
   * \param allMats A vector of Material information (all materials).
   * \param currentMat A Material object for the current material.
   * \param n_topo A Conduit node containing the intermediate input topology.
   * \param n_fields A Conduit node containing the intermediate input fields.
   * \param n_options MIR options.
   * \param n_newTopo[out] A Conduit node to contain the new topology.
   * \param n_newCoordset[out] A Conduit node to contain the new coordset.
   * \param n_newFields[out] A Conduit node to contain the new fields.
   *
   * \note This algorithm uses a ClipField with a MaterialIntersector that gives
   *       it the ability to access nodal volume fraction fields and make intersection
   *       decisions with that data.
   */
  template <typename ITopologyView, typename ICoordsetView>
  void iteration(int iter,
                 const ITopologyView &topoView,
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
    AXOM_ANNOTATE_SCOPE(axom::fmt::format("iteration {}", iter));

    const std::string colorField("__equiz__colors");

#if defined(AXOM_DEBUG_EQUIZ)
    //--------------------------------------------------------------------------
    //
    // Save the iteration inputs.
    //
    //--------------------------------------------------------------------------
    {
      AXOM_ANNOTATE_SCOPE("Saving input");
      conduit::Node n_mesh_input;
      n_mesh_input[n_topo.path()].set_external(n_topo);
      n_mesh_input[n_coordset.path()].set_external(n_coordset);
      n_mesh_input[n_fields.path()].set_external(n_fields);

      // save
      std::stringstream ss1;
      ss1 << "debug_equiz_input_iter." << iter;
      conduit::relay::io::blueprint::save_mesh(n_mesh_input, ss1.str(), "hdf5");
    }
#endif

    //--------------------------------------------------------------------------
    //
    // Make material intersector.
    //
    //--------------------------------------------------------------------------
    using ConnectivityType = typename ITopologyView::ConnectivityType;
    using IntersectorType =
      detail::MaterialIntersector<ConnectivityType, MatsetView::MaxMaterials>;

    IntersectorType intersector;
    int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    const int nmats = static_cast<int>(allMats.size());
    axom::Array<int> matNumberDevice(nmats, nmats, allocatorID),
      matIndexDevice(nmats, nmats, allocatorID);
    {
      AXOM_ANNOTATE_SCOPE("Intersector setup");
      // Populate intersector, including making a number:index map
      axom::Array<int> matNumber, matIndex;
      for(int index = 0; index < nmats; index++)
      {
        // Add a matvf view to the intersector.
        const std::string matFieldName = nodalFieldName(allMats[index].number);
        auto matVFView = bputils::make_array_view<MaterialVF>(
          n_fields.fetch_existing(matFieldName + "/values"));
        intersector.addMaterial(matVFView);

        matNumber.push_back(allMats[index].number);
        matIndex.push_back(index);
      }
      // Sort indices by matNumber.
      std::sort(matIndex.begin(), matIndex.end(), [&](auto idx1, auto idx2) {
        return matNumber[idx1] < matNumber[idx2];
      });
      std::sort(matNumber.begin(), matNumber.end());

      // Store the number:index map into the intersector. The number:index map lets us
      // ask for the field index for a material number, allowing scattered material
      // numbers to be used in the matset.
      axom::copy(matNumberDevice.data(), matNumber.data(), sizeof(int) * nmats);
      axom::copy(matIndexDevice.data(), matIndex.data(), sizeof(int) * nmats);
      intersector.setMaterialNumbers(matNumberDevice.view());
      intersector.setMaterialIndices(matIndexDevice.view());

      // Store the current zone material ids and current material number into the intersector.
      intersector.setZoneMaterialID(bputils::make_array_view<MaterialID>(
        n_fields.fetch_existing(zonalMaterialIDName() + "/values")));
      intersector.setCurrentMaterial(currentMat.number);
    }

    //--------------------------------------------------------------------------
    //
    // Make clip options
    //
    //--------------------------------------------------------------------------
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

    //--------------------------------------------------------------------------
    //
    // Clip the topology using the material intersector.
    //
    //--------------------------------------------------------------------------
    {
      AXOM_ANNOTATE_SCOPE("Clipping");

      using ClipperType =
        axom::mir::clipping::ClipField<ExecSpace, ITopologyView, ICoordsetView, IntersectorType>;
      ClipperType clipper(topoView, coordsetView, intersector);
      clipper.execute(n_topo,
                      n_coordset,
                      n_fields,
                      options,
                      n_newTopo,
                      n_newCoordset,
                      n_newFields);
    }

    //--------------------------------------------------------------------------
    //
    // Update zoneMaterialID based on color field.
    //
    //--------------------------------------------------------------------------
    {
      AXOM_ANNOTATE_SCOPE("Update zonalMaterialID");

      const auto colorView = bputils::make_array_view<ConnectivityType>(
        n_newFields.fetch_existing(colorField + "/values"));
      const auto nzonesNew = colorView.size();

      // Get zonalMAterialID field so we can make adjustments.
      conduit::Node &n_zonalMaterialID =
        n_newFields.fetch_existing(zonalMaterialIDName() + "/values");
      auto zonalMaterialID =
        bputils::make_array_view<MaterialID>(n_zonalMaterialID);
      const int currentMatNumber = currentMat.number;
      axom::for_all<ExecSpace>(
        nzonesNew,
        AXOM_LAMBDA(auto zoneIndex) {
          // Color the part we want with the current material.
          if(colorView[zoneIndex] == 1)
          {
            zonalMaterialID[zoneIndex] = currentMatNumber;
          }
        });
    }

#if defined(AXOM_DEBUG_EQUIZ)
    //--------------------------------------------------------------------------
    //
    // Save the clip results.
    //
    //--------------------------------------------------------------------------
    {
      AXOM_ANNOTATE_SCOPE("Saving output");
      conduit::Node mesh;
      mesh[n_newTopo.path()].set_external(n_newTopo);
      mesh[n_newCoordset.path()].set_external(n_newCoordset);
      mesh[n_newFields.path()].set_external(n_newFields);

      // print
      //printNode(mesh);
      //std::cout.flush();

      // save
      std::stringstream ss;
      ss << "debug_equiz_output_iter." << iter;
      conduit::relay::io::blueprint::save_mesh(mesh, ss.str(), "hdf5");
    }
#endif

    // We do not want the color field to survive into the next iteration.
    n_newFields.remove(colorField);
  }

  /**
   * \brief Build a new matset with only clean zones, representing the MIR output.
   *
   * \param n_matset n_matset The Conduit node that contains the input matset.
   * \param[inout] n_newFields The Conduit node that contains the fields for the MIR output.
   * \param[out] n_newMatset The node that contains the new matset.
   */
  void buildNewMatset(const conduit::Node &n_matset,
                      conduit::Node &n_newFields,
                      conduit::Node &n_newMatset) const
  {
    namespace bputils = axom::mir::utilities::blueprint;
    AXOM_ANNOTATE_SCOPE("buildNewMatset");

    // Get the zonalMaterialID field that has our new material ids.
    conduit::Node &n_zonalMaterialID =
      n_newFields[zonalMaterialIDName() + "/values"];
    auto zonalMaterialID =
      bputils::make_array_view<MaterialID>(n_zonalMaterialID);
    const auto nzones = n_zonalMaterialID.dtype().number_of_elements();

    // Copy some information from the old matset to the new one.
    if(n_matset.has_child("topology"))
    {
      n_newMatset["topology"].set(n_matset.fetch_existing("topology"));
    }
    if(n_matset.has_child("material_map"))
    {
      n_newMatset["material_map"].set(n_matset.fetch_existing("material_map"));
    }

    // Make new nodes in the matset.
    conduit::Node &n_material_ids = n_newMatset["material_ids"];
    conduit::Node &n_volume_fractions = n_newMatset["volume_fractions"];
    conduit::Node &n_sizes = n_newMatset["sizes"];
    conduit::Node &n_offsets = n_newMatset["offsets"];
    conduit::Node &n_indices = n_newMatset["indices"];

    bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;
    n_material_ids.set_allocator(c2a.getConduitAllocatorID());
    n_volume_fractions.set_allocator(c2a.getConduitAllocatorID());
    n_sizes.set_allocator(c2a.getConduitAllocatorID());
    n_offsets.set_allocator(c2a.getConduitAllocatorID());
    n_indices.set_allocator(c2a.getConduitAllocatorID());

    // We'll store the output matset in the same types as the input matset.
    using MIntType = typename MatsetView::IndexType;
    using MFloatType = typename MatsetView::FloatType;
    n_material_ids.set(
      conduit::DataType(bputils::cpp2conduit<MIntType>::id, nzones));
    n_volume_fractions.set(
      conduit::DataType(bputils::cpp2conduit<MFloatType>::id, nzones));
    n_sizes.set(conduit::DataType(bputils::cpp2conduit<MIntType>::id, nzones));
    n_offsets.set(conduit::DataType(bputils::cpp2conduit<MIntType>::id, nzones));
    n_indices.set(conduit::DataType(bputils::cpp2conduit<MIntType>::id, nzones));

    auto material_ids_view = bputils::make_array_view<MIntType>(n_material_ids);
    auto volume_fractions_view =
      bputils::make_array_view<MFloatType>(n_volume_fractions);
    auto sizes_view = bputils::make_array_view<MIntType>(n_sizes);
    auto offsets_view = bputils::make_array_view<MIntType>(n_offsets);
    auto indices_view = bputils::make_array_view<MIntType>(n_indices);

    // Fill in the new matset data arrays.
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(auto zoneIndex) {
        material_ids_view[zoneIndex] =
          static_cast<MIntType>(zonalMaterialID[zoneIndex]);
        volume_fractions_view[zoneIndex] = 1;
        sizes_view[zoneIndex] = 1;
        offsets_view[zoneIndex] = static_cast<MIntType>(zoneIndex);
        indices_view[zoneIndex] = static_cast<MIntType>(zoneIndex);
      });
  }

private:
  TopologyView m_topologyView;
  CoordsetView m_coordsetView;
  MatsetView m_matsetView;
};

}  // end namespace mir
}  // end namespace axom

#endif
