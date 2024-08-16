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
using MaterialIDArray = axom::Array<int>;
using MaterialIDView = axom::ArrayView<int>;
using MaterialVF = float;
using MaterialVFArray = axom::Array<MaterialVF>;
using MaterialVFView = axom::ArrayView<MaterialVF>;

constexpr static int NULL_MATERIAL = -1;
constexpr static MaterialVF NULL_MATERIAL_VF = 0.f;

/**
 * \brief Populate values for a new field based on a material's volume fraction.
 */
template <typename ExecSpace, typename MatsetView, typename OriginalElementsView>
struct MatsetToField
{
  using MaterialIndex = typename MatsetView::MaterialIndex;
  using FloatType = typename MatsetView::FloatType;

  void execute(const MatsetView &matsetView,
               MaterialIndex mat,
               MaterialVFView vfView)
  {
    const auto nzones = vfView.size();
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(auto zoneIndex) {
        FloatType vf {};
        matsetView.zoneContainsMaterial(zoneIndex, mat, vf);
        vfView[zoneIndex] = static_cast<MaterialVF>(vf);
      });
  }

  void execute(const MatsetView &matsetView,
               MaterialIndex mat,
               OriginalElementsView origElementsView,
               MaterialVFView vfView)
  {
    const auto nzones = vfView.size();
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(auto zoneIndex) {
        FloatType vf {};
        const auto origIndex = origElementsView[zoneIndex];
        matsetView.zoneContainsMaterial(origIndex, mat, vf);
        vfView[zoneIndex] = static_cast<MaterialVF>(vf);
      });
  }
};

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

/**
 * \brief This class is an intersection policy compatible with ClipField. It
 *        helps determine clip cases and weights using material-aware logic.
 *
 * \tparam MAXELEMENTS the max number of nodes in any zone we'll encounter.
 */
template <typename ConnectivityT, typename RelationViewType, int MAXELEMENTS = 8>
class MaterialIntersector
{
public:
  using ConnectivityType = ConnectivityT;
  using ConnectivityView = axom::ArrayView<ConnectivityType>;
  using RelationView = RelationViewType;

  MaterialIntersector(const RelationView &rv, const MaterialIDView &zoneMatIds, const MaterialVFView &zoneVFs, const MaterialVFView &nodalCurrentVF) : 
    m_view(rv, zoneMatIds, zoneVFs, nodalCurrentVF)
  { }

  /**
   * \brief This is a view class for MatsetIntersector that can be used in device code.
   */
  struct View
  {
    using VFList = StaticArray<float, MAXELEMENTS>;

    /// Constructor
    View() = default;

    /// Constructor
    View(const RelationView &rv, const MaterialIDView &zoneMatIds, const MaterialVFView &zoneVFs, const MaterialVFView &nodalCurrentVF) :
      m_relationView(rv), m_zoneMaterialIDView(zoneMatIds), m_zoneMaterialVFView(zoneVFs),
      m_nodalCurrentMaterialVFView(nodalCurrentVF)
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
    axom::IndexType determineClipCase(axom::IndexType zoneIndex, const ConnectivityView &nodeIdsView) const
    {
      VFList vf1, vf2;
      getVolumeFractionsAtNodes(m_zoneMaterialIDView[zoneIndex], nodeIdsView, vf1);
      getVolumeFractionsCurrent(nodeIdsView, vf2);

      axom::IndexType clipcase = 0;
      const auto n = nodeIdsView.size();
      for(IndexType i = 0; i < n; i++)
      {
        clipcase |= (vf2[i] > vf1[i]) ? (1 << i) : 0;
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
      // Get the volume fractions for mat1, mat2 at the edge endpoints id0, id1.
      ConnectivityType nodeIds[2] = {id0, id1};
      ConnectivityView nodeIdsView(nodeIds, 2);
      VFList vf1, vf2;
      getVolumeFractionsAtNodes(m_zoneMaterialIDView[zoneIndex], nodeIdsView, vf1);
      getVolumeFractionsCurrent(nodeIdsView, vf2);

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
     * \brief Return node-averaged volume fractions for the specified material at the given node ids.
     *
     * \param materialNumber The material number for which we want nodal volume fractions.
     * \param nodeIds A view containing the node ids where we want nodal volume fractions.
     * \param[out] vfs A list of output volume fractions.
     */
    AXOM_HOST_DEVICE
    void getVolumeFractionsAtNodes(int materialNumber, const ConnectivityView &nodeIds, VFList &vfs) const
    {
      for(axom::IndexType i = 0; i < nodeIds.size(); i++)
      {
        const auto nid = nodeIds[i];

        // Use the relation to average VFs for this material to the node.
        auto size = m_relationView.m_sizesView[nid];
        using SizeType = decltype(size);
        const auto offset = m_relationView.m_offsetsView[nid];

        MaterialVF vfSum{};
        for(SizeType j = 0; j < size; j++)
        {
          const auto zoneIndex = m_relationView.m_zonesView[offset + j];
          vfSum += (m_zoneMaterialIDView[zoneIndex] == materialNumber) ? m_zoneMaterialVFView[zoneIndex] : 0.f;
        }

        const MaterialVF nodeVF = vfSum / static_cast<MaterialVF>(size);
        vfs.push_back(nodeVF);
      }
    }

    AXOM_HOST_DEVICE
    void getVolumeFractionsCurrent(const ConnectivityView &nodeIds, VFList &vfs) const
    {
      for(axom::IndexType i = 0; i < nodeIds.size(); i++)
      {
        vfs.push_back(m_nodalCurrentMaterialVFView[nodeIds[i]]);
      }
    }

    RelationView m_relationView {};

    MaterialIDView m_zoneMaterialIDView {};         //!< Zonal field of material ids added so far.
    MaterialVFView m_zoneMaterialVFView {};         //!< Zonal vf for the zone inherited from the material

    MaterialVFView m_nodalCurrentMaterialVFView {}; //!< Nodal field with current material VFs
    int            m_currentMaterial {0};           //!< Id of the current material
  };

  /**
   * \brief Initialize the object from options.
   * \param n_options The node that contains the options.
   * \param n_fields The node that contains fields.
   */
  void initialize(const conduit::Node &n_options, const conduit::Node &AXOM_UNUSED_PARAM(n_fields))
  {
    // We pass in the current material from the MIR algorithm via the options.
    if(n_options.has_path("currentMaterial"))
      m_view.m_currentMaterial = n_options["currentMaterial"].to_int();
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
  void printNode(const conduit::Node &n)
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
    conduit::Node &n_InputMatset = n_Input[n_matset.path()];

    // Make "empty" material data, kind of unofficial fields on the input mesh.
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    const auto nzones = m_topologyView.numberOfZones();
    axom::Array<int> zoneMaterialID(nzones, nzones, allocatorID);
    axom::Array<MaterialVF> zoneMaterialVF(nzones, nzones, allocatorID);
    auto zoneMaterialIDView = zoneMaterialID.view();
    auto zoneMaterialVFView = zoneMaterialVF.view();
    axom::for_all<ExecSpace>(nzones, AXOM_LAMBDA(auto zoneIndex)
    {
      zoneMaterialIDView[zoneIndex] = NULL_MATERIAL;
      zoneMaterialVFView[zoneIndex] = NULL_MATERIAL_VF;
    });

    // Iterate over the materials
    const auto matInfo = axom::mir::views::materials(n_matset);
    for(size_t i = 0; i < matInfo.size(); i++)
    {
      if(i == 0)
      {
#if 1
        // Print input mesh.
        std::cout << "-------------------------------------------------------------------" << std::endl;
        conduit::Node mesh;
        mesh[n_topo.path()].set_external(n_topo);
        mesh[n_coordset.path()].set_external(n_coordset);
        mesh[n_fields.path()].set_external(n_fields);
        mesh[n_matset.path()].set_external(n_matset);
        printNode(mesh);
        std::cout << "-------------------------------------------------------------------" << std::endl;
#endif
        // The first time through, we can use the supplied views.
        iteration<TopologyView, CoordsetView>(m_topologyView,
                                              m_coordsetView,

                                              zoneMaterialID,
                                              zoneMaterialVF,
                                              matInfo[i],

                                              n_topo,
                                              n_coordset,
                                              n_fields,

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
        // Move the outputs of the last iteration to the inputs of this iteration.
        n_InputTopo.move(n_newTopo);
        n_InputCoordset.move(n_newCoordset);
        n_InputFields.move(n_newFields);
        n_InputMatset.move(n_newMatset);

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

                zoneMaterialID,
                zoneMaterialVF,
                matInfo[i],

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

  /**
   * \brief Perform one round of material clipping.
   *
   * \tparam ITopologyView The topology view type for the intermediate topology.
   * \tparam ICoordsetView The topology view type for the intermediate coordset.
   */
  template <typename ITopologyView, typename ICoordsetView>
  void iteration(const ITopologyView &topoView,
                 const ICoordsetView &coordsetView,

                 MaterialIDArray &zoneMaterialID,
                 MaterialVFArray &zoneMaterialVF,

                 const axom::mir::views::Material &currentMat,

                 const conduit::Node &n_topo,
                 const conduit::Node &n_coordset,
                 const conduit::Node &n_fields,

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
    const std::string nodalFieldName("__equiz__nodal_current_material");
    const std::string zonalFieldName("__equiz__zonal_current_material");

    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    // Get the ID of a Conduit allocator that will allocate through Axom with device allocator allocatorID.
    utilities::blueprint::ConduitAllocateThroughAxom<ExecSpace> c2a;

    // Make a node to zone relation.
    conduit::Node relation;
    {
      bputils::NodeToZoneRelationBuilder<ExecSpace> rb;
      rb.execute(n_topo, n_coordset, relation);
      printNode(relation);
      std::cout.flush();
    }

    const auto nzones = topoView.numberOfZones();
    const auto nnodes = coordsetView.numberOfNodes();

    //--------------------------------------------------------------------------
    //
    // Make input fields we can modify
    //
    //--------------------------------------------------------------------------
    conduit::Node n_inputFields;
    for(conduit::index_t i = 0; i < n_fields.number_of_children(); i++)
    {
      n_inputFields[n_fields[i].name()].set_external(n_fields[i]);
    }
    
    //--------------------------------------------------------------------------
    //
    // Make a zonal field on the input mesh for the current material, pulling
    // out its VFs from the original matsetView.
    //
    // NOTE: The fields are free-floating because we just are making array views
    //       using that data.
    //--------------------------------------------------------------------------
    conduit::Node &n_zonalField = n_inputFields[zonalFieldName];
    n_zonalField["topology"] = n_topo.name();
    n_zonalField["association"] = "element";
    n_zonalField["values"].set_allocator(c2a.getConduitAllocatorID());
    n_zonalField["values"].set(conduit::DataType(bputils::cpp2conduit<MaterialVF>::id, nzones));
    auto currentZonalMaterialVFView = bputils::make_array_view<MaterialVF>(n_zonalField["values"]);

    using ConnectivityType = typename ITopologyView::ConnectivityType;
    using OriginalElementsView = axom::ArrayView<ConnectivityType>;

    {
      MatsetToField<ExecSpace, MatsetView, OriginalElementsView> m2f;
      if(n_fields.has_path("originalElements"))
      {
        auto originalElementsView = bputils::make_array_view<ConnectivityType>(n_fields["originalElements/values"]);
        m2f.execute(m_matsetView, currentMat.number, originalElementsView, currentZonalMaterialVFView);
      }
      else
      {
        m2f.execute(m_matsetView, currentMat.number, currentZonalMaterialVFView);
      }
    }

    // Make a nodal field for the current material by recentering.
    conduit::Node &n_nodalField = n_inputFields[nodalFieldName];
    {
      n_nodalField["topology"] = n_topo.name();
      n_nodalField["association"] = "vertex";
      n_nodalField["values"].set_allocator(c2a.getConduitAllocatorID());
      n_nodalField["values"].set(conduit::DataType(bputils::cpp2conduit<MaterialVF>::id, nnodes));
      bputils::RecenterField<ExecSpace> z2n;
      z2n.execute(n_zonalField, relation, n_nodalField);
    }
    auto currentNodalMaterialVFView = bputils::make_array_view<MaterialVF>(n_nodalField["values"]);

#if 1
    //--------------------------------------------------------------------------
    //
    // Save the iteration inputs.
    //
    //--------------------------------------------------------------------------
    conduit::Node n_mesh_input;
    n_mesh_input[n_topo.path()].set_external(n_topo);
    n_mesh_input[n_coordset.path()].set_external(n_coordset);
    n_mesh_input[n_fields.path()].set_external(n_inputFields);

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
    const auto zoneMaterialIDView = zoneMaterialID.view();
    const auto zoneMaterialVFView = zoneMaterialVF.view();
    axom::mir::views::IndexNode_to_ArrayView_same(relation["zones"], relation["sizes"], relation["offsets"],
      [&](auto zonesView, auto sizesView, auto offsetsView)
    {
      // Wrap the node to zone relation in a view.
      using IndexT = typename decltype(zonesView)::value_type;
      using RelationViewType = RelationView<IndexT>;
      RelationViewType rv;
      rv.m_zonesView = zonesView;
      rv.m_sizesView = sizesView;
      rv.m_offsetsView = offsetsView;

      // Make a material intersector.
      using IntersectorType = MaterialIntersector<ConnectivityType, RelationViewType>;
      IntersectorType intersector(rv, zoneMaterialIDView, zoneMaterialVFView, currentNodalMaterialVFView);

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
      // Tell the intersector which material we're currently working on.
      options["currentMaterial"] = currentMat.number;
      options["topology"] = n_options["topology"];

      // Clip the topology using the matset.
      {
        using ClipperType = axom::mir::clipping::ClipField<ExecSpace, ITopologyView, ICoordsetView, IntersectorType>;
        ClipperType clipper(topoView, coordsetView, intersector);
        clipper.execute(n_topo,
                        n_coordset,
                        n_inputFields,
                        options,
                        n_newTopo,
                        n_newCoordset,
                        n_newFields);
      }
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
    MaterialIDArray newZoneMaterialID(nzonesNew, nzonesNew, allocatorID);
    MaterialVFArray newZoneMaterialVF(nzonesNew, nzonesNew, allocatorID);
    auto newZoneMaterialIDView = newZoneMaterialID.view();
    auto newZoneMaterialVFView = newZoneMaterialVF.view();
    const int currentMatNumber = currentMat.number;
    axom::for_all<ExecSpace>(nzonesNew, AXOM_LAMBDA(auto zoneIndex)
    {
      const auto origIndex = originalElementsView[zoneIndex];
      if(colorView[zoneIndex] == 1)
      {
        // Color 1 - This zone should be the current material. We can grab its
        //           volume fraction from the original matset.
        newZoneMaterialIDView[zoneIndex] = currentMatNumber;
        typename MatsetView::FloatType vf;
        m_matsetView.zoneContainsMaterial(origIndex, currentMatNumber, vf);
        newZoneMaterialVFView[zoneIndex] = static_cast<MaterialVF>(vf);
      }
      else
      {
        // Color 0 - This zone was not the current material. Pass through what was there before.
        newZoneMaterialIDView[zoneIndex] = zoneMaterialIDView[origIndex];
        newZoneMaterialVFView[zoneIndex] = zoneMaterialVFView[origIndex];
      }
    });
    // Pass the new data arrays out.
    zoneMaterialID.swap(newZoneMaterialID);
    zoneMaterialVF.swap(newZoneMaterialVF);

#if 1
    //--------------------------------------------------------------------------
    //
    // Save the clip results.
    //
    //--------------------------------------------------------------------------

#define SAVE_ZONE_MATERIALS
#ifdef SAVE_ZONE_MATERIALS
    // Attach the updated arrays to the new fields so we can look at them.
    n_newFields["zoneMaterialID/association"] = "element";
    n_newFields["zoneMaterialID/topology"] = n_topo.name();
    n_newFields["zoneMaterialID/values"].set_external(zoneMaterialID.data(), zoneMaterialID.size());

    n_newFields["zoneMaterialVF/association"] = "element";
    n_newFields["zoneMaterialVF/topology"] = n_topo.name();
    n_newFields["zoneMaterialVF/values"].set_external(zoneMaterialVF.data(), zoneMaterialVF.size());
#endif
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

#ifdef SAVE_ZONE_MATERIALS
    n_newFields.remove("zoneMaterialID");
    n_newFields.remove("zoneMaterialVF");
#endif
#endif

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
