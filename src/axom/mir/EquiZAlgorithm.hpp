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
/**
 * \brief Populate values for a new field based on a material's volume fraction.
 */
template <typename ExecSpace, typename MatsetView>
struct MatsetToField
{
  using MaterialIndex = typename MatsetView::MaterialIndex;
  using FloatType = typename MatsetView::FloatType;

  void execute(const MatsetView &matsetView,
               MaterialIndex mat,
               axom::ArrayView<FloatType> &vfValues)
  {
    const auto nzones = vfValues.size();
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(auto zoneIndex) {
        FloatType vf {};
        matsetView.zoneContainsMaterial(zoneIndex, mat, vf);
        vfValues[zoneIndex] = vf;
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
 * \tparam MAXELEMENTS the max number of nodes in any zone we'll encounter.
 */
template <typename ConnectivityT, typename MatsetViewType, typename RelationViewType, int MAXELEMENTS = 8>
class MatsetIntersector
{
public:
  using ConnectivityType = ConnectivityT;
  using ConnectivityView = axom::ArrayView<ConnectivityType>;
  using MatsetView = MatsetViewType;
  using RelationView = RelationViewType;

  MatsetIntersector(const MatsetView &mv, const RelationView &rv) : m_view(mv, rv)
  { }

  /**
   * \brief This is a view class for MatsetIntersector that can be used in device code.
   */
  struct View
  {
    using FloatType = typename MatsetView::FloatType;
    using VFList = StaticArray<typename MatsetView::FloatType, MAXELEMENTS>;
//    constexpr static ConnectivityType NULL_MAT = -1;

    View() = default;

    View(const MatsetView &mv, const RelationView &rv) : m_matsetView(mv), m_relationView(rv)
    { }

    AXOM_HOST_DEVICE
    axom::IndexType determineClipCase(axom::IndexType zoneIndex, const ConnectivityView &nodeIds) const
    {
      VFList vf1, vf2;
      const auto dominantMaterial = m_matsetView.backgroundMaterial(zoneIndex, m_currentMaterial); //dominantMaterial(zoneIndex);
      getVolumeFractionsAtNodes(dominantMaterial, m_currentMaterial, nodeIds, vf1, vf2);

      axom::IndexType clipcase = 0;
      const auto n = nodeIds.size();
      for(IndexType i = 0; i < n; i++)
      {
        clipcase |= (vf2[i] > vf1[i]) ? (1 << i) : 0;
      }
#if 0
std::cout << "clipcase: zoneIndex=" << zoneIndex
          << ", dominantMaterial=" << dominantMaterial
          << ", m_currentMaterial=" << m_currentMaterial
          << ", nodeIds={";
      for(IndexType i = 0; i < n; i++)
      {
        std::cout << nodeIds[i] << ", ";
      }
std::cout << "}, vf1={";
      for(IndexType i = 0; i < vf1.size(); i++)
      {
        std::cout << vf1[i] << ", ";
      }
std::cout << "}, vf2={";
      for(IndexType i = 0; i < vf2.size(); i++)
      {
        std::cout << vf2[i] << ", ";
      }
std::cout << "}, clipcase=" << clipcase << std::endl;
#endif
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
//      const auto dominantMaterial = m_matsetView.dominantMaterial(zoneIndex);
      const auto dominantMaterial = m_matsetView.backgroundMaterial(zoneIndex, m_currentMaterial); //dominantMaterial(zoneIndex);
      getVolumeFractionsAtNodes(dominantMaterial, m_currentMaterial, nodeIdsView, vf1, vf2);

      //vf1[0] = axom::utilities::max(vf1[0], 0.f);
      //vf1[1] = axom::utilities::max(vf1[1], 0.f);
      //vf2[0] = axom::utilities::max(vf2[0], 0.f);
      //vf2[1] = axom::utilities::max(vf2[1], 0.f);

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
     */
    AXOM_HOST_DEVICE
    void getVolumeFractionAtNodes(int materialNumber, const ConnectivityView &nodeIds, VFList &vfs) const
    {
//std::cout << "{\ngetVolumeFractionAtNodes: mat=" << materialNumber << std::endl;
      for(axom::IndexType i = 0; i < nodeIds.size(); i++)
      {
        const auto nid = nodeIds[i];

        // Use the relation to average VFs for this material to the node.
        auto size = m_relationView.m_sizesView[nid];
        using SizeType = decltype(size);
        const auto offset = m_relationView.m_offsetsView[nid];
//std::cout << "\tnid=" << nid << ", size=" << size << ", offset=" << offset << std::endl;
        float vfSum{};
        for(SizeType j = 0; j < size; j++)
        {
          const auto zoneIndex = m_relationView.m_zonesView[offset + j];
          float vf{};
          m_matsetView.zoneContainsMaterial(zoneIndex, materialNumber, vf);
          vfSum += vf;
//std::cout <<"\t\tzoneIndex=" << zoneIndex << ", vf=" << vf << ", vfSum=" << vfSum << std::endl;
        }
        const float nodeVF = vfSum / static_cast<float>(size);
//std::cout << "\tnodeVF=" << nodeVF << std::endl; 
        vfs.push_back(nodeVF);
      }
//std::cout << "}" << std::endl;
    }

    AXOM_HOST_DEVICE
    void getVolumeFractionsAtNodes(int mat1, int mat2, const ConnectivityView &nodeIds, VFList &vf1, VFList &vf2) const
    {
#if 1
      getVolumeFractionAtNodes(mat1, nodeIds, vf1);
      getVolumeFractionAtNodes(mat2, nodeIds, vf2);
#else
      if(mat1 == NULL_MAT) // Q: do we need this?
      {
        for(axom::IndexType i = 0; i < nodeIds.size(); i++)
          vf1.push_back(-1);
      }
      else
      {
        getVolumeFractionAtNodes(mat1, nodeIds, vf1);
      }

      if(mat2 == NULL_MAT)
      {
        for(axom::IndexType i = 0; i < nodeIds.size(); i++)
          vf1.push_back(-1);
      }
      else
      {
        getVolumeFractionAtNodes(mat2, nodeIds, vf2);
      }
#endif
    }

    MatsetView m_matsetView {};
    RelationView m_relationView {};
    int          m_currentMaterial {0};
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
std::cout << "MIR getTopologyName: " << n_options["topology"].as_string() << std::endl;
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
        iteration<TopologyView, CoordsetView, MatsetView>(m_topologyView,
                                                          m_coordsetView,
                                                          m_matsetView,
                                                          matInfo[i],
                                                          n_topo,
                                                          n_coordset,
                                                          n_fields,
                                                          n_matset,
                                                          n_options_copy,
                                                          n_newTopo,
                                                          n_newCoordset,
                                                          n_newFields,
                                                          n_newMatset);

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

        const auto shape =
          n_InputTopo.fetch_existing("elements/shape").as_string();
        if(shape == "mixed")
        {
          // The data are now an unstructured view, probably a mixed shape view.
          // Dispatch to an appropriate topo view
          // clangformat-off
          views::dispatch_explicit_coordset(n_InputCoordset, [&](auto coordsetView) {
            using ICoordsetView = decltype(coordsetView);
            views::dispatch_unstructured_mixed_topology(
              n_InputTopo,
              [&](const auto &AXOM_UNUSED_PARAM(shape), auto topologyView) {
                using ITopologyView = decltype(topologyView);

                // The output of the first iteration was a unibuffer matset.
                using IMatsetView =
                  axom::mir::views::UnibufferMaterialView<int, float, 10>;
                IMatsetView matsetView;
                matsetView.set(
                  bputils::make_array_view<int>(n_InputMatset["material_ids"]),
                  bputils::make_array_view<float>(
                    n_InputMatset["volume_fractions"]),
                  bputils::make_array_view<int>(n_InputMatset["sizes"]),
                  bputils::make_array_view<int>(n_InputMatset["offsets"]),
                  bputils::make_array_view<int>(n_InputMatset["indices"]));

                // Do the next iteration.
                iteration<ITopologyView, ICoordsetView, IMatsetView>(
                  topologyView,
                  coordsetView,
                  matsetView,
                  matInfo[i],
                  n_InputTopo,
                  n_InputCoordset,
                  n_InputFields,
                  n_InputMatset,
                  n_options_copy,
                  n_newTopo,
                  n_newCoordset,
                  n_newFields,
                  n_newMatset);
              });
          });
          // clangformat-on
        }
        else
        {
          // NOTE: what if the topo output was all clean after the 1st round?
        }
      }

      // TODO: we need to build the new matset.
    }
  }

  /**
   * \brief Perform one round of material clipping.
   *
   * \tparam ITopologyView The topology view type for the intermediate topology.
   * \tparam ICoordsetView The topology view type for the intermediate coordset.
   * \tparam IMatsetView The topology view type for the intermediate matset.
   */
  template <typename ITopologyView, typename ICoordsetView, typename IMatsetView>
  void iteration(const ITopologyView &topoView,
                 const ICoordsetView &coordsetView,
                 const IMatsetView &matsetView,

                 const axom::mir::views::Material &currentMat,

                 const conduit::Node &n_topo,
                 const conduit::Node &n_coordset,
                 const conduit::Node &n_fields,
                 const conduit::Node &n_matset,

                 const conduit::Node &n_options,

                 conduit::Node &n_newTopo,
                 conduit::Node &n_newCoordset,
                 conduit::Node &n_newFields,
                 conduit::Node &n_newMatset)
  {
    namespace bputils = axom::mir::utilities::blueprint;
    namespace bpmeshutils = conduit::blueprint::mesh::utils;
    std::cout << "------------------------ start of iteration "
                 "--------------------------------\n";

    const std::string colorField("__equiz__colors");

    // Make a node to zone relation.
    conduit::Node relation;
    {
      bputils::NodeToZoneRelationBuilder<ExecSpace> rb;
      rb.execute(n_topo, n_coordset, relation);
      std::cout << "\nnodeToZoneRelation = \n";
      printNode(relation);
      std::cout.flush();
    }

    // Make a shallow copy of the fields that we can modify.
    conduit::Node tmpFields;
    for(conduit::index_t i = 0; i < n_fields.number_of_children(); i++)
    {
      tmpFields[n_fields[i].name()].set_external(n_fields[i]);
    }

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

      // Make a matset intersector.
      using ConnectivityType = typename ITopologyView::ConnectivityType;
      using IntersectorType = MatsetIntersector<ConnectivityType, IMatsetView, RelationViewType>;
      IntersectorType intersector(matsetView, rv);

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
                        tmpFields,
                        options,
                        n_newTopo,
                        n_newCoordset,
                        n_newFields);
      }

      // Make a new matset.
      auto colorView =
        bputils::make_array_view<int>(n_newFields[colorField + "/values"]);
      auto origElemView =
        bputils::make_array_view<ConnectivityType>(n_newFields["originalElements/values"]);
      makeNewMatset<IMatsetView, ConnectivityType>(matsetView,
                                 currentMat,
                                 colorView,
                                 origElemView,
                                 n_matset,
                                 n_newMatset);

      // Save clip results.
      conduit::Node mesh;
      mesh[n_newTopo.path()].set_external(n_newTopo);
      mesh[n_newCoordset.path()].set_external(n_newCoordset);
      mesh[n_newFields.path()].set_external(n_newFields);
      mesh[n_newMatset.path()].set_external(n_newMatset);
#if 1
      printNode(mesh);
      std::cout << "------------------------ end of iteration "
                   "--------------------------------\n";
      std::cout.flush();

      std::stringstream ss;
      ss << "debug_equiz_" << currentMat.number;
      conduit::relay::io::blueprint::save_mesh(mesh, ss.str(), "hdf5");
#endif

      n_newFields.remove("originalElements");

    });
  }

  /**
   * \brief Make a new matset for the output dataset.
   *
   * \param curentMatsetView The view for the current matset.
   * \param currentMat The current material information.
   * \param colorView A view for the color field, which says how the dataset was clipped.
   * \param originalElementsView A view for the field that contains the original element number.
   * \param n_matset The input matset (pre-MIR).
   * \param n_newMatset The output matset.
   */
  template <typename IMatsetView, typename ConnectivityType>
  void makeNewMatset(  // by value
    const IMatsetView currentMatsetView,
    const axom::mir::views::Material currentMat,
    const axom::ArrayView<int> colorView,
    const axom::ArrayView<ConnectivityType> originalElementsView,
    // by reference
    const conduit::Node &n_matset,
    conduit::Node &n_newMatset)
  {
    namespace bputils = axom::mir::utilities::blueprint;
    using reduce_policy =
      typename axom::execution_space<ExecSpace>::reduce_policy;
    using IntType = typename IMatsetView::IndexType;
    using FloatType = typename IMatsetView::FloatType;
    const int intTypeID = bputils::cpp2conduit<IntType>::id;
    const int floatTypeID = bputils::cpp2conduit<FloatType>::id;

    // Copy the material map to the new matset.
    if(n_matset.has_child("material_map"))
    {
      bputils::copy<ExecSpace>(n_newMatset["material_map"],
                               n_matset.fetch_existing("material_map"));
    }
    if(n_matset.has_child("topology"))
    {
      n_newMatset["topology"].set(n_matset.fetch_existing("topology"));
    }

    // Get the ID of a Conduit allocator that will allocate through Axom with device allocator allocatorID.
    utilities::blueprint::ConduitAllocateThroughAxom<ExecSpace> c2a;

    // Make sizes, offsets for the new matset.
    const auto nzones = colorView.size();
    conduit::Node &n_sizes = n_newMatset["sizes"];
    conduit::Node &n_offsets = n_newMatset["offsets"];
    n_sizes.set_allocator(c2a.getConduitAllocatorID());
    n_offsets.set_allocator(c2a.getConduitAllocatorID());
    n_sizes.set(conduit::DataType(intTypeID, nzones));
    n_offsets.set(conduit::DataType(intTypeID, nzones));

    // Count how many slots are needed for the new material. Make sizes.
    RAJA::ReduceSum<reduce_policy, axom::IndexType> sum(0);
    auto sizesView = bputils::make_array_view<IntType>(n_sizes);
    std::cout << "makeNewMatset: nzones=" << nzones << ", currentMatNumber=" << currentMat.number << std::endl;
    const int currentMatNumber = currentMat.number;
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(auto zoneIndex) {
        axom::IndexType nmatsThisZone = 0;

        const auto origIndex = originalElementsView[zoneIndex];
        nmatsThisZone = currentMatsetView.numberOfMaterials(origIndex);
        //int dominantMatNumber = currentMatsetView.dominantMaterial(origIndex);
        const auto dominantMatNumber = m_matsetView.backgroundMaterial(origIndex, currentMatNumber);

        if(dominantMatNumber != currentMatNumber && nmatsThisZone > 1)
        {
          int matToRemove = (colorView[zoneIndex] == 0) ? currentMatNumber : dominantMatNumber;
          if(currentMatsetView.zoneContainsMaterial(origIndex, matToRemove))
          {
            nmatsThisZone--;
          }
        }
#if !defined(AXOM_DEVICE_CODE)
        if(zoneMatters(origIndex))
        {
          std::cout << "\tzone " << zoneIndex
                    << ": origElem=" << origIndex
                    << ", color=" << colorView[zoneIndex]
                    << ", dominant=" << dominantMatNumber
                    << ", nmats=" << nmatsThisZone << std::endl;
        }
#endif

        sizesView[zoneIndex] = nmatsThisZone;
        sum += nmatsThisZone;
      });

    // Make offsets.
    auto offsetsView = bputils::make_array_view<IntType>(n_offsets);
    axom::exclusive_scan<ExecSpace>(sizesView, offsetsView);

    // Make new matset data.
    const auto arraySize = sum.get();
    std::cout << "arraySize=" << arraySize << std::endl;
    conduit::Node &n_material_ids = n_newMatset["material_ids"];
    conduit::Node &n_volume_fractions = n_newMatset["volume_fractions"];
    conduit::Node &n_indices = n_newMatset["indices"];
    n_material_ids.set_allocator(c2a.getConduitAllocatorID());
    n_volume_fractions.set_allocator(c2a.getConduitAllocatorID());
    n_indices.set_allocator(c2a.getConduitAllocatorID());
    n_material_ids.set(conduit::DataType(intTypeID, arraySize));
    n_volume_fractions.set(conduit::DataType(floatTypeID, arraySize));
    n_indices.set(conduit::DataType(intTypeID, arraySize));

    auto matidsView = bputils::make_array_view<IntType>(n_material_ids);
    auto vfView = bputils::make_array_view<FloatType>(n_volume_fractions);
    auto indicesView = bputils::make_array_view<IntType>(n_indices);

    // Fill in the indices used.
    axom::for_all<ExecSpace>(
      arraySize,
      AXOM_LAMBDA(auto index) { indicesView[index] = index; });

    // Fill in the material data.
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(auto zoneIndex) {
        auto offset = offsetsView[zoneIndex];

        const auto origIndex = originalElementsView[zoneIndex];
        //int dominantMatNumber = currentMatsetView.dominantMaterial(origIndex);
        const auto dominantMatNumber = m_matsetView.backgroundMaterial(origIndex, currentMatNumber);

        // Get the material ids and volume fractions from the original material's zone.
        typename IMatsetView::IDList ids {};
        typename IMatsetView::VFList vfs {};
        currentMatsetView.zoneMaterials(origIndex, ids, vfs);

        if(dominantMatNumber != currentMatNumber && ids.size() > 1)
        {
          int matToRemove = (colorView[zoneIndex] == 0) ? currentMatNumber : dominantMatNumber;
          // Copy it through without colorToRemove.
          FloatType vfTotal {};
          for(axom::IndexType i = 0; i < ids.size(); i++)
          {
            vfTotal += (ids[i] != matToRemove) ? vfs[i] : 0;
          }
          for(axom::IndexType i = 0; i < ids.size(); i++)
          {
            if(ids[i] != matToRemove)
            {
              matidsView[offset] = ids[i];
              vfView[offset] = vfs[i] / vfTotal;

#if !defined(AXOM_DEVICE_CODE)
              if(zoneMatters(origIndex))
              {
std::cout << "\tzone " << zoneIndex
          << ": origElem=" << origIndex
          << ", color=" << colorView[zoneIndex]
          << ", current=" << currentMatNumber
          << ", dominant=" << dominantMatNumber
          << ", ids.size=" << ids.size()
          << ", offset=" << offset
          << ", id=" << ids[i]
          << ", vf=" << vfView[offset]
          << ", matToRemove=" << matToRemove
          << std::endl;
              }
#endif
              offset++;
            }
          }
        }
        else
        {
          // Copy it through.
          for(axom::IndexType i = 0; i < ids.size(); i++)
          {
            matidsView[offset] = ids[i];
            vfView[offset] = vfs[i];

#if !defined(AXOM_DEVICE_CODE)
              if(zoneMatters(origIndex))
              {
std::cout << "\tzone " << zoneIndex
          << ": origElem=" << origIndex
          << ", color=" << colorView[zoneIndex]
          << ", dominant=" << dominantMatNumber
          << ", offset=" << offset
          << ", id=" << ids[i]
          << ", vf=" << vfView[offset]
          << " (pass through)"
          << std::endl;
              }
#endif
            offset++;
          }
        }

/****
        if(colorView[zoneIndex] == COLOR0)
        {
#if !defined(AXOM_DEVICE_CODE)
          std::cout << "\tzone " << zoneIndex << ": offset=" << offset
                    << ", id=" << currentMatNumber << ", vf=1" << std::endl;
#endif
          matidsView[offset] = currentMatNumber;
          vfView[offset] = 1;
        }
        else
        {
          // Get the material ids and volume fractions from the original material's zone.
          const auto origIndex = originalElementsView[zoneIndex];
          typename IMatsetView::IDList ids {};
          typename IMatsetView::VFList vfs {};
          currentMatsetView.zoneMaterials(origIndex, ids, vfs);

          if(ids.size() > 1)
          {
            // The original zone was mixed. Filter out currentMatNumber

            FloatType vfTotal {};
            for(axom::IndexType i = 0; i < ids.size(); i++)
            {
              vfTotal += (ids[i] != currentMatNumber) ? vfs[i] : 0;
            }

            // Fill in the new materials.
            for(axom::IndexType i = 0; i < ids.size(); i++)
            {
              if(ids[i] != currentMatNumber)
              {
                matidsView[offset] = ids[i];
                vfView[offset] = vfs[i] / vfTotal;

#if !defined(AXOM_DEVICE_CODE)
std::cout << "\tzone " << zoneIndex << ": origElem=" << origIndex << ", offset=" << offset << ", id=" << ids[i] << ", vf=" << vfView[offset] << std::endl;
#endif
                offset++;
              }
            }
          }
          else
          {
            matidsView[offset] = ids[0];
            vfView[offset] = 1;
          }
        }
****/
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
