// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_EQUIZ_ALGORITHM_HPP_
#define AXOM_MIR_EQUIZ_ALGORITHM_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/mir.hpp"
#include "axom/slic.hpp"

// Include these directly for now.
#include "axom/mir/views/MaterialView.hpp"
#include "axom/mir/MIRAlgorithm.hpp"
#include "axom/mir/RecenterField.hpp"
#include "axom/mir/NodeToZoneRelationBuilder.hpp"
#include "axom/mir/ZoneListBuilder.hpp"
#include "axom/mir/ExtractZones.hpp"

#include <conduit/conduit.hpp>

#include <algorithm>
#include <string>

#if defined(AXOM_USE_RAJA)
  #include <RAJA/RAJA.hpp>
#endif

// Uncomment to save inputs and outputs.
#define AXOM_DEBUG_EQUIZ

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

      axom::IndexType clipcase = 0;
      const auto n = nodeIdsView.size();
      for(IndexType i = 0; i < n; i++)
      {
        const auto nid = nodeIdsView[i];
#if defined(AXOM_DEVICE_CODE)
        assert(nid >= 0 && nid < m_matvfViews[0].size());
#else
        SLIC_ASSERT_MSG(nid >= 0 && nid < m_matvfViews[0].size(),
                        axom::fmt::format("Node id {} is not in range [0, {}).",
                                          nid,
                                          m_matvfViews[0].size()));
#endif
        // clang-format off
        MaterialVF vf1 = (backgroundIndex != INVALID_INDEX) ? m_matvfViews[backgroundIndex][nid] : NULL_MATERIAL_VF;
        MaterialVF vf2 = (m_currentMaterialIndex != INVALID_INDEX) ? m_matvfViews[m_currentMaterialIndex][nid] : 0;
        // clang-format on

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

#if defined(AXOM_DEVICE_CODE)
      assert(id0 >= 0 && id0 < m_matvfViews[0].size());
      assert(id1 >= 0 && id1 < m_matvfViews[0].size());
#else
      SLIC_ASSERT_MSG(id0 >= 0 && id0 < m_matvfViews[0].size(),
                      axom::fmt::format("Node id {} is not in range [0, {}).",
                                        id0,
                                        m_matvfViews[0].size()));
      SLIC_ASSERT_MSG(id1 >= 0 && id1 < m_matvfViews[0].size(),
                      axom::fmt::format("Node id {} is not in range [0, {}).",
                                        id1,
                                        m_matvfViews[0].size()));
#endif

      // Get the volume fractions for mat1, mat2 at the edge endpoints id0, id1.
      MaterialVF vf1[2], vf2[2];
      // clang-format off
      vf1[0] = (backgroundIndex != INVALID_INDEX) ? m_matvfViews[backgroundIndex][id0] : NULL_MATERIAL_VF;
      vf1[1] = (backgroundIndex != INVALID_INDEX) ? m_matvfViews[backgroundIndex][id1] : NULL_MATERIAL_VF;
      vf2[0] = (m_currentMaterialIndex != INVALID_INDEX) ? m_matvfViews[m_currentMaterialIndex][id0] : 0;
      vf2[1] = (m_currentMaterialIndex != INVALID_INDEX) ? m_matvfViews[m_currentMaterialIndex][id1] : 0;
      // clang-format on

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
    inline int matNumberToIndex(int matNumber) const
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

    void setCurrentMaterial(int matNumber, int matNumberIndex)
    {
      m_currentMaterial = matNumber;
      m_currentMaterialIndex = matNumberIndex;
    }

    axom::StaticArray<MaterialVFView, MAXMATERIALS>
      m_matvfViews {};  //!< Array of volume fraction views
    axom::ArrayView<int> m_matNumbersView {};  //!< Sorted array of material numbers.
    axom::ArrayView<int> m_matIndicesView {};  //!< Array of indices into m_matvfViews for the material numbers.
    axom::ArrayView<int> m_zoneMatNumberView {};  //!< Contains the current material number that owns each zone.
    int m_currentMaterial {};  //!< The current material.
    int m_currentMaterialIndex {};  //!< The current material's index in the m_matvfViews.
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

  void setCurrentMaterial(int matNumber, int matNumberIndex)
  {
    m_view.setCurrentMaterial(matNumber, matNumberIndex);
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

/**
 * \brief Merge multiple unstructured Blueprint meshes
 *
 * \note The input meshes must currently contain a single coordset/topology/matset.
 */
template <typename ExecSpace>
class MergeMeshes
{
public:
  struct MeshInput
  {
    conduit::Node                   *m_input{nullptr};
    axom::ArrayView<axom::IndexType> m_nodeMapView{};
    axom::ArrayView<axom::IndexType> m_nodeSliceView{};
  };

  void execute(const std::vector<MeshInput> &inputs, conduit::Node &output) const
  {
    AXOM_ANNOTATE_SCOPE("MergeMeshes");
    SLIC_ASSERT(validInputs(inputs));

    if(inputs.size() == 1)
    {
      singleInput(inputs, output);
    }
    else if(inputs.size() > 1)
    {
      mergeInputs(inputs, output);
    }
  }
private:
  bool validInputs(const std::vector<MeshInput> &inputs) const
  {
    for(size_t i = 0; i < inputs.size(); i++)
    {
      if(inputs[i].m_input == nullptr)
        return false;

      const char *keys[] = {"coordsets", "topologies", "matsets"};
      for(int k = 0; k < 3; k++)
      {
        if(inputs[i].m_input->has_path(keys[k]))
        {
          const conduit::Node &n = inputs[i].m_input->fetch_existing(keys[k]);
          if(n.number_of_children() > 1)
            return false;
        }
      }

      inputs[i].m_input->fetch_existing("topologies");
    }
    return true;
  }

  void singleInput(const std::vector<MeshInput> &inputs, conduit::Node &output) const
  {
    namespace bputils = axom::mir::utilities::blueprint;
    bputils::copy<ExecSpace>(output, *(inputs[0].m_input));
  }

  void mergeInputs(const std::vector<MeshInput> &inputs, conduit::Node &output) const
  {
    mergeCoordset(inputs, output);
    mergeTopology(inputs, output);
    //mergeFields(inputs, output);
    //mergeMatset(inputs, output);
  }

  void mergeCoordset(const std::vector<MeshInput> &inputs, conduit::Node &output) const
  {
    namespace bputils = axom::mir::utilities::blueprint;
    AXOM_ANNOTATE_SCOPE("mergeCoordset");
    const axom::IndexType totalNodes = countNodes(inputs);
    conduit::Node &n_newCoordsets = output["coordsets"];
    conduit::Node *n_newValuesPtr = nullptr;
    int nComps = 2;

    axom::IndexType offsets[3] = {0, 0, 0};
    const axom::IndexType n = static_cast<axom::IndexType>(inputs.size());
    for(axom::IndexType i = 0; i < n; i++)
    {
      const conduit::Node &coordsets = inputs[i].m_input->fetch_existing("coordsets");
      const conduit::Node &n_srcCoordset = coordsets[0];
      const conduit::Node &n_srcValues = n_srcCoordset.fetch_existing("values");

      const auto type = n_srcCoordset.fetch_existing("type").as_string();
      SLIC_ASSERT(type == "explicit");

      // Make all of the components the first time.
      if(i == 0)
      {
        bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;

        conduit::Node &n_newCoordset = n_newCoordsets[n_srcCoordset.name()];
        n_newCoordset["type"] = "explicit";
        conduit::Node &n_newValues = n_newCoordset["values"];
        n_newValuesPtr = n_newCoordset.fetch_ptr("values");

        nComps = n_srcValues.number_of_children();
        for(int c = 0; c < nComps; c++)
        {
          const conduit::Node &n_srcComp = n_srcValues[c];
          conduit::Node &n_comp = n_newValues[n_srcComp.name()];
          n_comp.set_allocator(c2a.getConduitAllocatorID());
          n_comp.set(conduit::DataType(n_srcComp.dtype().id(), totalNodes));
        }
      }

      // Copy this input's coordinates into the new coordset.
      axom::mir::views::FloatNode_to_ArrayView(n_srcValues[0], [&](auto comp0)
      {
        using FloatType = typename decltype(comp0)::value_type;
        for(int c = 0; c < nComps; c++)
        {
          axom::IndexType size = 0, offset = offsets[c];

          const conduit::Node &n_srcComp = n_srcValues[c];
          conduit::Node &n_comp = n_newValuesPtr->child(c);
          auto srcCompView = bputils::make_array_view<FloatType>(n_srcComp);
          auto compView = bputils::make_array_view<FloatType>(n_comp);

          const auto nodeSliceView = inputs[i].m_nodeSliceView;
          if(nodeSliceView.size() > 0)
          {
            // Pull out specific nodes from the input.
            axom::for_all<ExecSpace>(nodeSliceView.size(), AXOM_LAMBDA(auto index)
            {
              const auto sliceIndex = nodeSliceView[index];
              compView[offset + index] = srcCompView[sliceIndex];
            });
            size = nodeSliceView.size();
          }
          else
          {
            // Pull out all nodes from the input.
            axom::for_all<ExecSpace>(srcCompView.size(), AXOM_LAMBDA(auto index)
            {
              compView[offset + index] = srcCompView[index];
            });
            size = srcCompView.size();
          }

          offsets[c] += size;
        }
      });
    }
  }

  axom::IndexType countNodes(const std::vector<MeshInput> &inputs) const
  {
    axom::IndexType nodeTotal = 0;
    axom::IndexType n = static_cast<axom::IndexType>(inputs.size());
    for(axom::IndexType i = 0; i < n; i++)
    {
      const conduit::Node &coordsets = inputs[i].m_input->fetch_existing("coordsets");
      const conduit::Node &coordset = coordsets[0];
      const auto type = coordset.fetch_existing("type").as_string();

      axom::IndexType nnodes = 0;
      if(inputs[i].m_nodeSliceView.size() > 0)
        nnodes = inputs[i].m_nodeSliceView.size();
      else
        nnodes = conduit::blueprint::mesh::utils::coordset::length(coordset);

      nodeTotal += nnodes;
    }
    return nodeTotal;
  }

  void countZones(const std::vector<MeshInput> &inputs, axom::IndexType &totalConnLength, axom::IndexType &totalZones) const
  {
    totalConnLength = 0;
    totalZones = 0;
    axom::IndexType n = static_cast<axom::IndexType>(inputs.size());
    for(axom::IndexType i = 0; i < n; i++)
    {
      const conduit::Node &n_topologies = inputs[i].m_input->fetch_existing("topologies");
      const conduit::Node &n_topo = n_topologies[0];
      const auto type = n_topo.fetch_existing("type").as_string();
      SLIC_ASSERT(type == "unstructured");

      const conduit::Node &n_conn = n_topo.fetch_existing("elements/connectivity");
      totalConnLength += n_conn.dtype().number_of_elements();

      const conduit::Node &n_size = n_topo.fetch_existing("elements/sizes");
      totalZones += n_size.dtype().number_of_elements();
    }
  }

  void mergeTopology(const std::vector<MeshInput> &inputs, conduit::Node &output) const
  {
    namespace bputils = axom::mir::utilities::blueprint;
    AXOM_ANNOTATE_SCOPE("mergeTopology");
    axom::IndexType totalConnLen = 0, totalZones = 0;
    countZones(inputs, totalConnLen, totalZones);  
    conduit::Node &n_newTopologies = output["topologies"];

    // Check whether there are mixed shapes.
    std::map<std::string, int> shape_map;
    const axom::IndexType n = static_cast<axom::IndexType>(inputs.size());
    for(axom::IndexType i = 0; i < n; i++)
    {
      const conduit::Node &n_topologies = inputs[i].m_input->fetch_existing("topologies");
      const conduit::Node &n_srcTopo = n_topologies[0];
      const auto type = n_srcTopo.fetch_existing("type").as_string();
      const auto shape = n_srcTopo.fetch_existing("elements/shape").as_string();
      SLIC_ASSERT(type == "unstructured");
      if(shape == "mixed")
      {
        const conduit::Node &n_shape_map = n_srcTopo.fetch_existing("elements/shape_map");
        for(int s = 0; s < n_shape_map.number_of_children(); s++)
        {
          const std::string sname = n_shape_map[i].name();
          const auto id = axom::mir::views::shapeNameToID(sname);
          SLIC_ASSERT(id == n_shape_map[i].to_int());
          shape_map[sname] = id;
        }
      }
      else
      {
        shape_map[shape] = axom::mir::views::shapeNameToID(shape);
      }
    }

    conduit::Node *n_newTopoPtr = nullptr;
    axom::IndexType connOffset = 0, sizesOffset = 0, shapesOffset = 0;
    for(axom::IndexType i = 0; i < n; i++)
    {
      const conduit::Node &n_topologies = inputs[i].m_input->fetch_existing("topologies");
      const conduit::Node &n_srcTopo = n_topologies[0];

      const std::string srcShape = n_srcTopo.fetch_existing("elements/shape").as_string();
      const conduit::Node &n_srcConn = n_srcTopo.fetch_existing("elements/connectivity");
      const conduit::Node &n_srcSizes = n_srcTopo.fetch_existing("elements/sizes");

      // Make all of the elements the first time.
      if(i == 0)
      {
        bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;

        conduit::Node &n_newTopo = n_newTopologies[n_srcTopo.name()];
        n_newTopoPtr = n_newTopologies.fetch_ptr(n_srcTopo.name());
        n_newTopo["type"] = "unstructured";
        n_newTopo["coordset"] = n_srcTopo["coordset"].as_string();

        conduit::Node &n_newConn = n_newTopo["elements/connectivity"];
        n_newConn.set_allocator(c2a.getConduitAllocatorID());
        n_newConn.set(conduit::DataType(n_srcConn.dtype().id(), totalConnLen));

        conduit::Node &n_newSizes = n_newTopo["elements/sizes"];
        n_newSizes.set_allocator(c2a.getConduitAllocatorID());
        n_newSizes.set(conduit::DataType(n_srcSizes.dtype().id(), totalZones));

        conduit::Node &n_newOffsets = n_newTopo["elements/offsets"];
        n_newOffsets.set_allocator(c2a.getConduitAllocatorID());
        n_newOffsets.set(conduit::DataType(n_srcConn.dtype().id(), totalZones));

        if(shape_map.size() > 1)
        {
          // Build a new shape map in the new topology.
          conduit::Node &n_shape_map = n_newTopo["shape_map"];
          for(auto it = shape_map.begin(); it != shape_map.end(); it++)
            n_shape_map[it->first] = it->second;

          conduit::Node &n_newShapes = n_newTopo["elements/shapes"];
          n_newShapes.set_allocator(c2a.getConduitAllocatorID());
          n_newShapes.set(conduit::DataType(n_srcConn.dtype().id(), totalZones));
        }
      }

      // Copy this input's connectivity into the new topology.
      axom::mir::views::IndexNode_to_ArrayView(n_srcConn, [&](auto srcConnView)
      {
        using ConnType = typename decltype(srcConnView)::value_type;
        conduit::Node &n_newConn = n_newTopoPtr->fetch_existing("elements/connectivity");
        auto connView = bputils::make_array_view<ConnType>(n_newConn);

        if(inputs[i].m_nodeMapView.size() > 0)
        {
          // Copy all zones from the input but map the nodes to new values.
          const auto nodeMapView = inputs[i].m_nodeMapView;
          axom::for_all<ExecSpace>(srcConnView.size(), AXOM_LAMBDA(auto index)
          {
            const auto nodeId = srcConnView[index];
            const auto newNodeId = nodeMapView[nodeId];
            connView[connOffset + index] = newNodeId;
          });
        }
        else
        {
          // Copy all zones from the input.
          axom::for_all<ExecSpace>(srcConnView.size(), AXOM_LAMBDA(auto index)
          {
            connView[connOffset + index] = srcConnView[index];
          });
        }
        connOffset += srcConnView.size();
      });

      // Copy this input's sizes into the new topology.
      axom::mir::views::IndexNode_to_ArrayView(n_srcSizes, [&](auto srcSizesView)
      {
        using ConnType = typename decltype(srcSizesView)::value_type;
        conduit::Node &n_newSizes = n_newTopoPtr->fetch_existing("elements/sizes");
        auto sizesView = bputils::make_array_view<ConnType>(n_newSizes);

        // Copy all sizes from the input.
        axom::for_all<ExecSpace>(srcSizesView.size(), AXOM_LAMBDA(auto index)
        {
          sizesView[sizesOffset + index] = srcSizesView[index];
        });

        sizesOffset += srcSizesView.size();
      });

      // Copy shape information if it exists.
      if(n_srcTopo.has_path("elements/shapes"))
      {
        const conduit::Node &n_srcShapes = n_srcTopo.fetch_existing("elements/shapes");

        axom::mir::views::IndexNode_to_ArrayView(n_srcShapes, [&](auto srcShapesView)
        {
          using ConnType = typename decltype(srcShapesView)::value_type;
          conduit::Node &n_newShapes = n_newTopoPtr->fetch_existing("elements/sizes");
          auto shapesView = bputils::make_array_view<ConnType>(n_newShapes);

          // Copy all sizes from the input.
          axom::for_all<ExecSpace>(srcShapesView.size(), AXOM_LAMBDA(auto index)
          {
            shapesView[shapesOffset + index] = srcShapesView[index];
          });

          shapesOffset += srcShapesView.size();
        });
      }
      else
      {
        // Fill in shape information.
        conduit::Node &n_newShapes = n_newTopoPtr->fetch_existing("elements/sizes");
        axom::mir::views::IndexNode_to_ArrayView(n_newShapes, [&](auto shapesView)
        {
          const int shapeId = axom::mir::views::shapeNameToID(srcShape);
          axom::for_all<ExecSpace>(shapesView.size(), AXOM_LAMBDA(auto index)
          {
            shapesView[shapesOffset + index] = shapeId;
          });
          shapesOffset += shapesView.size();
        });
      }
    }

    // Make new offsets from the sizes.
    conduit::Node &n_newSizes = n_newTopoPtr->fetch_existing("elements/sizes");
    axom::mir::views::IndexNode_to_ArrayView(n_newSizes, [&](auto sizesView)
    {
      using ConnType = typename decltype(sizesView)::value_type;
      conduit::Node &n_newOffsets = n_newTopoPtr->fetch_existing("elements/offsets");
      auto offsetsView = bputils::make_array_view<ConnType>(n_newOffsets);
      axom::exclusive_scan<ExecSpace>(sizesView, offsetsView);
    });
  }
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

#if 1
    // Come up with lists of clean/mixed zones.
    axom::Array<axom::IndexType> cleanZones, mixedZones;
    bputils::ZoneListBuilder<ExecSpace, TopologyView, MatsetView> zlb(m_topologyView, m_matsetView);
    if(n_options.has_child("selectedZones"))
    {
      auto selectedZonesView = bputils::make_array_view<axom::IndexType>(n_options.fetch_existing("selectedZones"));
      zlb.execute(m_coordsetView.numberOfNodes(), selectedZonesView, cleanZones, mixedZones);
    }
    else
    {
      zlb.execute(m_coordsetView.numberOfNodes(), cleanZones, mixedZones);
    }
    SLIC_ASSERT((cleanZones.size() + mixedZones.size()) == m_topologyView.numberOfZones());
std::cout << "cleanZones.size: " << cleanZones.size() << std::endl;
std::cout << "mixedZones.size: " << mixedZones.size() << std::endl;

    if(cleanZones.size() > 0 && mixedZones.size() > 0)
    {
      // Bundle the inputs so we can make the clean mesh
      conduit::Node n_ezInput;
      n_ezInput[n_topo.path()].set_external(n_topo);
      n_ezInput[n_coordset.path()].set_external(n_coordset);
      n_ezInput[n_fields.path()].set_external(n_fields);
      n_ezInput[n_matset.path()].set_external(n_matset);

std::cout << "Making clean mesh" << std::endl;
      // Make the clean mesh.
      bputils::ExtractZonesAndMatset<ExecSpace, TopologyView, CoordsetView, MatsetView> ez(m_topologyView, m_coordsetView, m_matsetView);
      conduit::Node n_ezopts, n_cleanOutput;
      n_ezopts["topology"] = n_topo.name();
      ez.execute(cleanZones, n_ezInput, n_ezopts, n_cleanOutput);

      // Get the number of nodes in the clean mesh.
      const conduit::Node &n_cleanCoordset = n_cleanOutput[n_coordset.path()];
      axom::IndexType numCleanNodes = n_cleanCoordset["values/x"].dtype().number_of_elements();

#if 1
std::cout << "Saving clean mesh" << std::endl;
      bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;

      AXOM_ANNOTATE_BEGIN("saveClean");
      conduit::relay::io::blueprint::save_mesh(n_cleanOutput, "clean", "hdf5");
      AXOM_ANNOTATE_END("saveClean");

      // Add a field to the original nodes.
      AXOM_ANNOTATE_BEGIN("origNodes");
      // Gather the inputs into a single root but replace the fields with
      // a new node to which we can add additional fields.
      conduit::Node n_root;
      n_root[n_coordset.path()].set_external(n_coordset);
      n_root[n_topo.path()].set_external(n_topo);
      n_root[n_matset.path()].set_external(n_matset);
      conduit::Node &n_root_coordset = n_root[n_coordset.path()];
      conduit::Node &n_root_topo = n_root[n_topo.path()];
      conduit::Node &n_root_matset = n_root[n_matset.path()];
      // Link in original fields.
      conduit::Node &n_root_fields = n_root["fields"];
      for(conduit::index_t i = 0; i < n_fields.number_of_children(); i++)
      {
        n_root_fields[n_fields[i].name()].set_external(n_fields[i]);
      }
      // Add a new field.
      conduit::Node &n_orig_nodes = n_root_fields["__equiz_original_node"];
      n_orig_nodes["topology"] = n_topo.name();
      n_orig_nodes["association"] = "vertex";
      n_orig_nodes["values"].set_allocator(c2a.getConduitAllocatorID());
      n_orig_nodes["values"].set(conduit::DataType::float32(m_coordsetView.numberOfNodes()));
      auto origNodesView = bputils::make_array_view<float>(n_orig_nodes["values"]);
      axom::for_all<ExecSpace>(m_coordsetView.numberOfNodes(), AXOM_LAMBDA(auto index)
      {
        origNodesView[index] = static_cast<float>(index);
      });
//printNode(n_root);
      // If there are fields in the options, make sure the new field is handled too.
      if(n_options_copy.has_child("fields"))
        n_options_copy["fields/__equiz_original_node"] = "__equiz_original_node";
      AXOM_ANNOTATE_END("origNodes");

//std::cout << "Making mixed mesh" << std::endl;
//      // Make a mixed mesh.
//      conduit::Node n_mixedOutput;
//      n_ezopts["topology"] = n_topo.name();
//      ez.execute(mixedZones, n_ezInput, n_ezopts, n_mixedOutput);
//std::cout << "Saving mixed mesh" << std::endl;
//      conduit::relay::io::blueprint::save_mesh(n_mixedOutput, "mixed", "hdf5");

#endif

std::cout << "Process mixed zones..." << std::endl;

      // Process the mixed part of the mesh. We select just the mixed zones.
      n_options_copy["selectedZones"].set_external(mixedZones.data(), mixedZones.size());
      n_options_copy["newNodesField"] = "__equiz_new_nodes";
printNode(n_options_copy);
      processMixedZones(n_root_topo, n_root_coordset, n_root_fields, n_root_matset,
                        n_options_copy,
                        n_newTopo, n_newCoordset, n_newFields, n_newMatset);
#if 1
      // Look over the nodes in the output and find the new ones.
      SLIC_ASSERT(n_newFields.has_child("__equiz_original_node"));
      AXOM_ANNOTATE_BEGIN("identifyOriginalNodes");
      const conduit::Node &n_output_orig_nodes = n_newFields["__equiz_original_node/values"];
      auto numOutputNodes = n_output_orig_nodes.dtype().number_of_elements();
std::cout << "numOutputNodes: " << numOutputNodes << std::endl;
      auto outputOrigNodesView = bputils::make_array_view<float>(n_output_orig_nodes);

      const conduit::Node &n_new_nodes_values = n_newFields["__equiz_new_nodes/values"];
      auto newNodesView = bputils::make_array_view<float>(n_new_nodes_values);

      // Make mask/offset for which nodes we need to keep from the mixed output.
      const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
      axom::Array<int> mask(numOutputNodes, numOutputNodes, allocatorID);
      axom::Array<int> maskOffset(numOutputNodes, numOutputNodes, allocatorID);
      auto maskView = mask.view();
      auto maskOffsetsView = mask.view();
      using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;
      RAJA::ReduceSum<reduce_policy, int> mask_reduce(0);
      axom::for_all<ExecSpace>(numOutputNodes, AXOM_LAMBDA(auto index)
      {
        const int ival = (newNodesView[index] > 0.f) ? 1 : 0;
        maskView[index] = ival;
        mask_reduce += ival;
      });
      axom::exclusive_scan<ExecSpace>(maskView, maskOffsetsView);
      const auto numNewNodes = mask_reduce.get();

      // Make a list of indices that we need to slice out of the node arrays.
      axom::Array<axom::IndexType> nodeSlice(numNewNodes, numNewNodes, allocatorID);
      auto nodeSliceView = nodeSlice.view();
      axom::for_all<ExecSpace>(numOutputNodes, AXOM_LAMBDA(auto index)
      {
        if(maskView[index] > 0)
        {
          nodeSliceView[maskOffsetsView[index]] = index;
        }
      });

      // Make a node map for mapping mixed connectivity into combined node numbering.
      axom::Array<axom::IndexType> nodeMap(numOutputNodes, numOutputNodes, allocatorID);
      auto nodeMapView = nodeMap.view();
      axom::for_all<ExecSpace>(numOutputNodes, AXOM_LAMBDA(auto index)
      {
        if(maskView[index] == 0)
        {
          nodeMapView[index] = static_cast<axom::IndexType>(outputOrigNodesView[index]);
        }
        else
        {
          nodeMapView[index] = numCleanNodes + maskOffsetsView[index];
        }
      });

std::cout << "Nodes that are unique to mixed output: " << numNewNodes << std::endl;
      AXOM_ANNOTATE_END("identifyOriginalNodes");

      //--------------------------------------------------------------------------
      //
      // Save the MIR output with the mask.
      //
      //--------------------------------------------------------------------------
      conduit::Node n_tmpMask;
      n_tmpMask[n_topo.path()].set_external(n_newTopo);
      n_tmpMask[n_coordset.path()].set_external(n_newCoordset);
      n_tmpMask[n_fields.path()].set_external(n_newFields);
      n_tmpMask[n_matset.path()].set_external(n_newMatset);
//      printNode(n_tmpMask);
      conduit::relay::io::blueprint::save_mesh(n_tmpMask,
                                               "debug_equiz_mask",
                                               "hdf5");

//      n_newFields.remove("__equiz_original_node");
//      n_newFields.remove("__equiz_mask");

#if 1
      conduit::Node n_merged;
      {
        AXOM_ANNOTATE_SCOPE("merge");
        // Merge the clean and mixed.
        using MeshInput = typename detail::MergeMeshes<ExecSpace>::MeshInput;
        std::vector<MeshInput> inputs;
        MeshInput i0;
        i0.m_input = &n_cleanOutput;
        inputs.push_back(i0);

        MeshInput i1;
        i1.m_input = &n_tmpMask;
        i1.m_nodeMapView = nodeMapView;
        i1.m_nodeSliceView = nodeSliceView;
        inputs.push_back(i1);

        detail::MergeMeshes<ExecSpace> mm;
        mm.execute(inputs, n_merged);
      }
      printNode(n_merged);
      conduit::relay::io::blueprint::save_mesh(n_merged,
                                               "debug_equiz_merged",
                                               "hdf5");
#endif

    
#endif
      // Merge the 2 parts together.
    }
    else if(cleanZones.size() == 0 && mixedZones.size() > 0)
    {
std::cout << "All zones are mixed..." << std::endl;
       // Only mixed zones.
      processMixedZones(n_topo, n_coordset, n_fields, n_matset,
                        n_options_copy,
                        n_newTopo, n_newCoordset, n_newFields, n_newMatset);
    }
    else if(cleanZones.size() > 0 && mixedZones.size() == 0)
    {
std::cout << "All zones are clean..." << std::endl;
      // There were no mixed zones. We can copy the input to the output.

      // Add an originalZones array.
    }
#else
std::cout << "Normal handling" << std::endl;
    // Handle all zones.
    processMixedZones(n_topo, n_coordset, n_fields, n_matset,
                      n_options_copy,
                      n_newTopo, n_newCoordset, n_newFields, n_newMatset);
#endif
  }

  void processMixedZones(const conduit::Node &n_topo,
                         const conduit::Node &n_coordset,
                         const conduit::Node &n_fields,
                         const conduit::Node &n_matset,
                         conduit::Node &n_options,
                         conduit::Node &n_newTopo,
                         conduit::Node &n_newCoordset,
                         conduit::Node &n_newFields,
                         conduit::Node &n_newMatset) const
  {
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

                                              n_options,

                                              n_newTopo,
                                              n_newCoordset,
                                              n_newFields);

        // In later iterations, we do not want to pass selectedZones through
        // since they are only valid on the current input topology. Also, if they
        // were passed then the new topology only has those selected zones.
        if(n_options.has_child("selectedZones"))
        {
          n_options.remove("selectedZones");
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
        // clang-format off
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

                                                      n_options,

                                                      n_newTopo,
                                                      n_newCoordset,
                                                      n_newFields);
            });
        });
        // clang-format on
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
      AXOM_ANNOTATE_SCOPE("relation");
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
std::cout << "makeNodeCenteredVFs: nzones=" << nzones << std::endl;
std::cout << "makeNodeCenteredVFs: nnodes=" << nnodes << std::endl;
    {
      AXOM_ANNOTATE_SCOPE("zonal");
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
      }
    }

    {
      AXOM_ANNOTATE_SCOPE("recenter");
      for(const auto &mat : mixedMats)
      {
        const int matNumber = mat.number;
        const std::string zonalName = zonalFieldName(matNumber);
        conduit::Node &n_zonalField = n_fields[zonalName];

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
                         const axom::mir::views::MaterialInformation
                           &AXOM_UNUSED_PARAM(mixedMats)
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
                 conduit::Node &n_newFields) const
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
      // Get the current material's index in the number:index map.
      int currentMatIndex = 0;
      for(axom::IndexType i = 0; i < matNumber.size(); i++)
      {
        if(matNumber[i] == currentMat.number)
        {
          currentMatIndex = matIndex[i];
          break;
        }
      }

      // Store the number:index map into the intersector. The number:index map lets us
      // ask for the field index for a material number, allowing scattered material
      // numbers to be used in the matset.
      axom::copy(matNumberDevice.data(), matNumber.data(), sizeof(int) * nmats);
      axom::copy(matIndexDevice.data(), matIndex.data(), sizeof(int) * nmats);
      intersector.setMaterialNumbers(matNumberDevice.view());
      intersector.setMaterialIndices(matIndexDevice.view());
      intersector.setCurrentMaterial(currentMat.number, currentMatIndex);

      // Store the current zone material ids and current material number into the intersector.
      intersector.setZoneMaterialID(bputils::make_array_view<MaterialID>(
        n_fields.fetch_existing(zonalMaterialIDName() + "/values")));
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
    if(n_options.has_child("fields"))
    {
      // Pass along fields, if present.
      options["fields"].set_external(n_options.fetch_existing("fields"));
    }
    if(n_options.has_child("newNodesField"))
    {
      // Pass along newNodesField, if present.
      options["newNodesField"] = n_options.fetch_existing("newNodesField").as_string();
    }
    options["topology"] = n_options["topology"];

    //--------------------------------------------------------------------------
    //
    // Clip the topology using the material intersector.
    //
    //--------------------------------------------------------------------------
    {
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

      const auto colorView = bputils::make_array_view<int>(
        n_newFields.fetch_existing(colorField + "/values"));
      const auto nzonesNew = colorView.size();

      // Get zonalMaterialID field so we can make adjustments.
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
