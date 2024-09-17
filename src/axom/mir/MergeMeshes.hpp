// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_MERGE_MESHES_HPP_
#define AXOM_MIR_MERGE_MESHES_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/mir.hpp"
#include "axom/slic.hpp"

#include <conduit/conduit.hpp>

#include <string>

#if defined(AXOM_USE_RAJA)
  #include <RAJA/RAJA.hpp>
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
 * \brief A mesh input containing a Blueprint mesh and some mapping array views.
 */
struct MeshInput
{
  conduit::Node                   *m_input{nullptr}; //!< Pointer to Blueprint mesh.
  axom::ArrayView<axom::IndexType> m_nodeMapView{};  //!< Map for mesh nodeIds to nodeIds in final mesh.
  axom::ArrayView<axom::IndexType> m_nodeSliceView{};//!< Node ids to be extracted and added to final mesh.
};

/**
 * \brief Merge multiple unstructured Blueprint meshes through MeshInput.
 *
 * \note The input meshes must currently contain a single coordset/topology/matset.
 */
template <typename ExecSpace>
class MergeMeshes
{
public:
  /**
   * \brief Merge the input Blueprint meshes into a single Blueprint mesh.
   *
   * \param inputs A vector of inputs to be merged.
   * \param options A Node containing algorithm options.
   * \param[out] output The node that will contain the merged mesh.
   */
  void execute(const std::vector<MeshInput> &inputs, const conduit::Node &options, conduit::Node &output) const
  {
    AXOM_ANNOTATE_SCOPE("MergeMeshes");
    bool ok = validInputs(inputs, options);
    if(!ok)
    {
      SLIC_ASSERT_MSG(ok, "Unsupported inputs were provided.");
      return;
    }

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
  /**
   * \brief Check that the mesh inputs are valid and meet constraints. There must
   *        be 1 coordset/topology/matset. The coordset must be explicit and the
   *        topology must be unstructured and for now, non-polyhedral.
   *
   * \param inputs The mesh inputs.
   *
   * \return True if the inputs appear to be valid; False otherwise.
   */
  bool validInputs(const std::vector<MeshInput> &inputs, const conduit::Node &options) const
  {
    std::string topoName;
    if(options.has_child("topology"))
      topoName = options.fetch_existing("topology").as_string();

    for(size_t i = 0; i < inputs.size(); i++)
    {
      if(inputs[i].m_input == nullptr)
        return false;

      // If we did not specify which topology, make sure that there is only 1.
      const char *keys[] = {"coordsets", "topologies", "matsets"};
      if(topoName.empty())
      {
        for(int k = 0; k < 3; k++)
        {
          if(inputs[i].m_input->has_path(keys[k]))
          {
            const conduit::Node &n = inputs[i].m_input->fetch_existing(keys[k]);
            if(n.number_of_children() > 1)
              return false;
          }
        }
      }

      const conduit::Node &n_coordsets = inputs[i].m_input->fetch_existing("coordsets");
      const conduit::Node &n_coordset = n_coordsets[0];
      if(n_coordset["type"].as_string() != "explicit")
      {
         return false;
      }

      const conduit::Node &n_topologies = inputs[i].m_input->fetch_existing("topologies");
      const conduit::Node *n_topo = nullptr;
      if(topoName.empty())
        n_topo = n_topologies.child_ptr(0);
      else
        n_topo = n_topologies.fetch_ptr(topoName);
      if(n_topo->operator[]("type").as_string() != "unstructured")
      {
         return false;
      }
      // For now
      if(n_topo->operator[]("elements/shape").as_string() == "polyhedral")
      {
         return false;
      }

      // Require no nodeMap/nodeSlice or that they both be present.
      if(!((inputs[i].m_nodeMapView.size() == 0 && inputs[i].m_nodeSliceView.size() == 0) ||
           (inputs[i].m_nodeMapView.size() > 0 && inputs[i].m_nodeSliceView.size() > 0)))
      {
         return false;
      }
    }
    return true;
  }

  /**
   * \brief Merge a single input (copy it to the output).
   *
   * \param inputs A vector of inputs to be merged.
   * \param[out] output The node that will contain the merged mesh. 
   */
  void singleInput(const std::vector<MeshInput> &inputs, conduit::Node &output) const
  {
    namespace bputils = axom::mir::utilities::blueprint;
    bputils::copy<ExecSpace>(output, *(inputs[0].m_input));
  }

  /**
   * \brief Merge a multiple inputs.
   *
   * \param inputs A vector of inputs to be merged.
   * \param[out] output The node that will contain the merged mesh. 
   */
  void mergeInputs(const std::vector<MeshInput> &inputs, conduit::Node &output) const
  {
    mergeCoordset(inputs, output);
    mergeTopology(inputs, output);
    //mergeFields(inputs, output);
    //mergeMatset(inputs, output);
  }

  /**
   * \brief Merge multiple coordsets into a single coordset. No node merging takes place.
   *
   * \param inputs A vector of inputs to be merged.
   * \param[out] output The node that will contain the merged coordset. 
   */
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

  axom::IndexType countNodes(const std::vector<MeshInput> &inputs, size_t index) const
  {
    SLIC_ASSERT(index < inputs.size());

    const conduit::Node &coordsets = inputs[index].m_input->fetch_existing("coordsets");
    const conduit::Node &coordset = coordsets[0];
    const auto type = coordset.fetch_existing("type").as_string();

    axom::IndexType nnodes = 0;
    if(inputs[index].m_nodeSliceView.size() > 0)
      nnodes = inputs[index].m_nodeSliceView.size();
    else
      nnodes = conduit::blueprint::mesh::utils::coordset::length(coordset);

    return nnodes;
  }

  axom::IndexType countNodes(const std::vector<MeshInput> &inputs) const
  {
    axom::IndexType nodeTotal = 0;
    for(size_t i = 0; i < inputs.size(); i++)
    {
      nodeTotal += countNodes(inputs, i);
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
          const std::string sname = n_shape_map[s].name();
          const auto id = axom::mir::views::shapeNameToID(sname);
          SLIC_ASSERT(id == n_shape_map[s].to_int());
          shape_map[sname] = id;
        }
      }
      else
      {
        shape_map[shape] = axom::mir::views::shapeNameToID(shape);
      }
    }

    conduit::Node *n_newTopoPtr = nullptr;
    axom::IndexType connOffset = 0, sizesOffset = 0, shapesOffset = 0, coordOffset = 0;
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
          n_newTopo["elements/shape"] = "mixed";

          // Build a new shape map in the new topology.
          conduit::Node &n_shape_map = n_newTopo["elements/shape_map"];
          for(auto it = shape_map.begin(); it != shape_map.end(); it++)
            n_shape_map[it->first] = it->second;

          conduit::Node &n_newShapes = n_newTopo["elements/shapes"];
          n_newShapes.set_allocator(c2a.getConduitAllocatorID());
          n_newShapes.set(conduit::DataType(n_srcConn.dtype().id(), totalZones));
        }
        else
        {
          n_newTopo["elements/shape"] = shape_map.begin()->first;
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
          // The supplied nodeMap is assumed to be a mapping from the current
          // node connectivity to the merged node connectivity.
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
          // Copy all zones from the input. Map the nodes to the new values.
          axom::for_all<ExecSpace>(srcConnView.size(), AXOM_LAMBDA(auto index)
          {
            connView[connOffset + index] = coordOffset + srcConnView[index];
          });
        }
        connOffset += srcConnView.size();
        coordOffset += countNodes(inputs, static_cast<size_t>(i));
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
          conduit::Node &n_newShapes = n_newTopoPtr->fetch_existing("elements/shapes");
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
        // Fill in shape information. There is no source shape information. Use
        // sizes to get the number of zones.
        const conduit::Node &n_srcSizes = n_srcTopo.fetch_existing("elements/sizes");
        axom::IndexType nz = n_srcSizes.dtype().number_of_elements();
        conduit::Node &n_newShapes = n_newTopoPtr->fetch_existing("elements/shapes");
        axom::mir::views::IndexNode_to_ArrayView(n_newShapes, [&](auto shapesView)
        {
          const int shapeId = axom::mir::views::shapeNameToID(srcShape);
          axom::for_all<ExecSpace>(nz, AXOM_LAMBDA(auto index)
          {
            shapesView[shapesOffset + index] = shapeId;
          });
          shapesOffset += nz;
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

}  // end namespace blueprint
}  // end namespace utilities
}  // end namespace mir
}  // end namespace axom

#endif
