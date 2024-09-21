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

#include "axom/mir/views/dispatch_material.hpp"

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
  conduit::Node *m_input {nullptr};  //!< Pointer to Blueprint mesh.
  axom::ArrayView<axom::IndexType>
    m_nodeMapView {};  //!< Map for mesh nodeIds to nodeIds in final mesh.
  axom::ArrayView<axom::IndexType>
    m_nodeSliceView {};  //!< Node ids to be extracted and added to final mesh.
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
  void execute(const std::vector<MeshInput> &inputs,
               const conduit::Node &options,
               conduit::Node &output) const
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

// The following members are private (unless using CUDA)
#if !defined(__CUDACC__)
private:
#endif

  /**
   * \brief This struct contains information used when merging fields.
   */
  struct FieldInformation
  {
    std::string topology;
    std::string association;
    int dtype;
    std::vector<std::string> components;
  };

  /**
   * \brief Check that the mesh inputs are valid and meet constraints. There must
   *        be 1 coordset/topology/matset. The coordset must be explicit and the
   *        topology must be unstructured and for now, non-polyhedral.
   *
   * \param inputs The mesh inputs.
   *
   * \return True if the inputs appear to be valid; False otherwise.
   */
  bool validInputs(const std::vector<MeshInput> &inputs,
                   const conduit::Node &options) const
  {
    std::string topoName;
    if(options.has_child("topology"))
      topoName = options.fetch_existing("topology").as_string();

    for(size_t i = 0; i < inputs.size(); i++)
    {
      if(inputs[i].m_input == nullptr) return false;

      // If we did not specify which topology, make sure that there is only 1.
      const char *keys[] = {"coordsets", "topologies", "matsets"};
      if(topoName.empty())
      {
        for(int k = 0; k < 3; k++)
        {
          if(inputs[i].m_input->has_path(keys[k]))
          {
            const conduit::Node &n = inputs[i].m_input->fetch_existing(keys[k]);
            if(n.number_of_children() > 1) return false;
          }
        }
      }

      const conduit::Node &n_coordsets =
        inputs[i].m_input->fetch_existing("coordsets");
      const conduit::Node &n_coordset = n_coordsets[0];
      if(n_coordset["type"].as_string() != "explicit")
      {
        return false;
      }

      const conduit::Node &n_topologies =
        inputs[i].m_input->fetch_existing("topologies");
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
      if(!((inputs[i].m_nodeMapView.size() == 0 &&
            inputs[i].m_nodeSliceView.size() == 0) ||
           (inputs[i].m_nodeMapView.size() > 0 &&
            inputs[i].m_nodeSliceView.size() > 0)))
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
    mergeFields(inputs, output);
    mergeMatset(inputs, output);
  }

  /**
   * \brief Merge multiple coordsets into a single coordset. No node merging takes place.
   *
   * \param inputs A vector of inputs to be merged.
   * \param[out] output The node that will contain the output mesh.
   */
  void mergeCoordset(const std::vector<MeshInput> &inputs,
                     conduit::Node &output) const
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
      const conduit::Node &coordsets =
        inputs[i].m_input->fetch_existing("coordsets");
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
      axom::mir::views::FloatNode_to_ArrayView(n_srcValues[0], [&](auto comp0) {
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
            axom::for_all<ExecSpace>(
              nodeSliceView.size(),
              AXOM_LAMBDA(axom::IndexType index) {
                const auto sliceIndex = nodeSliceView[index];
                compView[offset + index] = srcCompView[sliceIndex];
              });
            size = nodeSliceView.size();
          }
          else
          {
            // Pull out all nodes from the input.
            axom::for_all<ExecSpace>(
              srcCompView.size(),
              AXOM_LAMBDA(axom::IndexType index) {
                compView[offset + index] = srcCompView[index];
              });
            size = srcCompView.size();
          }

          offsets[c] += size;
        }
      });
    }
  }

  axom::IndexType countNodes(const std::vector<MeshInput> &inputs,
                             size_t index) const
  {
    SLIC_ASSERT(index < inputs.size());

    const conduit::Node &coordsets =
      inputs[index].m_input->fetch_existing("coordsets");
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
      const auto nnodes = countNodes(inputs, i);
      nodeTotal += nnodes;
    }
    return nodeTotal;
  }

  axom::IndexType countZones(const std::vector<MeshInput> &inputs,
                             size_t index) const
  {
    const conduit::Node &n_topologies =
      inputs[index].m_input->fetch_existing("topologies");
    const conduit::Node &n_topo = n_topologies[0];

    const conduit::Node &n_size = n_topo.fetch_existing("elements/sizes");
    axom::IndexType nzones = n_size.dtype().number_of_elements();
    return nzones;
  }

  void countZones(const std::vector<MeshInput> &inputs,
                  axom::IndexType &totalConnLength,
                  axom::IndexType &totalZones) const
  {
    totalConnLength = 0;
    totalZones = 0;
    axom::IndexType n = static_cast<axom::IndexType>(inputs.size());
    for(axom::IndexType i = 0; i < n; i++)
    {
      const conduit::Node &n_topologies =
        inputs[i].m_input->fetch_existing("topologies");
      const conduit::Node &n_topo = n_topologies[0];
      const auto type = n_topo.fetch_existing("type").as_string();
      SLIC_ASSERT(type == "unstructured");

      const conduit::Node &n_conn =
        n_topo.fetch_existing("elements/connectivity");
      const auto connLength = n_conn.dtype().number_of_elements();
      totalConnLength += connLength;

      const conduit::Node &n_size = n_topo.fetch_existing("elements/sizes");
      const auto nzones = n_size.dtype().number_of_elements();
      totalZones += nzones;
    }
  }

  /**
   * \brief Merge multiple topologies into a single topology.
   *
   * \param inputs A vector of inputs to be merged.
   * \param[out] output The node that will contain the output mesh.
   */
  void mergeTopology(const std::vector<MeshInput> &inputs,
                     conduit::Node &output) const
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
      const conduit::Node &n_topologies =
        inputs[i].m_input->fetch_existing("topologies");
      const conduit::Node &n_srcTopo = n_topologies[0];
      const auto type = n_srcTopo.fetch_existing("type").as_string();
      const auto shape = n_srcTopo.fetch_existing("elements/shape").as_string();
      SLIC_ASSERT(type == "unstructured");
      if(shape == "mixed")
      {
        const conduit::Node &n_shape_map =
          n_srcTopo.fetch_existing("elements/shape_map");
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
    axom::IndexType connOffset = 0, sizesOffset = 0, shapesOffset = 0,
                    coordOffset = 0;
    for(axom::IndexType i = 0; i < n; i++)
    {
      const conduit::Node &n_topologies =
        inputs[i].m_input->fetch_existing("topologies");
      const conduit::Node &n_srcTopo = n_topologies[0];

      const std::string srcShape =
        n_srcTopo.fetch_existing("elements/shape").as_string();
      const conduit::Node &n_srcConn =
        n_srcTopo.fetch_existing("elements/connectivity");
      const conduit::Node &n_srcSizes =
        n_srcTopo.fetch_existing("elements/sizes");

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
      axom::mir::views::IndexNode_to_ArrayView(n_srcConn, [&](auto srcConnView) {
        using ConnType = typename decltype(srcConnView)::value_type;
        conduit::Node &n_newConn =
          n_newTopoPtr->fetch_existing("elements/connectivity");
        auto connView = bputils::make_array_view<ConnType>(n_newConn);

        if(inputs[i].m_nodeMapView.size() > 0)
        {
          // Copy all zones from the input but map the nodes to new values.
          // The supplied nodeMap is assumed to be a mapping from the current
          // node connectivity to the merged node connectivity.
          const auto nodeMapView = inputs[i].m_nodeMapView;
          axom::for_all<ExecSpace>(
            srcConnView.size(),
            AXOM_LAMBDA(axom::IndexType index) {
              const auto nodeId = srcConnView[index];
              const auto newNodeId = nodeMapView[nodeId];
              connView[connOffset + index] = newNodeId;
            });
        }
        else
        {
          // Copy all zones from the input. Map the nodes to the new values.
          axom::for_all<ExecSpace>(
            srcConnView.size(),
            AXOM_LAMBDA(axom::IndexType index) {
              connView[connOffset + index] = coordOffset + srcConnView[index];
            });
        }
        connOffset += srcConnView.size();
        coordOffset += countNodes(inputs, static_cast<size_t>(i));
      });

      // Copy this input's sizes into the new topology.
      axom::mir::views::IndexNode_to_ArrayView(n_srcSizes, [&](auto srcSizesView) {
        using ConnType = typename decltype(srcSizesView)::value_type;
        conduit::Node &n_newSizes =
          n_newTopoPtr->fetch_existing("elements/sizes");
        auto sizesView = bputils::make_array_view<ConnType>(n_newSizes);

        // Copy all sizes from the input.
        axom::for_all<ExecSpace>(
          srcSizesView.size(),
          AXOM_LAMBDA(axom::IndexType index) {
            sizesView[sizesOffset + index] = srcSizesView[index];
          });

        sizesOffset += srcSizesView.size();
      });

      // Copy shape information if it exists.
      if(n_srcTopo.has_path("elements/shapes"))
      {
        const conduit::Node &n_srcShapes =
          n_srcTopo.fetch_existing("elements/shapes");

        axom::mir::views::IndexNode_to_ArrayView(
          n_srcShapes,
          [&](auto srcShapesView) {
            using ConnType = typename decltype(srcShapesView)::value_type;
            conduit::Node &n_newShapes =
              n_newTopoPtr->fetch_existing("elements/shapes");
            auto shapesView = bputils::make_array_view<ConnType>(n_newShapes);

            // Copy all sizes from the input.
            axom::for_all<ExecSpace>(
              srcShapesView.size(),
              AXOM_LAMBDA(axom::IndexType index) {
                shapesView[shapesOffset + index] = srcShapesView[index];
              });

            shapesOffset += srcShapesView.size();
          });
      }
      else
      {
        // Fill in shape information. There is no source shape information. Use
        // sizes to get the number of zones.
        const conduit::Node &n_srcSizes =
          n_srcTopo.fetch_existing("elements/sizes");
        axom::IndexType nz = n_srcSizes.dtype().number_of_elements();
        conduit::Node &n_newShapes =
          n_newTopoPtr->fetch_existing("elements/shapes");
        axom::mir::views::IndexNode_to_ArrayView(n_newShapes, [&](auto shapesView) {
          const int shapeId = axom::mir::views::shapeNameToID(srcShape);
          axom::for_all<ExecSpace>(
            nz,
            AXOM_LAMBDA(axom::IndexType index) {
              shapesView[shapesOffset + index] = shapeId;
            });
          shapesOffset += nz;
        });
      }
    }

    // Make new offsets from the sizes.
    conduit::Node &n_newSizes = n_newTopoPtr->fetch_existing("elements/sizes");
    axom::mir::views::IndexNode_to_ArrayView(n_newSizes, [&](auto sizesView) {
      using ConnType = typename decltype(sizesView)::value_type;
      conduit::Node &n_newOffsets =
        n_newTopoPtr->fetch_existing("elements/offsets");
      auto offsetsView = bputils::make_array_view<ConnType>(n_newOffsets);
      axom::exclusive_scan<ExecSpace>(sizesView, offsetsView);
    });
  }

  /**
   * \brief Merge fields that exist on the various mesh inputs. Zero-fill values
   *        where a field does not exist in an input.
   *
   * \param inputs A vector of inputs to be merged.
   * \param[out] output The node that will contain the output mesh.
   */
  void mergeFields(const std::vector<MeshInput> &inputs, conduit::Node &output) const
  {
    namespace bputils = axom::mir::utilities::blueprint;
    AXOM_ANNOTATE_SCOPE("mergeFields");
    axom::IndexType totalNodes = countNodes(inputs);
    axom::IndexType totalConnLen = 0, totalZones = 0;
    countZones(inputs, totalConnLen, totalZones);
    const axom::IndexType n = static_cast<axom::IndexType>(inputs.size());

    // Determine whether any inputs have fields.
    bool hasFields = false;
    for(axom::IndexType i = 0; i < n; i++)
    {
      hasFields |= inputs[i].m_input->has_child("fields");
    }

    if(hasFields)
    {
      // Make field information in case some inputs do not have the field.
      std::map<std::string, FieldInformation> fieldInfo;
      for(axom::IndexType i = 0; i < n; i++)
      {
        if(inputs[i].m_input->has_child("fields"))
        {
          const conduit::Node &n_fields =
            inputs[i].m_input->fetch_existing("fields");
          for(conduit::index_t c = 0; c < n_fields.number_of_children(); c++)
          {
            const conduit::Node &n_field = n_fields[c];
            const conduit::Node &n_values = n_field.fetch_existing("values");
            FieldInformation fi;
            fi.topology = n_field.fetch_existing("topology").as_string();
            fi.association = n_field.fetch_existing("association").as_string();
            if(n_values.number_of_children() > 0)
            {
              for(conduit::index_t comp = 0; comp < n_values.number_of_children();
                  comp++)
              {
                fi.components.push_back(n_values[comp].name());
                fi.dtype = n_values[comp].dtype().id();
              }
            }
            else
            {
              fi.dtype = n_values.dtype().id();
            }
            fieldInfo[n_field.name()] = fi;
          }
        }
      }

      // Make new fields
      bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;
      conduit::Node &n_newFields = output["fields"];
      for(auto it = fieldInfo.begin(); it != fieldInfo.end(); it++)
      {
        conduit::Node &n_newField = n_newFields[it->first];
        n_newField["association"] = it->second.association;
        n_newField["topology"] = it->second.topology;
        conduit::Node &n_values = n_newField["values"];
        if(it->second.components.empty())
        {
          // Scalar
          conduit::Node &n_values = n_newField["values"];
          n_values.set_allocator(c2a.getConduitAllocatorID());
          const std::string srcPath("fields/" + it->first + "/values");
          if(it->second.association == "element")
          {
            n_values.set(conduit::DataType(it->second.dtype, totalZones));
            copyZonal(inputs, n_values, srcPath);
          }
          else if(it->second.association == "vertex")
          {
            n_values.set(conduit::DataType(it->second.dtype, totalNodes));
            copyNodal(inputs, n_values, srcPath);
          }
        }
        else
        {
          // Vector
          for(size_t ci = 0; ci < it->second.components.size(); ci++)
          {
            conduit::Node &n_comp = n_values[it->second.components[ci]];
            n_comp.set_allocator(c2a.getConduitAllocatorID());
            const std::string srcPath("fields/" + it->first + "/values/" +
                                      it->second.components[ci]);
            if(it->second.association == "element")
            {
              n_comp.set(conduit::DataType(it->second.dtype, totalZones));
              copyZonal(inputs, n_comp, srcPath);
            }
            else if(it->second.association == "vertex")
            {
              n_comp.set(conduit::DataType(it->second.dtype, totalNodes));
              copyNodal(inputs, n_comp, srcPath);
            }
          }
        }
      }
    }
  }

  /**
   * \brief Copy zonal field data into a Conduit node.
   *
   * \param inputs A vector of inputs to be merged.
   * \param[out] n_values The node will be populated with data values from the field inputs.
   * \param srcPath The path to the source data in each input node.
   */
  void copyZonal(const std::vector<MeshInput> &inputs,
                 conduit::Node &n_values,
                 const std::string &srcPath) const
  {
    axom::IndexType offset = 0;
    for(size_t i = 0; i < inputs.size(); i++)
    {
      const auto nzones = countZones(inputs, i);

      if(inputs[i].m_input->has_path(srcPath))
      {
        const conduit::Node &n_src_values =
          inputs[i].m_input->fetch_existing(srcPath);
        axom::mir::views::Node_to_ArrayView(n_src_values,
                                            n_values,
                                            [&](auto srcView, auto destView) {
                                              axom::for_all<ExecSpace>(
                                                nzones,
                                                AXOM_LAMBDA(axom::IndexType index) {
                                                  destView[offset + index] =
                                                    srcView[index];
                                                });
                                            });
      }
      else
      {
        axom::mir::views::Node_to_ArrayView(n_values, [&](auto destView) {
          axom::for_all<ExecSpace>(
            nzones,
            AXOM_LAMBDA(axom::IndexType index) { destView[offset + index] = 0; });
        });
      }
      offset += nzones;
    }
  }

#if 0
#define AXOM_NODE_TO_ARRAYVIEW1(node1, view1, CODE) \
axom::mir::views::Node_to_ArrayView(node1, [&](auto view1){CODE});
#else
#define AXOM_NODE_TO_ARRAYVIEW1(node1, view1, CODE) \
  switch(node1.dtype().id()) \
  { \
  case conduit::DataType::INT8_ID: \
   { axom::ArrayView<conduit::int8> view1(static_cast<conduit::int8 *>(node1.data_ptr()), node1.dtype().number_of_elements()); \
     CODE \
   } break; \
  case conduit::DataType::INT16_ID: { \
     axom::ArrayView<conduit::int16> view1(static_cast<conduit::int16 *>(node1.data_ptr()), node1.dtype().number_of_elements()); \
     CODE \
   } break; \
  case conduit::DataType::INT32_ID: \
   { axom::ArrayView<conduit::int32> view1(static_cast<conduit::int32 *>(node1.data_ptr()), node1.dtype().number_of_elements()); \
     CODE \
   } break; \
  case conduit::DataType::INT64_ID: { \
     axom::ArrayView<conduit::int64> view1(static_cast<conduit::int64 *>(node1.data_ptr()), node1.dtype().number_of_elements()); \
     CODE \
   } break; \
  case conduit::DataType::UINT8_ID: \
   { axom::ArrayView<conduit::uint8> view1(static_cast<conduit::uint8 *>(node1.data_ptr()), node1.dtype().number_of_elements()); \
     CODE \
   } break; \
  case conduit::DataType::UINT16_ID: { \
     axom::ArrayView<conduit::uint16> view1(static_cast<conduit::uint16 *>(node1.data_ptr()), node1.dtype().number_of_elements()); \
     CODE \
   } break; \
  case conduit::DataType::UINT32_ID: \
   { axom::ArrayView<conduit::uint32> view1(static_cast<conduit::uint32 *>(node1.data_ptr()), node1.dtype().number_of_elements()); \
     CODE \
   } break; \
  case conduit::DataType::UINT64_ID: { \
     axom::ArrayView<conduit::uint64> view1(static_cast<conduit::uint64 *>(node1.data_ptr()), node1.dtype().number_of_elements()); \
     CODE \
   } break; \
  case conduit::DataType::FLOAT32_ID: \
   { axom::ArrayView<conduit::float32> view1(static_cast<conduit::float32 *>(node1.data_ptr()), node1.dtype().number_of_elements()); \
     CODE \
   } break; \
  case conduit::DataType::FLOAT64_ID: { \
     axom::ArrayView<conduit::float64> view1(static_cast<conduit::float64 *>(node1.data_ptr()), node1.dtype().number_of_elements()); \
     CODE \
   } break; \
  }

#define AXOM_NODE_TO_ARRAYVIEW_SAME2(node1, node2, view1, view2, CODE) \
  switch(node1.dtype().id()) \
  { \
  case conduit::DataType::INT8_ID: {\
     axom::ArrayView<conduit::int8> view1(static_cast<conduit::int8 *>(const_cast<void *>(node1.data_ptr())), node1.dtype().number_of_elements()); \
     axom::ArrayView<conduit::int8> view2(static_cast<conduit::int8 *>(const_cast<void *>(node2.data_ptr())), node2.dtype().number_of_elements()); \
     CODE \
   } break; \
  case conduit::DataType::INT16_ID: { \
     axom::ArrayView<conduit::int16> view1(static_cast<conduit::int16 *>(const_cast<void *>(node1.data_ptr())), node1.dtype().number_of_elements()); \
     axom::ArrayView<conduit::int16> view2(static_cast<conduit::int16 *>(const_cast<void *>(node2.data_ptr())), node2.dtype().number_of_elements()); \
     CODE \
   } break; \
  case conduit::DataType::INT32_ID: {\
     axom::ArrayView<conduit::int32> view1(static_cast<conduit::int32 *>(const_cast<void *>(node1.data_ptr())), node1.dtype().number_of_elements()); \
     axom::ArrayView<conduit::int32> view2(static_cast<conduit::int32 *>(const_cast<void *>(node2.data_ptr())), node2.dtype().number_of_elements()); \
     CODE \
   } break; \
  case conduit::DataType::INT64_ID: { \
     axom::ArrayView<conduit::int64> view1(static_cast<conduit::int64 *>(const_cast<void *>(node1.data_ptr())), node1.dtype().number_of_elements()); \
     axom::ArrayView<conduit::int64> view2(static_cast<conduit::int64 *>(const_cast<void *>(node2.data_ptr())), node2.dtype().number_of_elements()); \
     CODE \
   } break; \
  case conduit::DataType::UINT8_ID: {\
     axom::ArrayView<conduit::uint8> view1(static_cast<conduit::uint8 *>(const_cast<void *>(node1.data_ptr())), node1.dtype().number_of_elements()); \
     axom::ArrayView<conduit::uint8> view2(static_cast<conduit::uint8 *>(const_cast<void *>(node2.data_ptr())), node2.dtype().number_of_elements()); \
     CODE \
   } break; \
  case conduit::DataType::UINT16_ID: { \
     axom::ArrayView<conduit::uint16> view1(static_cast<conduit::uint16 *>(const_cast<void *>(node1.data_ptr())), node1.dtype().number_of_elements()); \
     axom::ArrayView<conduit::uint16> view2(static_cast<conduit::uint16 *>(const_cast<void *>(node2.data_ptr())), node2.dtype().number_of_elements()); \
     CODE \
   } break; \
  case conduit::DataType::UINT32_ID: { \
     axom::ArrayView<conduit::uint32> view1(static_cast<conduit::uint32 *>(const_cast<void *>(node1.data_ptr())), node1.dtype().number_of_elements()); \
     axom::ArrayView<conduit::uint32> view2(static_cast<conduit::uint32 *>(const_cast<void *>(node2.data_ptr())), node2.dtype().number_of_elements()); \
     CODE \
   } break; \
  case conduit::DataType::UINT64_ID: { \
     axom::ArrayView<conduit::uint64> view1(static_cast<conduit::uint64 *>(const_cast<void *>(node1.data_ptr())), node1.dtype().number_of_elements()); \
     axom::ArrayView<conduit::uint64> view2(static_cast<conduit::uint64 *>(const_cast<void *>(node2.data_ptr())), node2.dtype().number_of_elements()); \
     CODE \
   } break; \
  case conduit::DataType::FLOAT32_ID: {\
     axom::ArrayView<conduit::float32> view1(static_cast<conduit::float32 *>(const_cast<void *>(node1.data_ptr())), node1.dtype().number_of_elements()); \
     axom::ArrayView<conduit::float32> view2(static_cast<conduit::float32 *>(const_cast<void *>(node2.data_ptr())), node2.dtype().number_of_elements()); \
     CODE \
   } break; \
  case conduit::DataType::FLOAT64_ID: { \
     axom::ArrayView<conduit::float64> view1(static_cast<conduit::float64 *>(const_cast<void *>(node1.data_ptr())), node1.dtype().number_of_elements()); \
     axom::ArrayView<conduit::float64> view2(static_cast<conduit::float64 *>(const_cast<void *>(node2.data_ptr())), node2.dtype().number_of_elements()); \
     CODE \
   } break; \
  }
#endif
  /**
   * \brief Copy nodal field data into a Conduit node.
   *
   * \param inputs A vector of inputs to be merged.
   * \param[out] n_values The node will be populated with data values from the field inputs.
   * \param srcPath The path to the source data in each input node.
   */
  void copyNodal(const std::vector<MeshInput> &inputs,
                 conduit::Node &n_values,
                 const std::string &srcPath) const
  {
    axom::IndexType offset = 0;
    for(size_t i = 0; i < inputs.size(); i++)
    {
      const auto nnodes = countNodes(inputs, i);

      if(inputs[i].m_input->has_path(srcPath))
      {
        const conduit::Node &n_src_values =
          inputs[i].m_input->fetch_existing(srcPath);

        axom::mir::views::Node_to_ArrayView(
          n_src_values,
          n_values,
          [&](auto srcView, auto destView) {

          copyNodal_copy(inputs[i].m_nodeSliceView, srcView, destView, nnodes, offset);

          });
      }
      else
      {
#if 1
        axom::mir::views::Node_to_ArrayView(n_values, [&](auto destView) {
          copyNodal_fill(destView, nnodes, offset);
        });
#else
        axom::mir::views::Node_to_ArrayView(n_values, [&](auto destView) {
          axom::for_all<ExecSpace>(
            nnodes,
            AXOM_LAMBDA(axom::IndexType index) { destView[offset + index] = 0; });
        });
#endif
      }
      offset += nnodes;
    }
  }

  template <typename SrcViewType, typename DestViewType>
  void copyNodal_copy(axom::ArrayView<axom::IndexType> nodeSliceView, SrcViewType srcView, DestViewType destView, axom::IndexType nnodes, axom::IndexType offset) const
  {
            if(nodeSliceView.empty())
            {
              axom::for_all<ExecSpace>(
                nnodes,
                AXOM_LAMBDA(axom::IndexType index) {
                  destView[offset + index] = srcView[index];
                });
            }
            else
            {
              axom::for_all<ExecSpace>(
                nnodes,
                AXOM_LAMBDA(axom::IndexType index) {
                  const auto nodeId = nodeSliceView[index];
                  destView[offset + index] = srcView[nodeId];
                });
            }
   }

  template <typename DestViewType>
  void copyNodal_fill(DestViewType destView, axom::IndexType nnodes, axom::IndexType offset) const
  {
    axom::for_all<ExecSpace>(
           nnodes,
           AXOM_LAMBDA(axom::IndexType index) { destView[offset + index] = 0; });
  }

  /**
   * \brief Merge matsets that exist on the various mesh inputs.
   *
   * \param inputs A vector of inputs to be merged.
   * \param[out] output The node that will contain the output mesh.
   */
  void mergeMatset(const std::vector<MeshInput> &inputs, conduit::Node &output) const
  {
    AXOM_ANNOTATE_SCOPE("mergeMatset");
    namespace bputils = axom::mir::utilities::blueprint;
    using reduce_policy =
      typename axom::execution_space<ExecSpace>::reduce_policy;
    bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;

    // Make a pass through the inputs and make a list of the material names.
    bool hasMatsets = false, defaultMaterial = false;
    int nmats = 0;
    std::map<std::string, int> allMats;
    std::string matsetName, topoName;
    for(size_t i = 0; i < inputs.size(); i++)
    {
      if(inputs[i].m_input->has_path("matsets"))
      {
        conduit::Node &n_matsets = inputs[i].m_input->fetch_existing("matsets");
        conduit::Node &n_matset = n_matsets[0];
        matsetName = n_matset.name();
        topoName = n_matset.fetch_existing("topology").as_string();
        auto matInfo = axom::mir::views::materials(n_matset);
        for(const auto &info : matInfo)
        {
          if(allMats.find(info.name) == allMats.end())
          {
            allMats[info.name] = nmats++;
          }
        }
        hasMatsets = true;
      }
      else
      {
        defaultMaterial = true;
      }
    }

    if(hasMatsets)
    {
      // One or more inputs did not have a matset.
      if(defaultMaterial) allMats["default"] = nmats++;

      // Make a pass through the matsets to determine the overall storage.
      axom::IndexType totalZones = 0, totalMatCount = 0;
      int itype, ftype;
      {
        AXOM_ANNOTATE_SCOPE("sizes");
        for(size_t i = 0; i < inputs.size(); i++)
        {
          const auto nzones = countZones(inputs, i);
          totalZones += nzones;

          if(inputs[i].m_input->has_path("matsets"))
          {
            conduit::Node &n_matsets =
              inputs[i].m_input->fetch_existing("matsets");
            conduit::Node &n_matset = n_matsets[0];
            axom::IndexType matCount = 0;
            axom::mir::views::dispatch_material(n_matset, [&](auto matsetView) {
              // Figure out the types to use for storing the data.
              using IType = typename decltype(matsetView)::IndexType;
              using FType = typename decltype(matsetView)::FloatType;
              itype = bputils::cpp2conduit<IType>::id;
              ftype = bputils::cpp2conduit<FType>::id;

              RAJA::ReduceSum<reduce_policy, axom::IndexType> matCount_reduce(0);
              axom::for_all<ExecSpace>(
                nzones,
                AXOM_LAMBDA(axom::IndexType zoneIndex) {
                  const auto nmats = matsetView.numberOfMaterials(zoneIndex);
                  matCount_reduce += nmats;
                });
              matCount = matCount_reduce.get();
            });
            totalMatCount += matCount;
          }
          else
          {
            totalMatCount += nzones;
          }
        }
      }

      // Allocate
      AXOM_ANNOTATE_BEGIN("allocate");
      conduit::Node &n_newMatset = output["matsets/" + matsetName];
      n_newMatset["topology"] = topoName;
      conduit::Node &n_volume_fractions = n_newMatset["volume_fractions"];
      n_volume_fractions.set_allocator(c2a.getConduitAllocatorID());
      n_volume_fractions.set(conduit::DataType(ftype, totalMatCount));

      conduit::Node &n_material_ids = n_newMatset["material_ids"];
      n_material_ids.set_allocator(c2a.getConduitAllocatorID());
      n_material_ids.set(conduit::DataType(itype, totalMatCount));

      conduit::Node &n_sizes = n_newMatset["sizes"];
      n_sizes.set_allocator(c2a.getConduitAllocatorID());
      n_sizes.set(conduit::DataType(itype, totalZones));

      conduit::Node &n_offsets = n_newMatset["offsets"];
      n_offsets.set_allocator(c2a.getConduitAllocatorID());
      n_offsets.set(conduit::DataType(itype, totalZones));

      conduit::Node &n_indices = n_newMatset["indices"];
      n_indices.set_allocator(c2a.getConduitAllocatorID());
      n_indices.set(conduit::DataType(itype, totalMatCount));
      AXOM_ANNOTATE_END("allocate");

      {
        AXOM_ANNOTATE_SCOPE("populate");

        // Make material_map.
        conduit::Node &n_material_map = n_newMatset["material_map"];
        for(auto it = allMats.begin(); it != allMats.end(); it++)
          n_material_map[it->first] = it->second;

        // Populate
        axom::mir::views::IndexNode_to_ArrayView_same(
          n_material_ids,
          n_sizes,
          n_offsets,
          n_indices,
          [&](auto materialIdsView,
              auto sizesView,
              auto offsetsView,
              auto indicesView) {
            axom::mir::views::FloatNode_to_ArrayView(
              n_volume_fractions,
              [&](auto volumeFractionsView) {
                // Fill in sizes array.
                axom::IndexType zOffset = 0;
                for(size_t i = 0; i < inputs.size(); i++)
                {
                  const auto nzones = countZones(inputs, i);

                  if(inputs[i].m_input->has_child("matsets"))
                  {
                    conduit::Node &n_matsets =
                      inputs[i].m_input->fetch_existing("matsets");
                    conduit::Node &n_matset = n_matsets[0];

                    axom::mir::views::dispatch_material(
                      n_matset,
                      [&](auto matsetView) {
                        axom::for_all<ExecSpace>(
                          nzones,
                          AXOM_LAMBDA(axom::IndexType zoneIndex) {
                            sizesView[zOffset + zoneIndex] =
                              matsetView.numberOfMaterials(zoneIndex);
                          });
                      });
                  }
                  else
                  {
                    axom::for_all<ExecSpace>(
                      nzones,
                      AXOM_LAMBDA(axom::IndexType zoneIndex) {
                        sizesView[zOffset + zoneIndex] = 1;
                      });
                  }
                  zOffset += nzones;
                }
                // Make offsets.
                axom::exclusive_scan<ExecSpace>(sizesView, offsetsView);

                // Make indices.
                axom::for_all<ExecSpace>(
                  totalMatCount,
                  AXOM_LAMBDA(axom::IndexType index) { indicesView[index] = index; });

                // Fill in material info.
                zOffset = 0;
                for(size_t i = 0; i < inputs.size(); i++)
                {
                  const auto nzones = countZones(inputs, i);

                  if(inputs[i].m_input->has_child("matsets"))
                  {
                    conduit::Node &n_matsets =
                      inputs[i].m_input->fetch_existing("matsets");
                    conduit::Node &n_matset = n_matsets[0];

                    axom::mir::views::dispatch_material(
                      n_matset,
                      [&](auto matsetView) {
                        using IDList = typename decltype(matsetView)::IDList;
                        using VFList = typename decltype(matsetView)::VFList;
                        using MatID = typename decltype(matsetView)::IndexType;

                        // Make some maps for renumbering material numbers.
                        const auto localMaterialMap =
                          axom::mir::views::materials(n_matset);
                        std::map<MatID, MatID> localToAll;
                        for(const auto &info : localMaterialMap)
                        {
                          MatID matno = allMats[info.name];
                          localToAll[info.number] = matno;
                        }
                        std::vector<MatID> localVec, allVec;
                        for(auto it = localToAll.begin(); it != localToAll.end();
                            it++)
                        {
                          localVec.push_back(it->first);
                          allVec.push_back(it->second);
                        }
                        // Put maps on device.
                        const int allocatorID =
                          axom::execution_space<ExecSpace>::allocatorID();
                        axom::Array<MatID> local(localVec.size(),
                                                 localVec.size(),
                                                 allocatorID);
                        axom::Array<MatID> all(allVec.size(),
                                               allVec.size(),
                                               allocatorID);
                        axom::copy(local.data(),
                                   localVec.data(),
                                   sizeof(MatID) * local.size());
                        axom::copy(all.data(),
                                   allVec.data(),
                                   sizeof(MatID) * all.size());
                        const auto localView = local.view();
                        const auto allView = all.view();

                        axom::for_all<ExecSpace>(
                          nzones,
                          AXOM_LAMBDA(axom::IndexType zoneIndex) {
                            // Get this zone's materials.
                            IDList ids;
                            VFList vfs;
                            matsetView.zoneMaterials(zoneIndex, ids, vfs);

                            // Store the materials in the new material.
                            const auto zoneStart =
                              offsetsView[zOffset + zoneIndex];
                            for(axom::IndexType mi = 0; mi < ids.size(); mi++)
                            {
                              const auto destIndex = zoneStart + mi;
                              volumeFractionsView[destIndex] = vfs[mi];

                              // Get the index of the material number in the local map.
                              const auto mapIndex =
                                axom::mir::utilities::bsearch(ids[mi], localView);
                              assert(mapIndex != -1);
                              // We'll store the all materials number.
                              const auto allMatno = allView[mapIndex];
                              materialIdsView[destIndex] = allMatno;
                            }
                          });
                      });
                  }
                  else
                  {
                    const int dmat = allMats["default"];
                    axom::for_all<ExecSpace>(
                      nzones,
                      AXOM_LAMBDA(axom::IndexType zoneIndex) {
                        const auto zoneStart = offsetsView[zOffset + zoneIndex];
                        volumeFractionsView[zoneStart] = 1;
                        materialIdsView[zoneStart] = dmat;
                      });
                  }
                  zOffset += nzones;
                }
              });
          });
      }
    }  // if hasMatsets
  }
};

}  // end namespace blueprint
}  // end namespace utilities
}  // end namespace mir
}  // end namespace axom

#endif
