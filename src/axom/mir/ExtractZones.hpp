// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_EXTRACT_ZONES_HPP
#define AXOM_MIR_EXTRACT_ZONES_HPP

#include <axom/core.hpp>
#include <axom/mir.hpp>
#include <axom/mir/MatsetSlicer.hpp> // Needed to get MatsetSlicer

namespace axom
{
namespace mir
{
namespace utilities
{
namespace blueprint
{

/**
 * \brief Make a new topology and coordset by extracting certain zones from the input mesh.
 *
 * \tparam ExecSpace The execution space where the algorithm will run.
 * \tparam TopologyView The topology view type on which the algorithm will run.
 */
template <typename ExecSpace, typename TopologyView, typename CoordsetView>
class ExtractZones
{
  using ConnectivityType = typename TopologyView::ConnectivityType;
  using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;
public:
  using SelectedZonesView = axom::ArrayView<axom::IndexType>;

  /**
   * \brief Constructor
   *
   * \param topoView The input topology view.
   * \param coordsetView The input coordset view.
   * \param matsetView The input matset view.
   */
  ExtractZones(const TopologyView &topoView, const CoordsetView &coordsetView) :
    m_topologyView(topoView), m_coordsetView(coordsetView)
  {
  }

  /**
   * \brief Select zones from the input mesh by id and output them in the output mesh.
   *
   * \param selectedZonesView A view that contains the selected zone ids.
   * \param n_input The input mesh.
   * \param n_options The input options.
   * \param[out] n_output The output mesh.
   *
   * \note The \a n_options node contains a "topology" string that is selects the
   *       name of the topology to extract.
   */
  void execute(const SelectedZonesView &selectedZonesView,
               const conduit::Node &n_input,
               const conduit::Node &n_options,
               conduit::Node &n_output) const
  {
    AXOM_ANNOTATE_SCOPE("ExtractZones");
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    AXOM_ANNOTATE_BEGIN("nodeMap");
    // We need to figure out which nodes to keep.
    const auto nnodes = m_coordsetView.numberOfNodes();
    axom::Array<int> mask(nnodes, nnodes, allocatorID);
    auto maskView = mask.view();
    mask.fill(0);

    // Mark all the selected zones' nodes as 1. Multiple threads may write 1 to the same node.
    RAJA::ReduceSum<reduce_policy, int> connsize_reduce(0);
    m_topologyView.template for_selected_zones<ExecSpace>(selectedZonesView, AXOM_LAMBDA(auto AXOM_UNUSED_PARAM(szIndex), auto AXOM_UNUSED_PARAM(zoneIndex), const auto &zone)
    {
      const int nids = zone.numberOfNodes();
      for(int i = 0; i < nids; i++)
      {
        const auto nodeId = zone.getId(i);
        maskView[nodeId] = 1;
      }
      connsize_reduce += nids;
    });

    // Count the used nodes.
    RAJA::ReduceSum<reduce_policy, int> mask_reduce(0);
    axom::for_all<ExecSpace>(nnodes, AXOM_LAMBDA(auto index)
    {
      mask_reduce += maskView[index];
    });
    const int newNumNodes = mask_reduce.get();

    // Make a compact list of nodes.
    axom::Array<int> maskOffsets(nnodes, nnodes, allocatorID);
    auto maskOffsetsView = maskOffsets.view();
    axom::exclusive_scan<ExecSpace>(maskView, maskOffsetsView);

    // Make an array of original node ids that we can use to "slice" the nodal data.
    axom::Array<ConnectivityType> old2new(nnodes, nnodes, allocatorID);
    axom::Array<axom::IndexType> nodeSlice(newNumNodes, newNumNodes, allocatorID);
    auto old2newView = old2new.view();
    auto nodeSliceView = nodeSlice.view();
    axom::for_all<ExecSpace>(nnodes, AXOM_LAMBDA(auto index)
    {
      if(maskView[index] > 0)
      {
        nodeSliceView[maskOffsetsView[index]] = index;
        old2newView[index] = maskOffsetsView[index];
      }
    });
    AXOM_ANNOTATE_END("nodeMap");

    // Now make a new output topology.
    const conduit::Node &n_topologies = n_input.fetch_existing("topologies");
    const std::string topoName = topologyName(n_input, n_options);
    const conduit::Node &n_topo = n_topologies.fetch_existing(topoName);

    conduit::Node &n_newTopo = n_output["topologies/" + topoName];
    const auto newConnSize = connsize_reduce.get();
    makeTopology(selectedZonesView, newConnSize, old2newView, n_topo, n_newTopo);

    // Make a new coordset.
    SliceData nSlice;
    nSlice.m_indicesView = nodeSliceView;
    const std::string coordsetName = n_topo.fetch_existing("coordset").as_string();
    const conduit::Node &n_coordset = n_input.fetch_existing("coordsets/" + coordsetName);
    conduit::Node &n_newCoordset = n_output["coordsets/" + coordsetName];
    makeCoordset(nSlice, n_coordset, n_newCoordset);

    // Make new fields.
    if(n_input.has_child("fields"))
    {
      const conduit::Node &n_fields = n_input.fetch_existing("fields");
      conduit::Node &n_newFields = n_output["fields"];
      SliceData zSlice;
      zSlice.m_indicesView = selectedZonesView;
      makeFields(nSlice, zSlice, n_fields, n_newFields);
    }
  }

protected:
  /**
   * \brief Make the output topology for just the selected zones.
   *
   * \param selectedZonesView A view that contains the ids of the zones to extract.
   * \param newConnSize The expected size of the new connectivity.
   * \param old2newView A view that lets us map old node numbers to new node numbers.
   * \param n_topo The input topology.
   * \param n_newTopo A node to contain the new topology.
   */
  void makeTopology(const SelectedZonesView &selectedZonesView,
                    int newConnSize,
                    const axom::ArrayView<ConnectivityType> &old2newView,
                    const conduit::Node &n_topo,
                    conduit::Node &n_newTopo) const
  {
    AXOM_ANNOTATE_SCOPE("makeTopology");
    namespace bputils = axom::mir::utilities::blueprint;
    bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;

    const std::string shape = outputShape(n_topo);
    if(shape == "polyhedron")
    {
      // TODO: Handle polyhedron shape.
      // NOTE: We could know whether we can have PH topos if the TopologyView handles PH zones. Maybe this method is separated out and partially specialized.
      SLIC_ERROR("Polyhedron is not handled yet.");
    }
    else
    {
      // The topology is not polyhedral. We'll output unstructured.
      n_newTopo["type"] = "unstructured";
      n_newTopo["coordset"] = n_topo["coordset"].as_string();
      n_newTopo["elements/shape"] = outputShape(n_topo);

      conduit::Node &n_conn = n_newTopo["elements/connectivity"];
      n_conn.set_allocator(c2a.getConduitAllocatorID());
      n_conn.set(conduit::DataType(cpp2conduit<ConnectivityType>::id, newConnSize));
      auto connView = bputils::make_array_view<ConnectivityType>(n_conn);

      conduit::Node &n_sizes = n_newTopo["elements/sizes"];
      n_sizes.set_allocator(c2a.getConduitAllocatorID());
      n_sizes.set(conduit::DataType(cpp2conduit<ConnectivityType>::id, selectedZonesView.size()));
      auto sizesView = bputils::make_array_view<ConnectivityType>(n_sizes);

      conduit::Node &n_offsets = n_newTopo["elements/offsets"];
      n_offsets.set_allocator(c2a.getConduitAllocatorID());
      n_offsets.set(conduit::DataType(cpp2conduit<ConnectivityType>::id, selectedZonesView.size()));
      auto offsetsView = bputils::make_array_view<ConnectivityType>(n_offsets);

      // Fill in sizes, offsets, connectivity.
      m_topologyView.template for_selected_zones<ExecSpace>(selectedZonesView, AXOM_LAMBDA(auto szIndex, auto AXOM_UNUSED_PARAM(zoneIndex), const auto &zone)
      {
        sizesView[szIndex] = zone.numberOfNodes();
      });
      axom::exclusive_scan<ExecSpace>(sizesView, offsetsView);
      const axom::ArrayView<ConnectivityType> deviceOld2NewView(old2newView);
      m_topologyView.template for_selected_zones<ExecSpace>(selectedZonesView, AXOM_LAMBDA(auto szIndex, auto AXOM_UNUSED_PARAM(zoneIndex), const auto &zone)
      {
        const int size = static_cast<int>(sizesView[szIndex]);
        const auto offset = offsetsView[szIndex];
        for(int i = 0; i < size; i++)
        {
          const auto oldNodeId = zone.getId(i);
          const auto newNodeId = deviceOld2NewView[oldNodeId];
          connView[offset + i] = newNodeId;
        }
      });

      // Handle shapes, if present.
      if(n_topo.has_path("elements/shapes"))
      {
        const conduit::Node &n_shapes = n_topo.fetch_existing("elements/shapes");
        auto shapesView = bputils::make_array_view<ConnectivityType>(n_shapes);

        conduit::Node &n_newShapes = n_newTopo["elements/shapes"];
        n_newShapes.set_allocator(c2a.getConduitAllocatorID());
        n_newShapes.set(conduit::DataType(cpp2conduit<ConnectivityType>::id, selectedZonesView.size()));
        auto newShapesView = bputils::make_array_view<ConnectivityType>(n_newShapes);

        const SelectedZonesView deviceSelectedZonesView(selectedZonesView);
        axom::for_all<ExecSpace>(selectedZonesView.size(), AXOM_LAMBDA(auto index)
        {
          newShapesView[index] = shapesView[deviceSelectedZonesView[index]];
        });
      }
    }
  }

  /**
   * \brief Make the new coordset using the blend data and the input coordset/coordsetview.
   *
   * \param blend The BlendData that we need to construct the new coordset.
   * \param n_coordset The input coordset, which is passed for metadata.
   * \param[out] n_newCoordset The new coordset.
   */
  void makeCoordset(const SliceData &slice,
                    const conduit::Node &n_coordset,
                    conduit::Node &n_newCoordset) const
  {
    AXOM_ANNOTATE_SCOPE("makeCoordset");
    axom::mir::utilities::blueprint::CoordsetSlicer<ExecSpace, CoordsetView> cs(m_coordsetView);
    n_newCoordset.reset();
    cs.execute(slice, n_coordset, n_newCoordset);
  }

  /**
   * \brief Make fields for the output mesh, as needed.
   *
   * \param nodeSlice Node slice information.
   * \param zoneSlice Zone slice information.
   * \param n_fields The input fields.
   * \param n_newFields The output fields.
   */
  void makeFields(const SliceData &nodeSlice, const SliceData &zoneSlice, const conduit::Node &n_fields, conduit::Node &n_newFields) const
  {
    AXOM_ANNOTATE_SCOPE("makeFields");
    for(conduit::index_t i = 0; i < n_fields.number_of_children(); i++)
    {
      const conduit::Node &n_field = n_fields[i];
      const std::string association = n_field["association"].as_string();
      conduit::Node &n_newField = n_newFields[n_field.name()];
      axom::mir::utilities::blueprint::FieldSlicer<ExecSpace> fs;
      if(association == "element")
      {
        fs.execute(zoneSlice, n_field, n_newField);
      }
      else if(association == "vertex")
      {
        fs.execute(nodeSlice, n_field, n_newField);
      }
    }
  }

  /**
   * \brief Get the topology name.
   *
   * \param n_input The input mesh.
   * \param n_options The options node that may contain a "topology" string.
   *
   * \return Returns the options topology name, if present. Otherwise, it returns the first topology name.
   */
  std::string topologyName(const conduit::Node &n_input, const conduit::Node &n_options) const
  {
    std::string name;
    if(n_options.has_path("topology"))
    {
      name = n_options["topology"].as_string();
    }
    else
    {
      const conduit::Node &n_topologies = n_input.fetch_existing("topologies");
      name = n_topologies[0].name();
    }
    return name;
  }

  /**
   * \brief Return the name of the output shape type.
   *
   * \param n_topo The input topology node.
   *
   * \return The name of the output shape.
   */
  std::string outputShape(const conduit::Node &n_topo) const
  {
    std::string shape;
    if(n_topo["type"].as_string() == "unstructured")
    {
      shape = n_topo["elements/shape"].as_string();
    }
    else
    {
      if(TopologyView::dimension() == 3)
      {
        shape = "hex";
      }
      else if(TopologyView::dimension() == 2)
      {
        shape = "quad";
      }
      else
      {
        shape = "line";
      }
    }
    return shape;
  }

private:
  TopologyView m_topologyView;
  CoordsetView m_coordsetView;
};

/**
 * \brief Make a new topology and coordset by extracting certain zones from the input mesh.
 *
 * \tparam ExecSpace The execution space where the algorithm will run.
 * \tparam TopologyView The topology view type on which the algorithm will run.
 * \tparam MatsetView The matset view type on which the algorithm will run.
 */
template <typename ExecSpace, typename TopologyView, typename CoordsetView, typename MatsetView>
class ExtractZonesAndMatset : public ExtractZones<ExecSpace, TopologyView, CoordsetView>
{
public:
  using SelectedZonesView = axom::ArrayView<axom::IndexType>;

  /**
   * \brief Constructor
   *
   * \param topoView The input topology view.
   * \param coordsetView The input coordset view.
   * \param matsetView The input matset view.
   */
  ExtractZonesAndMatset(const TopologyView &topoView, const CoordsetView &coordsetView, const MatsetView &matsetView) :
    ExtractZones<ExecSpace, TopologyView, CoordsetView>(topoView, coordsetView), m_matsetView(matsetView)
  {
  }

  /**
   * \brief Select zones from the input mesh by id and output them in the output mesh.
   *
   * \param selectedZonesView A view that contains the selected zone ids.
   * \param n_input The input mesh.
   * \param n_options The input options.
   * \param[out] n_output The output mesh.
   *
   * \note The \a n_options node contains a "topology" string that is selects the
   *       name of the topology to extract.
   */
  void execute(const SelectedZonesView &selectedZonesView,
               const conduit::Node &n_input,
               const conduit::Node &n_options,
               conduit::Node &n_output) const
  {
    AXOM_ANNOTATE_SCOPE("ExtractZonesAndMatset");

    // Call base class to handle mesh/coordset/fields
    ExtractZones<ExecSpace, TopologyView, CoordsetView>::execute(selectedZonesView, n_input, n_options, n_output);

    // Make new matset.
    const std::string topoName = ExtractZones<ExecSpace, TopologyView, CoordsetView>::topologyName(n_input, n_options);
    std::string mname = matsetName(n_input, topoName);
    if(!mname.empty())
    {
      const conduit::Node &n_matset = n_input.fetch_existing("matsets/" + mname);
      conduit::Node &n_newMatset = n_output["matsets/" + mname];   
      makeMatset(selectedZonesView, n_matset, n_newMatset);
    }
  }

private:
  /**
   * \brief Return the matset for the named topology.
   *
   * \param n_input The input mesh.
   * \param topoName The name of the topology.
   *
   * \return The name of the matset for the topology or an empty string if no matset was found.
   */
  std::string matsetName(const conduit::Node &n_input, const std::string &topoName) const
  {
    std::string matset;
    if(n_input.has_child("matsets"))
    {
      const conduit::Node &n_matsets = n_input.fetch_existing("matsets");
      for(conduit::index_t i = 0; i < n_matsets.number_of_children(); i++)
      {
        const conduit::Node &n_matset = n_matsets[i];
        if(n_matset["topology"].as_string() == topoName)
        {
          matset = n_matset.name();
          break;
        }
      }
    }
    return matset;
  }

  /**
   * \brief Make a new matset that covers the selected zones.
   *
   * \param selectedZonesView A view that contains the selected zones.
   * \param n_matset The input matset.
   * \param n_newMatset A node that will contain the new matset.
   */
  void makeMatset(const SelectedZonesView &selectedZonesView, const conduit::Node &n_matset, conduit::Node &n_newMatset) const
  {
    AXOM_ANNOTATE_SCOPE("makeMatset");
    MatsetSlicer<ExecSpace, MatsetView> ms;
    ms.execute(m_matsetView, selectedZonesView, n_matset, n_newMatset);
  }

  MatsetView   m_matsetView;
};



} // end namespace blueprint
} // end namespace utilities
} // end namespace mir
} // end namespace axom

#endif
