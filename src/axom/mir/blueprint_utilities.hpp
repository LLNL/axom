// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_BLUEPRINT_UTILITIES_HPP_
#define AXOM_MIR_BLUEPRINT_UTILITIES_HPP_

#include "axom/core/execution/execution_space.hpp"
#include "axom/core/Array.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/mir/views/dispatch_structured_topologies.hpp"

#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint.hpp>

namespace axom
{
namespace mir
{
namespace utilities
{
namespace blueprint
{

// TODO: Add in a routine to migrate a Conduit node to a new memory space.
// copy(const conduit::Node &src, conduit::Node &dest);

/**
 * \accelerated
 * \brief Make an unstructured representation of a structured topology.
 *
 * \tparam ExecSpace The execution space where the work will be done.
 *
 * \param
 *
 * \note There are blueprint methods for this sort of thing but this one is accelerated.
 */
template <typename ExecSpace>
void
to_unstructured(const conduit::Node &topo, const conduit::Node &coordset, const std::string &topoName, conduit::Node &mesh)
{
  const std::string type = topo.fetch_existing("type").as_string();
  const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();

  mesh["coordsets"][coordset.name()].set_external(coordset);
  conduit::Node &newtopo = mesh["topologies"][topoName];

  if(type == "unstructured")
  {
    newtopo.set_external(topo);
  }
  else
  {
    newtopo["type"] = "unstructured";
    conduit::Node &n_newconn = newtopo["elements/connectivity"];
    n_newconn.set_allocator(allocatorID);

    // Fill in the connectivity.
    views::dispatch_structured_topologies(topo, coordset, [&](const std::string &shape, auto &topoView)
    {
      const auto nzones = topoView.numberOfZones();
      int ptsPerZone = 2;
      if(shape == "quad")
        ptsPerZone = 4;
      else if(shape == "hex")
        ptsPerZone = 8;

      newtopo["elements/shape"] = shape;

      const auto connSize = nzones * ptsPerZone;
      n_newconn.set(conduit::DataType::index_t(connSize));
      axom::ArrayView<conduit::index_t> conn(reinterpret_cast<conduit::index_t *>(n_newconn.data_ptr()), connSize);
      auto conn_view = conn.view();
      topoView. template for_all_zones<ExecSpace>(AXOM_LAMBDA(auto zoneIndex, const auto &zone)
      {
        const auto start = zoneIndex * ptsPerZone;
        for(int i = 0; i < 4; i++)
          conn_view[start + i] = static_cast<conduit::index_t>(zone.getIds()[i]);
      });
    });
  }
}

#if 0
/**
 \brief This method slices a Conduit field using a node that contains index value 
        and stores the data in another node. The data in the nodes are assumed to
        be in the memory space suitable for the ExecSpace.

 \tparam ExecSpace The execution space.

 \param[in] field The Blueprint field to be sliced.
 \param[in] indices A Conduit node that contains an array of index values into the field.
 \param[out] output A node that contains the new sliced field.
        
 */
template <typename ExecSpace>
void
sliceField(const conduit::Node &field, const conduit::Node &indices, conduit::Node &output)
{
  // Start making an output field.
  output["topology"] = field["topology"].as_string();
  output["association"] = field["association"].as_string();
  conduit::Node &output_values = output["values"];

  const conduit::Node &field_values = field["values"];
  const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
  const axom::IndexType = indices.dtype().number_of_elements();

  // Slice the field using the indices and store values in the output field.
  if(field_values.number_of_children() > 0)
  {
    for(conduit::index_t i = 0; i < field_values.number_of_children(); i++)
    {
      const conduit::Node &componentNode = field_values[i];
      const std::string componentName(componentNode.name());

      // Allocate memory for the output view using the allocator for the ExecSpace.
      conduit::Node &destNode = output_values[componentName];
      destNode.set_allocator(allocatorID);
      destNode.set(conduit::DataType(componentNode.dtype().id(), numIndices));

      views::IndexNode_to_ArrayView(indices, [&](auto indicesView)
      {
        views::Node_to_ArrayView(componentNode, destNode, [&](auto fieldValuesView, auto destView)
        {
          axom::for_all<ExecSpace>(
            numIndices,
            AXOM_LAMBDA(axom::IndexType i) {
              destView[i] = fieldValuesView[indicesView[i]];
            });
        });
      });
    }
  }
  else
  {
    // Allocate memory for the output view using the allocator for the ExecSpace.
    output_values.set_allocator(allocatorID);
    output_values.set(conduit::DataType(field_values.dtype().id(), numIndices));

    IndexNodeToArrayView(indices, [&](auto indicesView)
    {
      NodeToArrayView(field_values, output_values, [&](auto fieldValuesView, auto destView)
      {
        axom::for_all<ExecSpace>(
          numIndices,
          AXOM_LAMBDA(axom::IndexType i) {
            destView[i] = fieldValuesView[indicesView[i]];
          });
      });
    });
  }
}

template <typename ExecSpace>
void
blendField(const conduit::Node &field, const conduit::Node &indices, const conduit::Node &weights, const conduit::Node &sizes, const conduit::Node &offsets, conduit::Node &output)
{
  // Start making an output field.
  output["topology"] = field["topology"].as_string();
  output["association"] = "element";
  conduit::Node &output_values = output["values"];

  const conduit::Node &field_values = field["values"];
  const auto n_sizes = static_cast<axom::IndexType>(sizes.dtype().number_of_elements());
  int execSpaceAllocatorID = axom::execution_space<ExecSpace>::allocatorID();

  // Slice the field using the indices and store values in the output field.
  if(field_values.number_of_children() > 0)
  {
    for(conduit::index_t i = 0; i < field_values.number_of_children(); i++)
    {
      const conduit::Node &componentNode = field_values[i];
      const std::string componentName(componentNode.name());

      // Allocate memory for the output view using the allocator for the ExecSpace.
      conduit::Node &destNode = output_values[componentName];
      destNode.set_allocator(execSpaceAllocatorID);
      destNode.set(conduit::DataType(componentNode.dtype().id(), n_sizes));

      IndexNodeToArrayView(indices, sizes, offsets, [&](auto indicesView, auto sizesView, auto offsetsView)
      {
        NodeToArrayView(componentNode, weights, destNode, [&](auto fieldValuesView, auto weightsView, auto destView)
        {
          axom::for_all<ExecSpace>(
           n_sizes,
           AXOM_LAMBDA(axom::IndexType i) {
             const auto offset = offsetsView[i];
             const auto size = sizesView[i];
             destView[i] = 0;
             for(int j = 0; j < size; j++)
             {
               const auto idx = indicesView[offset + j];
               destView[i] += fieldValuesView[idx] * weightsView[idx];
             }
           });
        });
      });
    }
  }
  else
  {
    // Allocate memory for the output view using the allocator for the ExecSpace.
    output_values.set_allocator(execSpaceAllocatorID);
    output_values.set(conduit::DataType(field_values.dtype().id(), n_sizes));

    IndexNodeToArrayView(indices, sizes, offsets, [&](auto indicesView, auto sizesView, auto offsetsView)
    {
      NodeToArrayView(field_values, weights, output_values, [&](auto fieldValuesView, auto weightsView, auto destView)
      {
        axom::for_all<ExecSpace>(
          n_sizes,
          AXOM_LAMBDA(axom::IndexType i) {
            const auto offset = offsetsView[i];
            const auto size = sizesView[i];
            destView[i] = 0;
            for(int j = 0; j < size; j++)
            {
              const auto idx = indicesView[offset + j];
              destView[i] += fieldValuesView[idx] * weightsView[idx];
            }
         });
      });
    });
  }
}

void conduit_move(const conduit::Node &src, conduit::Node &dest, int dest_allocator)
{
    if(src.number_of_children() > 0)
    {
        for(conduit::index_t i = 0; i < src.number_of_children()
        {
            conduit_move(src[i], dest[src[i].name()], allocator);
        }
    }
    else
    {
        if(src.dtype().number_of_elements() > 1)
        {
            // Allocate the node's memory in the right place.
            dest.reset();
            dest.set_allocator(allocator);
            dest.set(conduit::DataType(src.dtype().id(), src.dtype().number_of_elements()));

            // Copy the data to the destination node. Axom uses Umpire to manage that.
            if(src.is_compact())
                axom::copy(dest.data_ptr(), src.data_ptr(), src.dtype().bytes_compact());
            else
            {
                // NOTE: this assumes that src is on the host. Why would we have strided data on device?
                conduit::Node tmp;
                src.compact_to(tmp);
                axom::copy(dest.data_ptr(), tmp.data_ptr(), tmp.dtype().bytes_compact());
            }
        }
        else
        {
            // The node data fits in the node. It's on the host.
            dest.set(src);
        }
    }
}
#endif

} // end namespace blueprint
} // end namespace utilities
} // end namespace mir
} // end namespace axom

#endif
