// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


#ifndef AXOM_MIR_ALGORITHM_HPP_
#define AXOM_MIR_ALGORITHM_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/core/ArrayView.hpp"

#include <conduit/conduit.hpp>

#include <vector>
#include <string>

namespace axom
{

namespace mir
{

/**
 \brief Base class for Material Interface Reconstruction (MIR) algorithms.
 */
class MIRAlgorithm
{
public:
   MIRAlgorithm() = default;
   virtual ~MIRAlgorithm() = default;

   /**
    \brief Perform material interface reconstruction on the meshes supplied in the
           root node. Root can either be a mesh domain or a node that contains multiple
           domains.

    \param[in] root The root node that contains either a mesh or list of mesh
                    domains that contain a topology and matset to be used for MIR.
    \param[in] options A node that contains options that help govern MIR execution.

options:
  topology: main
  new_topology: mirtopo
  new_coordset: mircoords
  fields:
    - temperature
    - pressure
  zones: [0,1,6,9]
  mapping: 0

    The "topology" option specifies which topology we'll reconstruct. It must have an associated matset.
    "new_topology" is the name of the topology that will be created in the output node.
    "new_coordset" is the name of the new coordset that will be created in the output node. If it is not provided then the name of the topology's coordset will be used.  
    "fields" is the name of the fields to map to the new topology. If fields is specified but empty, no fields will be mapped. If fields is not present then all fields will be mapped.
    "zones" is a list of zone indices from the topology that need to be reconstructed. If not present then all zones will be considered.
    "mapping" indicates whether we should include an original_element_numbers field on the new topology to indicate where each new zone came from in the original topology.

    \param[out] output A node that will contain the new entities.

    */
   virtual void execute(const conduit::Node &root,
                        const conduit::Node &options,
                        conduit::Node &output);
protected:
   /**
    \brief Perform material interface reconstruction on a single domain. Derived classes
           must implement this method and any device-specific coding gets handled under it.

    \param[in] topo The Conduit node containing the topology that will be used for MIR.
    \param[in] coordset The Conduit node containing the topology's coordset.
    \param[in] options The Conduit node containing the options that help govern MIR execution.

    \param[out] new_topo A Conduit node that will contain the new topology.
    \param[out] new_coordset A Conduit node that will contain the new coordset.
    
    */
   virtual void execute(const conduit::Node &topo,
                        const conduit::Node &coordset,
                        const conduit::Node &options,
                        conduit::Node &new_topo,
                        conduit::Node &new_coordset) = 0;

   // Utility methods for derived types.
   void copyState(const conduit::Node &mesh, conduit::Node &destMesh) const;
   std::string topologyName(const conduit::Node &mesh, const conduit::Node &options) const;
   std::string newTopologyName(const conduit::Node &mesh, const conduit::Node &options) const;
   std::string newCoordsetName(const conduit::Node &mesh, const conduit::Node &options) const;
   std::vector<std::string> fieldNames(const conduit::Node &mesh, const conduit::Node &options) const;

   const conduit::Node &topology(const conduit::Node &input, const conduit::Node &options) const;
   const conduit::Node &matset(const conduit::Node &input, const conduit::Node &options) const;


   // TODO: method for mapping element field to new topo
   // TODO: method for mapping vertex field to new topo
};

#if 0
// Find a home for this stuff.


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
  output["association"] = "element";
  conduit::Node &output_values = output["values"];

  const conduit::Node &field_values = field["values"];
  const auto n_indices = static_cast<axom::IndexType>(indices.dtype().number_of_elements());
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
      destNode.set(conduit::DataType(componentNode.dtype().id(), n_indices));

      IndexNodeToArrayView(indices, [&](auto indicesView)
      {
        NodeToArrayView(componentNode, destNode, [&](auto fieldValuesView, auto destView)
        {
          axom::for_all<ExecSpace>(
            n_indices,
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
    output_values.set_allocator(execSpaceAllocatorID);
    output_values.set(conduit::DataType(field_values.dtype().id(), n_indices));

    IndexNodeToArrayView(indices, [&](auto indicesView)
    {
      NodeToArrayView(field_values, output_values, [&](auto fieldValuesView, auto destView)
      {
        axom::for_all<ExecSpace>(
          n_indices,
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

class ElviraMIRAlgorithm : public MIRAlgorithm
{
public:
  using RuntimePolicy = axom::runtime_policy::Policy;

  ElviraMIRAlgorithm() = default;
  virtual ~ElviraMIRAlgorithm() = default;

  void setExecPolicy(RuntimePolicy policy) { m_execPolicy = policy; }

protected:
   /// Implement the Elvira MIR algorithm on a single domain.
   virtual void execute(const conduit::Node &topo,
                        const conduit::Node &coordset,
                        const conduit::Node &options,
                        conduit::Node &new_topo,
                        conduit::Node &new_coordset) override;

   
   /// Implement the Elvira MIR algorithm on a single domain for a given ExecSpace.
   template <typename ExecSpace>
   void executeImpl(const conduit::Node &topo,
                    const conduit::Node &coordset,
                    const conduit::Node &options,
                    conduit::Node &new_topo,
                    conduit::Node &new_coordset);

   RuntimePolicy m_execPolicy{RuntimePolicy::seq};
};

// class EquiZMIRAlgorithm : public MIRAlgorithml


} // end namespace mir
} // end namespace axom

#endif
