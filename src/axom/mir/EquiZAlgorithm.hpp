// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_EQUIZ_ALGORITHM_HPP_
#define AXOM_MIR_EQUIZ_ALGORITHM_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/mir.hpp"

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
  
  void execute(const MatsetView &matsetView, MaterialIndex mat, axom::ArrayView<FloatType> &vfValues)
  {
    const auto nzones = vfValues.size();
    axom::for_all<ExecSpace>(nzones, AXOM_LAMBDA(auto zoneIndex)
    {
      FloatType vf{};
      matsetView.zoneContainsMaterial(zoneIndex, mat, vf);
      vfValues[zoneIndex] = vf;
    });
  }
};
#if 0
/**
 * \accelerated
 * \brief Implements Meredith's Equi-Z algorithm on the GPU using Blueprint inputs/outputs.
 */
template <typename ExecSpace, typename TopologyView, typename CoordsetView, typename MatsetView>
class EquiZAlgorithm : public axom::mir::MIRAlgorithm
{
public:
#if 0
  EquiZAlgorithm(const TopologyView &topoView, const CoordsetView &coordsetView, const MatsetView &matsetView) :
    m_topologyView(topoView), m_coordsetView(coordsetView), m_matsetView(matsetView)
  {
  }

  virtual ~EquiZAlgorithm() = default;

protected:
  virtual void executeDomain(const conduit::Node &n_topo,
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

    // Copy the options.
    conduit::Node n_options_copy;
    bputils::copy<ExecSpace>(n_options_copy, n_options);

    // Iterate over the materials 
    const auto matInfo = axom::mir::views::materials(n_matset);
    conduit::Node n_InputTopo, n_InputCoordset, n_InputFields, n_InputMatset;
    for(size_t i = 0; i < matInfo.size(); i++)
    {
      if(i == 0)
      {
        // The first time through, we can use the supplied views.
        // clangformat-off
        iteration<TopologyView, CoordsetView, MatsetView>(
          m_topologyView, m_coordsetView, m_matsetView,
          matInfo[i],
          n_topo, n_coordset, n_fields, n_matset,
          n_options_copy,
          n_newTopo, n_newCoordset, n_newFields, n_newMatset);
        // clangformat-on

        // In later iterations, we do not want to pass selectedZones through
        // since they are only valid on the current input topology. Also, if they
        // were passed then the new topology only has those selected zones.
        n_options_copy.remove("selectedZones");
      }
      else
      {
        // Move the outputs of the last iteration to the inputs of this iteration.
        n_InputTopo.move(n_newTopo);
        n_InputCoordset.move(n_newCoordset);
        n_InputFields.move(n_newFields);
        n_InputMatset.move(n_newMatset);

        const auto shape = n_newTopo.fetch_existing("elements/shape").as_string();
        if(shape == "mixed")
        {
          // The data are now an unstructured view, probably a mixed shape view.
          // Dispatch to an appropriate topo view
          // clangformat-off
          views::dispatch_explicit_coordset(n_InputCoordset, [&](auto coordsetView)
          {
            using ICoordsetView = decltype(coordsetView);
            views::dispatch_unstructured_mixed_topology(n_InputTopo, [&](auto topologyView)
            {
              using ITopologyView = decltype(topologyView);

              // Assume we made this type out of the first iteration.
              using IMatsetView = axom::mir::views::UnibufferMaterialView<int, float, 10>;
              IMatsetView matsetView;
              matsetView.set(
                bputils::make_array_view<int>(n_matset["material_ids"]),
                bputils::make_array_view<float>(n_matset["volume_fractions"]),
                bputils::make_array_view<int>(n_matset["sizes"]),
                bputils::make_array_view<int>(n_matset["offsets"]),
                bputils::make_array_view<int>(n_matset["indices"]));
             
//              this->template iteration<ITopologyView, ICoordsetView, MatsetView>(
// See if the compiler will infer the right args.
              iteration<ITopologyView, ICoordsetView, IMatsetView>(
                topologyView, coordsetView, matsetView,
                matInfo[i],
                n_topo, n_coordset, n_fields, n_matset,
                n_options_copy,
                n_newTopo, n_newCoordset, n_newFields, n_newMatset);
            });
          });
          // clangformat-on
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
    using FloatType = typename MatsetView::FloatType;
    constexpr auto floatTypeID = bputils::cpp2conduit<FloatType>::id;

    // Get the ID of a Conduit allocator that will allocate through Axom with device allocator allocatorID.
    bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;

    const std::string zoneCenteredField("__equiz__zoneMatField");
    const std::string nodeCenteredField("__equiz__nodeMatField");
    const std::string colorField("__equiz__colors");

    const auto nzones = topoView.numberOfZones();

    // Make a node to zone relation.
    conduit::Node relation;
    {
      bputils::NodeToZoneRelationBuilder<ExecSpace> rb;
      rb.execute(n_topo, relation); // <----------------------- Should this algorithm take a topo view?
    }

    // Make a shallow copy of the fields that we can modify.
    conduit::Node tmpFields;
    for(conduit::index_t i = 0; i < n_fields.number_of_children(); i++)
    {
      tmpFields[n_fields[i].name()].set_external(n_fields[i]);
    }

    // Make a zonal field for the current material's volume fractions.
    conduit::Node &zoneMatField = tmpFields[zoneCenteredField];
    zoneMatField["topology"] = n_topo.name();
    zoneMatField["association"] = "element";
    conduit::Node &zoneValues = zoneMatField["values"];
    zoneValues.set_allocator(c2a.getConduitAllocatorID());
    zoneValues.set(conduit::DataType(floatTypeID, nzones));

    // Populate the zonal volume fraction field from the material view.
    auto zoneMatFieldView = bputils::make_array_view<FloatType>(zoneValues);
    {
      MatsetToField<ExecSpace, IMatsetView> m2f;
      m2f.execute(matsetView, currentMat.number, zoneMatFieldView);
    }

    // Recenter the field to nodal using the node to zone relation.
    conduit::Node &nodeMatField = tmpFields[nodeCenteredField];
    {
      bputils::RecenterField<ExecSpace> recenter;
      recenter.execute(zoneMatField, relation, nodeMatField);
    }

/**
 NOTE - I am a little worried that I've not yet appreciated how the intersections are found
        in the EquiZ algorithm and that I'll need to add a new template parameter to ClipField
        that lets that behavior be customizable so I can override it here.
 */

    // Now, clip the topology using the nodal field.
    conduit::Node options;
    options["inside"] = 1;
    options["outside"] = 1;
    options["clipField"] = nodeCenteredField;
    options["clipValue"] = 0.5;
    options["colorField"] = colorField;
    if(n_options.has_child("selectedZones"))
    {
      // Pass selectedZones along in the clip options, if present.
      options["selectedZones"].set_external(n_options["selectedZones"]);
    }
    {
      axom::mir::clipping::ClipField<ExecSpace, ITopologyView, ICoordsetView> clipper(topoView, coordsetView);
      clipper.execute(n_topo, n_coordset, tmpFields, options, n_newTopo, n_newCoordset, n_newFields);
      // Q: Would it be better here to just use ClipFieldFilterDevice? We don't need all that flexibility but it might be better for linking since it would have been created already for ClipFieldFilter.
    }

  #if 1
    conduit::Node mesh;
    mesh[n_newTopo.path()].set_external(n_newTopo);
    mesh[n_newCoordset.path()].set_external(n_newCoordset);
    mesh[n_newFields.path()].set_external(n_newFields);
    mesh[n_newMatset.path()].set_external(n_newMatset);
    std::stringstream ss;
    ss << "debug_equiz_" << currentMat.number;
    conduit::relay::io::blueprint::save_mesh(mesh, ss.str(), "hdf5");
  #endif

    n_newFields.remove(zoneCenteredField);
    n_newFields.remove(nodeCenteredField);

    // Make a new matset.
    auto colorView = axom::mir::utilities::blueprint::make_array_view<int>(n_newFields[colorField +"/values"]);
    auto origElemView = axom::mir::utilities::blueprint::make_array_view<int>(n_newFields["originalElements/values"]);
    makeNewMatset<IMatsetView>(matsetView, currentMat, colorView, origElemView, n_matset, n_newMatset);

    n_newFields.remove("originalElements");
  }

  template <typename IMatsetView>
  void makeNewMatset(// by value
                     const IMatsetView currentMatsetView,
                     const axom::mir::views::Material currentMat,
                     const axom::ArrayView<int> colorView, 
                     const axom::ArrayView<int> originalElementsView,
                     // by reference
                     const conduit::Node &n_matset,
                     conduit::Node &n_newMatset)
  {
    namespace bputils = axom::mir::utilities::blueprint;
    using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;
    using IntType = typename IMatsetView::IndexType;
    using FloatType = typename IMatsetView::FloatType;
    const int intTypeID = bputils::cpp2conduit<IntType>::id;
    const int floatTypeID = bputils::cpp2conduit<FloatType>::id;

    // Copy the material map to the new matset.
    if(n_matset.has_child("material_map"))
    {
      n_newMatset.set(n_matset.fetch_existing("material_map"));
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
    auto sizesView = axom::mir::utilities::blueprint::make_array_view<IntType>(n_sizes);
    axom::for_all<ExecSpace>(nzones, AXOM_LAMBDA(auto zoneIndex)
    {
      axom::IndexType nmatsThisZone = 0;
      if(colorView[zoneIndex] == COLOR0)
      {
        nmatsThisZone = 1;
      }
      else
      {
        // These were the unselected parts of the zone after clipping.
        // When we make the new material, this set of zones remove the current material.
        nmatsThisZone = currentMatsetView.numberOfMaterials(zoneIndex) - 1;
      }
      sizesView[zoneIndex] = nmatsThisZone;
      sum += nmatsThisZone;
    });

    // Make offsets.
    auto offsetsView = axom::mir::utilities::blueprint::make_array_view<IntType>(n_offsets);
    axom::exclusive_scan<ExecSpace>(sizesView, offsetsView);

    // Make new matset data.
    const auto arraySize = sum.get();
    conduit::Node &n_material_ids = n_newMatset["material_ids"];
    conduit::Node &n_volume_fractions = n_newMatset["volume_fractions"];
    conduit::Node &n_indices = n_newMatset["indices"];
    n_material_ids.set_allocator(c2a.getConduitAllocatorID());
    n_volume_fractions.set_allocator(c2a.getConduitAllocatorID());
    n_indices.set_allocator(c2a.getConduitAllocatorID());
    n_material_ids.set(conduit::DataType(intTypeID, arraySize));
    n_volume_fractions.set(conduit::DataType(floatTypeID, arraySize));
    n_indices.set(conduit::DataType(intTypeID, arraySize));

    auto matidsView = axom::mir::utilities::blueprint::make_array_view<IntType>(n_material_ids);
    auto vfView = axom::mir::utilities::blueprint::make_array_view<FloatType>(n_volume_fractions);
    auto indicesView = axom::mir::utilities::blueprint::make_array_view<IntType>(n_material_ids);

    // Fill in the indices used.
    axom::for_all<ExecSpace>(arraySize, AXOM_LAMBDA(auto index)
    {
      indicesView[index] = index;
    });

    // Fill in the material data.
    const int currentMatNumber = currentMat.number;
    axom::for_all<ExecSpace>(nzones, AXOM_LAMBDA(auto zoneIndex)
    {
      auto offset = offsetsView[zoneIndex];
      if(colorView[zoneIndex] == COLOR0)
      {
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

        // Total up the materials, excluding the current material.
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
            offset++;
          }
        }
      }
    });
  }

private:
  TopologyView m_topologyView;
  CoordsetView m_coordsetView;
  MatsetView m_matsetView;
#endif
};
#endif
}  // end namespace mir
}  // end namespace axom

#endif
