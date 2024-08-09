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

namespace axom
{
namespace mir
{
#if 0
/**
 * \brief Populate a new field
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

// TODO we can 

/**
 * \accelerated
 * \brief Implements Meredith's Equi-Z algorithm on the GPU using Blueprint inputs/outputs.
 */
template <typename ExecSpace, typename TopologyView, typename CoordsetView, typename MatsetView>
class EquiZAlgorithm : public axom::mir::MIRAlgorithm
{
public:

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
    conduit::Node n_options_copy(n_options);

    // Iterate over the materials 
    const auto matInfo = axom::mir::views::materials(n_matset);
    conduit::Node n_InputTopo, n_InputCoordset, n_InputFields, n_InputMatset;
    for(size_t i = 0; i < matInfo.size(); i++)
    {
      if(i == 0)
      {
        // The first time through, we can use the supplied views.
        // clangformat-off
        iteration(
          m_topologyView, m_coordsetView, m_matsetView,
          matInfo[i],
          n_topo, n_coordset, n_fields, n_matset,
          n_options_copy,
          n_newTopo, n_newCoordset, n_newFields, n_newMatset);
        // clangformat-on

        // In these iterations, we do not want to pass selectedZones through
        // since they are only valid on the current input topology.
        n_options_copy.remove("selectedZones");
      }
      else
      {
        // Move the outputs of the last iteration to the inputs of this iteration.
        n_InputTopo.move(n_newTopo);
        n_InputCoordset.move(n_newCoordset);
        n_InputFields.move(n_newFields);
        n_InputMatset.move(n_newFields);

        const auto shape = n_newTopo.fetch_existing("elements/shape").as_string();
        if(shape == "mixed")
        {
          // The data are now an unstructured view, probably a mixed shape view.
          // Dispatch to an appropriate topo view
          // clangformat-off
          views::dispatch_explicit_coordset(n_InputCoordset, [&](auto coordsetView)
          {
            using ICoordSetView = decltype(coordsetView);
            views::dispatch_unstructured_mixed_topology(n_InputTopo, [&](auto topologyView)
            {
              using ITopologyView = decltype(topologyView);

              // Assume we made this type out of the first iteration.
              axom::mir::views::UnibufferMaterialView<int, float, 10> matsetView;
              matsetView.set(
                axom::mir::utilities::blueprint::make_array_view<int>(n_matset["material_ids"]),
                axom::mir::utilities::blueprint::make_array_view<float>(n_matset["volume_fractions"]),
                axom::mir::utilities::blueprint::make_array_view<int>(n_matset["sizes"]),
                axom::mir::utilities::blueprint::make_array_view<int>(n_matset["offsets"]),
                axom::mir::utilities::blueprint::make_array_view<int>(n_matset["indices"]));

//              this->template iteration<ITopologyView, ICoordsetView, MatsetView>(
// See if the compiler will infer the right args.
              iteration(
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
   */
  template <typename ITopologyView, typename ICoordsetView, typename IMatsetView>
  void iteration(const ITopoView &topoView,
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
    utilities::blueprint::ConduitAllocateThroughAxom<ExecSpace> c2a;

    const std::string zoneCenteredField("__equiz__zoneMatField");
    const std::string nodeCenteredField("__equiz__nodeMatField");

    const auto nzones = topoView.numberOfZones();

    // Make a node to zone relation.
    conduit::Node relation;
    axom::mir::NodeToZoneRelationBuilder<ExecSpace> builder;
    builder.execute(n_topo, relation); // <----------------------- Should this algorithm take a topo view?

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
    MatsetToField<ExecSpace, IMatsetView> m2f;
    auto zoneMatFieldView = make_array_view<FloatType>(zoneValues);
    m2f.execute(matsetView, currentMat.number, zoneMatFieldView);

    // Recenter the field to nodal using the node to zone relation.
    conduit::Node &nodeMatField = tmpFields[nodeCenteredField];
    axom::mir::RecenterField<ExecSpace> recenter;
    recenter.execute(zoneMatField, relation, nodeMatField);

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
    axom::mir::clipping::ClipField<ExecSpace, ITopologyView, ICoordsetView> clipper(topoView, coordsetView);
    clipper.execute(n_topo, n_coordset, tmpFields, options, n_newTopo, n_newCoordset, n_newFields);
    // Q: Would it be better here to just use ClipFieldFilterDevice? We don't need all that flexibility but it might be better for linking since it would have been created already for ClipFieldFilter.

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
  }

private:
  TopologyView m_topologyView;
  CoordsetView m_coordsetView;
  MatsetView m_matsetView;
};
//------------------------------------------------------------------------------
#endif
}  // end namespace mir
}  // end namespace axom

#endif
