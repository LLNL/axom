// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for internals.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_ELVIRA_ALGORITHM_HPP_
#define AXOM_MIR_ELVIRA_ALGORITHM_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/mir.hpp"
#include "axom/slic.hpp"

// Include these directly for now.
#include "axom/mir/MIRAlgorithm.hpp"
#include "axom/mir/utilities/ExtractZones.hpp"
#include "axom/mir/utilities/MergeMeshes.hpp"
#include "axom/mir/utilities/NodeToZoneRelationBuilder.hpp"
#include "axom/mir/utilities/RecenterField.hpp"
#include "axom/mir/utilities/ZoneListBuilder.hpp"
#include "axom/mir/views/dispatch_coordset.hpp"
#include "axom/mir/views/MaterialView.hpp"

#include <conduit/conduit.hpp>

// RAJA
#if defined(AXOM_USE_RAJA)
  #include "RAJA/RAJA.hpp"
#endif

#include <algorithm>
#include <string>

// Uncomment to save inputs and outputs.
// #define AXOM_ELVIRA_DEBUG

#if defined(AXOM_ELVIRA_DEBUG)
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

/*!
 * \accelerated
 * \brief Implements Elvira algorithm on the GPU using Blueprint inputs/outputs.
 */
template <typename ExecSpace, typename TopologyView, typename CoordsetView, typename MatsetView>
class ElviraAlgorithm : public axom::mir::MIRAlgorithm
{
  using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;

public:
  using ConnectivityType = typename TopologyView::ConnectivityType;

  /*!
   * \brief Constructor
   *
   * \param topoView The topology view to use for the input data.
   * \param coordsetView The coordset view to use for the input data.
   * \param matsetView The matset view to use for the input data.
   */
  EquiZAlgorithm(const TopologyView &topoView,
                 const CoordsetView &coordsetView,
                 const MatsetView &matsetView)
    : axom::mir::MIRAlgorithm()
    , m_topologyView(topoView)
    , m_coordsetView(coordsetView)
    , m_matsetView(matsetView)
  { }

  /// Destructor
  virtual ~EquiZAlgorithm() = default;

// The following members are protected (unless using CUDA)
#if !defined(__CUDACC__)
protected:
#endif

#if defined(AXOM_EQUIZ_DEBUG)
  void printNode(const conduit::Node &n) const
  {
    conduit::Node options;
    options["num_children_threshold"] = 10000;
    options["num_elements_threshold"] = 10000;
    n.to_summary_string_stream(std::cout, options);
  }
#endif

  /*!
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
    AXOM_ANNOTATE_SCOPE("ElviraAlgorithm");

    // Copy the options.
    conduit::Node n_options_copy;
    bputils::copy<ExecSpace>(n_options_copy, n_options);
    n_options_copy["topology"] = n_topo.name();

#if defined(AXOM_ELVIRA_DEBUG)
    // Save the MIR input.
    conduit::Node n_tmpInput;
    n_tmpInput[n_topo.path()].set_external(n_topo);
    n_tmpInput[n_coordset.path()].set_external(n_coordset);
    n_tmpInput[n_fields.path()].set_external(n_fields);
    n_tmpInput[n_matset.path()].set_external(n_matset);
    conduit::relay::io::blueprint::save_mesh(n_tmpInput,
                                             "debug_elvira_input",
                                             "hdf5");
#endif
  }

  /*!
   * \brief Adds original ids field to supplied fields node.
   *
   * \param n_field The new field node.
   * \param topoName The topology name for the field.
   * \param association The field association.
   * \param nvalues The number of nodes in the field.
   *
   * \note This field is added to the mesh before feeding it through MIR so we will have an idea
   *       of which nodes are original nodes in the output. Blended nodes may not have good values
   *       but there is a mask field that can identify those nodes.
   */
  void addOriginal(conduit::Node &n_field,
                   const std::string &topoName,
                   const std::string &association,
                   axom::IndexType nvalues) const
  {
    AXOM_ANNOTATE_SCOPE("addOriginal");
    namespace bputils = axom::mir::utilities::blueprint;
    bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;

    // Add a new field for the original ids.
    n_field["topology"] = topoName;
    n_field["association"] = association;
    n_field["values"].set_allocator(c2a.getConduitAllocatorID());
    n_field["values"].set(
      conduit::DataType(bputils::cpp2conduit<ConnectivityType>::id, nvalues));
    auto view = bputils::make_array_view<ConnectivityType>(n_field["values"]);
    axom::for_all<ExecSpace>(
      nvalues,
      AXOM_LAMBDA(axom::IndexType index) {
        view[index] = static_cast<ConnectivityType>(index);
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
