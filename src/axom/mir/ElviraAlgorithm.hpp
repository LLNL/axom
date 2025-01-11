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

/*!
 * \accelerated
 * \brief Implements Elvira algorithm on the GPU using Blueprint inputs/outputs.
 */
template <typename ExecSpace, typename IndexPolicy, typename CoordsetView, typename MatsetView>
class ElviraAlgorithm2D : public axom::mir::MIRAlgorithm
{
  using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;

public:
  using TopologyView = axom::mir::views::StructuredTopologyView<IndexPolicy>;
  using ConnectivityType = typename TopologyView::ConnectivityType;

  // Make sure we only instantiate this algorithm for 2D meshes.
  static_assert(IndexPolicy::dimension() == 2, "IndexPolicy must be 2D");

  /*!
   * \brief Constructor
   *
   * \param topoView The topology view to use for the input data.
   * \param coordsetView The coordset view to use for the input data.
   * \param matsetView The matset view to use for the input data.
   */
  ElviraAlgorithm2D(const TopologyView &topoView,
                    const CoordsetView &coordsetView,
                    const MatsetView &matsetView)
    : axom::mir::MIRAlgorithm()
    , m_topologyView(topoView)
    , m_coordsetView(coordsetView)
    , m_matsetView(matsetView)
  { }

  /// Destructor
  virtual ~ElviraAlgorithm2D() = default;

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
    AXOM_ANNOTATE_SCOPE("ElviraAlgorithm2D");

    // Copy the options to make sure they are in the right memory space.
    conduit::Node n_options_copy;
    bputils::copy<ExecSpace>(n_options_copy, n_options);
    n_options_copy["topology"] = n_topo.name();


  }

private:
  TopologyView m_topologyView;
  CoordsetView m_coordsetView;
  MatsetView m_matsetView;
};

}  // end namespace mir
}  // end namespace axom

#endif
