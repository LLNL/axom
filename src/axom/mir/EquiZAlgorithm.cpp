// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/mir/EquiZAlgorithm.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/mir/views/StructuredTopologyView.hpp"
#include "axom/mir/views/dispatch_coordset.hpp"
#include "axom/mir/views/dispatch_topology.hpp"
#include "axom/mir/views/dispatch_material.hpp"
#include "axom/mir/views/dispatch_utilities.hpp"
#include "axom/mir/clipping/ClipTableManager.hpp"
#include "axom/mir/NodeToZoneRelationBuilder.hpp"

#include <conduit/conduit_blueprint.hpp>

// RAJA
#if defined(AXOM_USE_RAJA)
  #include "RAJA/RAJA.hpp"
#endif

#if 0
// clang-format off
#if defined (AXOM_USE_RAJA) && defined (AXOM_USE_UMPIRE)
  using seq_exec = axom::SEQ_EXEC;

  #if defined(AXOM_USE_OPENMP)
    using omp_exec = axom::OMP_EXEC;
  #else
    using omp_exec = seq_exec;
  #endif

  #if defined(AXOM_USE_CUDA)
    constexpr int CUDA_BLOCK_SIZE = 256;
    using cuda_exec = axom::CUDA_EXEC<CUDA_BLOCK_SIZE>;
  #else
    using cuda_exec = seq_exec;
  #endif

  #if defined(AXOM_USE_HIP)
    constexpr int HIP_BLOCK_SIZE = 64;
    using hip_exec = axom::HIP_EXEC<HIP_BLOCK_SIZE>;
  #else
    using hip_exec = seq_exec;
  #endif
#endif
// clang-format on
#endif

namespace axom
{
namespace mir
{
void EquiZAlgorithm::execute(const conduit::Node &topo,
                             const conduit::Node &coordset,
                             const conduit::Node &matset,
                             const conduit::Node &options,
                             conduit::Node &new_topo,
                             conduit::Node &new_coordset,
                             conduit::Node &new_matset)
{
#if 0
  #if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  switch(m_execPolicy)
  {
    #if defined(AXOM_USE_OPENMP)
  case RuntimePolicy::omp:
    executeImpl<omp_exec>(topo, coordset, matset, options, new_topo, new_coordset, new_matset);
    break;
    #endif
    #if defined(AXOM_USE_CUDA)
  case RuntimePolicy::cuda:
    executeImpl<cuda_exec>(topo, coordset, matset, options, new_topo, new_coordset, new_matset);
    break;
    #endif
    #if defined(AXOM_USE_HIP)
  case RuntimePolicy::hip:
    executeImpl<hip_exec>(topo, coordset, matset, options, new_topo, new_coordset, new_matset);
    break;
    #endif
  default:
    // Falls through
  case RuntimePolicy::seq:
    executeImpl<seq_exec>(topo, coordset, matset, options, new_topo, new_coordset, new_matset);
    break;
  }
  #endif
#endif
}

template <typename ExecSpace>
void EquiZAlgorithm::executeImpl(const conduit::Node &topo,
                                 const conduit::Node &coordset,
                                 const conduit::Node &matset,
                                 const conduit::Node &options,
                                 conduit::Node &new_topo,
                                 conduit::Node &new_coordset,
                                 conduit::Node &new_matset)
{
  // TODO: migrate data for topo, coordset to appropriate memory space if needed.
#if 0
  views::dispatch_coordset(coordset, [&](auto &coordsetView)
  {
    views::dispatch_topology<views::select_dimensions(2,3)>(topo, coordset, [&](const std::string &shape, auto &topoView)
    {
      views::dispatch_material(matset, [&](auto &matsetView)
      {
        EquiZAlgorithm<ExecSpace, decltype(topoView), decltype(coordsetView), decltype(matsetView)> mir;
        mir.execute(topoView, coordsetView, matsetView, options);
      });
    });
  });
#endif
#if 0
  if(options.has_path("zones"))
  {
    const conduit::Node &n_zones = options.fetch_existing("zones");

/// NOTE: since each inner dispatch could be a lot of code, should I just make a zones array for the case where zones is not provided?

    // Operate on a list of zones.
    views::dispatch_material(matset, [&](auto &matsetView)
    {
      views::IndexNode_to_ArrayView(n_zones, [&](auto zonesView)
      {
        views::dispatch_topology<views::select_dimensions(2,3)>(topo, coordset, [&](const std::string &shape, auto &topoView)
        {
          views::dispatch_coordset(coordset, [&](auto &coordsetView)
          {
            // Create the clipping tables for the topo dimension.
            axom::mir::clipping::ClipTableManager<ExecSpace> clipManager;
            clipManager.load(topoView.dimension());

            const auto matinfo = views::materials(matset);
            for(const auto &mat : matinfo)
            {
              const auto matID = mat.number;

              axom::mir::utilities::NodeToZoneRelationBuilder<ExecSpace> nz;
              conduit::Node rel;
              nz.execute(topo, rel);

              views::IndexNode_to_ArrayView_same(rel["zones"], rel["sizes"], rel["offsets"], [&](auto relZonesView, auto relSizesView, auto relOffsetsView)
              {
                // We need to get views for the various shape types.
                using TableView = typename axom::mir::clipping::ClipTableManager<ExecSpace>::Table::View; 
                axom::StackArray<TableView, ST_MAX> tables;
                tables[ST_TET] = clipManager[ST_TET].view();

                topoView. template for_selected_zones<ExecSpace>(zonesView, AXOM_LAMBDA(auto zoneIndex, const auto &zone)
                {          
                
                });
              });
            }
          }); // dispatch_matset
        }); // dispatch_coordset
      }); // dispatch_topology
    });
  }
  else
  {
    // Operate on all zones.
    views::dispatch_coordset(coordset, [&](auto &coordsetView)
    {
      views::dispatch_topology<views::select_dimensions(2,3)>(topo, coordset, [&](const std::string &shape, auto &topoView)
      {
        views::dispatch_material(matset, [&](auto &matsetView)
        {
          topoView. template for_all_zones<ExecSpace>(AXOM_LAMBDA(auto zoneIndex, const auto &zone)
          {          

          });
        });
      });
    });

  }
#endif
}

}  // namespace mir
}  // namespace axom
