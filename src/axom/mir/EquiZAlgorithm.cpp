// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/mir/EquiZAlgorithm.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/mir/views/StructuredTopologyView.hpp"
#include "axom/mir/views/dispatch_coordset.hpp"
#include "axom/mir/views/dispatch_topology.hpp"

#include <conduit/conduit_blueprint.hpp>

// RAJA
#if defined(AXOM_USE_RAJA)
  #include "RAJA/RAJA.hpp"
#endif

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

void EquiZAlgorithm::execute(const conduit::Node &topo,
                                 const conduit::Node &coordset,
                                 const conduit::Node &options,
                                 conduit::Node &new_topo,
                                 conduit::Node &new_coordset)
{
#if defined (AXOM_USE_RAJA) && defined (AXOM_USE_UMPIRE)
  switch(m_execPolicy)
  {
  #if defined(AXOM_USE_OPENMP)
  case RuntimePolicy::omp:
    executeImpl<omp_exec>(topo, coordset, options, new_topo, new_coordset);
    break;
  #endif
  #if defined(AXOM_USE_CUDA)
  case RuntimePolicy::cuda:
    executeImpl<cuda_exec>(topo, coordset, options, new_topo, new_coordset);
    break;
  #endif
  #if defined(AXOM_USE_HIP)
  case RuntimePolicy::hip:
    executeImpl<hip_exec>(topo, coordset, options, new_topo, new_coordset);
    break;
  #endif
  default:
    // Falls through
  case RuntimePolicy::seq:
    executeImpl<seq_exec>(topo, coordset, options, new_topo, new_coordset);
    break;
  }
#endif
}

template <typename ExecSpace>
void EquiZAlgorithm::executeImpl(const conduit::Node &topo,
                                 const conduit::Node &coordset,
                                 const conduit::Node &options,
                                 conduit::Node &new_topo,
                                 conduit::Node &new_coordset)
{
  if(options.has_path("zones"))
  {
    const conduit::Node &zones = options.fetch_existing("zones");
    const auto nzones = zones.dtype().number_of_elements();

/// NOTE: since each inner dispatch could be a lot of code, should I just make a zones array for the case where zones is not provided?

    // Operate on a list of zones.
    views::IndexNode_To_ArrayView(options["zones"], [&](auto zonesView)
    {
      views::dispatch_coordset(coordset, [&](auto &coordsetView)
      {
        views::dispatch_topology(topo, coordset, [&](const std::string &shape, auto &topoView)
        {
          topoView. template for_selected_zones<ExecSpace>(zonesView, AXOM_LAMBDA(auto zoneIndex, const auto &zone)
          {          
  
          });
        });
      });
    });
  }
  else
  {
    // Operate on all zones.
    views::dispatch_coordset(coordset, [&](auto &coordsetView)
    {
      views::dispatch_topology(topo, coordset, [&](const std::string &shape, auto &topoView)
      {
        topoView. template for_all_zones<ExecSpace>(AXOM_LAMBDA(auto zoneIndex, const auto &zone)
        {          

        });
      });
    });

  }
}

} // namespace mir
} // namespace axom
