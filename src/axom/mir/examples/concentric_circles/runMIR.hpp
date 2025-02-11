// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_EXAMPLES_CONCENTRIC_CIRCLES_RUNMIR_HPP
#define AXOM_MIR_EXAMPLES_CONCENTRIC_CIRCLES_RUNMIR_HPP
#include "axom/config.hpp"
#include "axom/core.hpp"  // for axom macros
#include "axom/slic.hpp"
#include "axom/mir.hpp"  // for Mir classes & functions

template <typename ExecSpace>
int runMIR(const conduit::Node &hostMesh,
           const conduit::Node &options,
           conduit::Node &hostResult)
{
  AXOM_ANNOTATE_BEGIN("runMIR");

  namespace bputils = axom::mir::utilities::blueprint;
  using namespace axom::mir::views;
  SLIC_INFO(axom::fmt::format("Using policy {}",
                              axom::execution_space<ExecSpace>::name()));

  // Check materials.
  constexpr int MAXMATERIALS = 20;
  auto materialInfo = materials(hostMesh["matsets/mat"]);
  if(materialInfo.size() >= MAXMATERIALS)
  {
    SLIC_WARNING(
      axom::fmt::format("To use more than {} materials, recompile with "
                        "larger MAXMATERIALS value.",
                        MAXMATERIALS));
    return -4;
  }

  conduit::Node deviceMesh;
  {
    AXOM_ANNOTATE_SCOPE("host->device");
    bputils::copy<ExecSpace>(deviceMesh, hostMesh);
  }

  // _equiz_mir_start
  // Make views (we know beforehand which types to make)
  using CoordsetView = ExplicitCoordsetView<float, 2>;
  CoordsetView coordsetView(
    bputils::make_array_view<float>(deviceMesh["coordsets/coords/values/x"]),
    bputils::make_array_view<float>(deviceMesh["coordsets/coords/values/y"]));

  using TopoView = UnstructuredTopologySingleShapeView<QuadShape<int>>;
  TopoView topoView(bputils::make_array_view<int>(
    deviceMesh["topologies/mesh/elements/connectivity"]));

  using MatsetView = UnibufferMaterialView<int, float, MAXMATERIALS>;
  MatsetView matsetView;
  matsetView.set(
    bputils::make_array_view<int>(deviceMesh["matsets/mat/material_ids"]),
    bputils::make_array_view<float>(deviceMesh["matsets/mat/volume_fractions"]),
    bputils::make_array_view<int>(deviceMesh["matsets/mat/sizes"]),
    bputils::make_array_view<int>(deviceMesh["matsets/mat/offsets"]),
    bputils::make_array_view<int>(deviceMesh["matsets/mat/indices"]));

  using MIR =
    axom::mir::EquiZAlgorithm<ExecSpace, TopoView, CoordsetView, MatsetView>;
  MIR m(topoView, coordsetView, matsetView);
  conduit::Node deviceResult;
  m.execute(deviceMesh, options, deviceResult);
  // _equiz_mir_end

  {
    AXOM_ANNOTATE_SCOPE("device->host");
    bputils::copy<axom::SEQ_EXEC>(hostResult, deviceResult);
  }

  return 0;
}

// Prototypes.
int runMIR_seq(const conduit::Node &mesh,
               const conduit::Node &options,
               conduit::Node &result);
int runMIR_omp(const conduit::Node &mesh,
               const conduit::Node &options,
               conduit::Node &result);
int runMIR_cuda(const conduit::Node &mesh,
                const conduit::Node &options,
                conduit::Node &result);
int runMIR_hip(const conduit::Node &mesh,
               const conduit::Node &options,
               conduit::Node &result);

#endif
