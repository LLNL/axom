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

} // end namespace blueprint
} // end namespace utilities
} // end namespace mir
} // end namespace axom

#endif
