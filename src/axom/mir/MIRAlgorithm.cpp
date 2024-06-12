// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/MIRAlgorithm.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/slic.hpp"

#include <conduit_blueprint_mesh.hpp>
#include <conduit_blueprint_mesh_utils.hpp>

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

namespace axom
{

namespace quest
{

void MIRAlgorithm::execute(const conduit::Node &root,
                           const conduit::Node &options,
                           conduit::Node &output)
{
  auto domains = conduit::blueprint::mesh::domains(root);
  if(domains.size() > 1)
  {
    // Handle multiple domains
    for(const auto &dom_ptr : domains)
    {
      const conduit::Node &dom = *dom_ptr;
      const std::string topoName = topologyName(dom, options);
      const std::string newTopoName = newTopologyName(dom, options);
      const std::string newCSName = newCoordsetName(dom, options);

      conduit::Node &newDomain = output.append();
      const conduit::Node &topologies = dom.fetch_existing("topologies");
      const conduit::Node &topo = topologies.fetch_existing(topoName);
      const conduit::Node *cset = conduit::blueprint::mesh::utils::find_reference_node(topo, "coordset");

      conduit::Node &newTopo = newDomain["topologies/" + newTopoName];
      conduit::Node &newCoordset = newDomain["coordsets/" + newCSName];
      copyState(dom, newDomain);
      execute(topo, *cset, options, newTopo, newCoordset);
    }
  }
  else if(domains.size() > 0)
  {
    // Handle single domain
    const conduit::Node &dom = *domains[0];

    const std::string topoName = topologyName(dom, options);
    const std::string newTopoName = newTopologyName(dom, options);
    const std::string newCSName = newCoordsetName(dom, options);

    const conduit::Node &topologies = dom.fetch_existing("topologies");
    const conduit::Node &topo = topologies.fetch_existing(topoName);
    const conduit::Node *cset = conduit::blueprint::mesh::utils::find_reference_node(topo, "coordset");

    conduit::Node &newTopo = output["topologies/" + newTopoName];
    conduit::Node &newCoordset = output["coordsets/" + newCSName];
    copyState(dom, output);
    execute(topo, *cset, options, newTopo, newCoordset);
  }
}

void
MIRAlgorithm::copyState(const conduit::Node &mesh, conduit::Node &destMesh) const
{
  if(mesh.has_path("state"))
    destMesh["state"].set(mesh["state"]);
}

std::string
MIRAlgorithm::topologyName(const conduit::Node &mesh, const conduit::Node &options) const
{
  std::string topoName;
  if(options.has_path("topology"))
    topoName = options.fetch_existing("topology").as_string();
  else if(mesh.has_path("topologies"))
  {
    const conduit::Node &topologies = mesh.fetch_existing("topologies");
    topoName = topologies[0].name();
  }
  return topoName;
}

std::string
MIRAlgorithm::newTopologyName(const conduit::Node &mesh, const conduit::Node &options) const
{
  std::string topoName;
  if(options.has_path("new_topology"))
    topoName = options.fetch_existing("new_topology").as_string();
  else
    topoName = topologyName(mesh, options);
  return topoName;
}

std::string
MIRAlgorithm::newCoordsetName(const conduit::Node &mesh, const conduit::Node &options) const
{
  std::string csetName;
  if(options.has_path("new_coordset"))
    csetName = options.fetch_existing("new_coordset").as_string();
  else
  {
    std::string topoName = topologyName(mesh, options);
    const conduit::Node &topologies = mesh.fetch_existing("topologies");
    const conduit::Node &topo = topologies.fetch_existing(topoName);
    csetName = topo.fetch_existing("coordset").as_string();
  }

  return csetName;
}

const conduit::Node &MIRAlgorithm::topology(const conduit::Node &mesh, const conduit::Node &options) const
{
  const std::string topoName = topologyName(mesh, options);
  const conduit::Node &topologies = mesh.fetch_existing("topologies");
  return topologies.fetch_existing(topoName);
}

const conduit::Node &MIRAlgorithm::matset(const conduit::Node &mesh, const conduit::Node &options) const
{
  const std::string topoName = topologyName(mesh, options);
  const conduit::Node &matsets = mesh.fetch_existing("matsets");
  for(conduit::index_t i = 0; i < matsets.number_of_children(); i++)
  {
    const conduit::Node &matset = matsets[i];
    if(matset["topology"].as_string() == topoName)
      return matset;
  }
  // We did not find one. TODO: throw exception.
  // return first to eliminate compiler warning.
  return matsets[0];
}

std::vector<std::string>
MIRAlgorithm::fieldNames(const conduit::Node &mesh, const conduit::Node &options) const
{
  std::vector<std::string> names;
  if(options.has_path("fields"))
  {
    const conduit::Node &fields = options["fields"];
    if(fields.number_of_children() > 0)
    {
      for(conduit::index_t i = 0; i < fields.number_of_children(); i++)
        names.push_back(fields[i].name());
    }
    else
    {
      if(fields.dtype().is_char8_str())
        names.push_back(fields.as_string());
    }
  }
  else if(mesh.has_child("fields"))
  {
    const conduit::Node &fields = mesh.fetch_existing("fields");
    for(conduit::index_t i = 0; i < fields.number_of_children(); i++)
      names.push_back(fields[i].name());
  }
  return names;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
namespace elvira
{
} // end namespace elvira

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void ElviraMIRAlgorithm::execute(const conduit::Node &topo,
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

template <typename ExecSpace, typename FuncType>
void for_all_zones(const conduit::Node &topo, FuncType &&func)
{
  conduit::index_t dims[3] = {1, 1, 1};
  conduit::blueprint::mesh::utils::topology::logical_dims(topo, dims, 3);
  const conduit::index_t dimension = conduit::blueprint::mesh::topology::dims(topo);
  const conduit::index_t nzones = conduit::blueprint::mesh::topology::length(topo);

  if(dimension == 1)
  {
    // Line elements
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(axom::IndexType zoneIndex)
      {
        conduit::index_t ids[2];
        ids[0] = zoneIndex;
        ids[1] = zoneIndex + 1;
        func(zoneIndex, ids, 2);
      });
  }
  else if(dimension == 2)
  {
    // Quad elements
    const conduit::index_t nx = dims[0] + 1;
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(axom::IndexType zoneIndex)
      {
        // zoneIndex to i,j
        const conduit::index_t i = zoneIndex % dims[0];
        const conduit::index_t j = zoneIndex / dims[0];
        // Make node ids
        const conduit::index_t jnx = j * nx;
        const conduit::index_t j1nx = jnx + nx;
        conduit::index_t ids[4];
        ids[0] = jnx + i;
        ids[1] = jnx + i + 1;
        ids[2] = j1nx + i + 1;
        ids[3] = j1nx + i;
        func(zoneIndex, ids, 4);
      });
  }
  else if(dimension == 3)
  {
    // Hex elements
    const conduit::index_t nx = dims[0] + 1;
    const conduit::index_t ny = dims[1] + 1;
    const conduit::index_t nz = dims[2] + 1;
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(axom::IndexType zoneIndex)
      {
        // zoneIndex to i,j,k
        const conduit::index_t i = zoneIndex % dims[0];
        const conduit::index_t j = (zoneIndex % (dims[0] * dims[1])) / dims[0];
        const conduit::index_t k = zoneIndex / (dims[0] * dims[1]);

        // Make node ids
        const conduit::index_t knxny  = k * nx * ny;
        const conduit::index_t k1nxny = (k + 1) * nx * ny;
        const conduit::index_t jnx = j * nx;
        const conduit::index_t j1nx = jnx + nx;
        conduit::index_t ids[8];
        ids[0] = knxny  + jnx  + i;
        ids[1] = knxny  + jnx  + i + 1;
        ids[2] = knxny  + j1nx + i + 1;
        ids[3] = knxny  + j1nx + i;
        ids[4] = k1nxny + jnx  + i;
        ids[5] = k1nxny + jnx  + i + 1;
        ids[6] = k1nxny + j1nx + i + 1;
        ids[7] = k1nxny + j1nx + i;

        func(zoneIndex, ids, 4);
      });
  }
  else
  {
//    CONDUIT_ERROR("Unsupported dimension given to traverse_structured "
//        << dimension << ".");
  }
}

template <typename ExecSpace>
void ElviraMIRAlgorithm::executeImpl(const conduit::Node &topo,
                                     const conduit::Node &coordset,
                                     const conduit::Node &options,
                                     conduit::Node &new_topo,
                                     conduit::Node &new_coordset)
{
  // I'm not sure some parts of Elvira work for unstructured topologies. (stencils)
  const auto type = topo["type"].as_string();
  if(type == "unstructured")
  {
    SLIC_ERROR("ElviraMIRAlgorithm cannot be applied to unstructured topologies.");
    return;
  }

  if(options.has_path("zones"))
  {
    const conduit::Node &zones = options.fetch_existing("zones");
    const auto nzones = zones.dtype().number_of_elements();
    detail::IndexNodeToArrayView(options["zones"], [&](auto zonesView)
    {
      axom::for_all<ExecSpace>(
          nzones,
          AXOM_LAMBDA(axom::IndexType i) {
            // Do something with zonesView[i].

            // Get number of materials for zone.
            // for each material, make a plane eq for its interface.

            // for each material, intersect zone with the plane. Make PH fragment and volume, add fragment to new_topo

         });
    });
  }
  else
  {
    foreach_coordset_type // rectilinear, structured, uniform - the problem becomse passing ArrayViews to the lambda. Could we pass a StackArray of ArrayViews?
    {
      // Foreach zone
      for_all_zones<ExecSpace>(topo, AXOM_LAMBDA(conduit::index_t zoneIndex, const conduit::index_t *ids, conduit::index_t nids)
      {
        // on-device

        // might want to pass in the ijk coordinate to deal with rectilinear coordsets. Or perhaps it is best to pass in an array of the coordinate values that we pull out as part of the zone iteration.

      });
    });
  }
}


} // namespace quest
} // namespace axom
