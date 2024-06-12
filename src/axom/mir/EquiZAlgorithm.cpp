#include "axom/mir/EquiZAlgorithm.hpp"


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
namespace mir
{

#if 0
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
#endif

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
#if 0
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
#endif
  }
  else
  {
#if 0
    foreach_coordset_type // rectilinear, structured, uniform - the problem becomse passing ArrayViews to the lambda. Could we pass a StackArray of ArrayViews?
    {
      // Foreach zone
      for_all_zones<ExecSpace>(topo, AXOM_LAMBDA(conduit::index_t zoneIndex, const conduit::index_t *ids, conduit::index_t nids)
      {
        // on-device

        // might want to pass in the ijk coordinate to deal with rectilinear coordsets. Or perhaps it is best to pass in an array of the coordinate values that we pull out as part of the zone iteration.

      });
    });
#endif
  }
}

} // namespace mir
} // namespace axom
