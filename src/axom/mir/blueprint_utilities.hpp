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
#include "axom/mir/views/dispatch_structured_topology.hpp"
#include "axom/mir/views/NodeArrayView.hpp"

#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint.hpp>

#if defined(AXOM_USE_RAJA)
  #include <RAJA/RAJA.hpp>
#endif

#include <limits>
#include <utility>
#include <string>

namespace axom
{
namespace mir
{
namespace utilities
{
namespace blueprint
{
//------------------------------------------------------------------------------
/**
 * \brief This class provides a couple of type traits that let us map C++ types
 *        to types / values useful in Conduit.
 */
template <typename T>
struct cpp2conduit
{
  static constexpr conduit::index_t type = conduit::DataType::EMPTY_ID;
};

template <>
struct cpp2conduit<conduit::int8>
{
  using type = conduit::int8;
  static constexpr conduit::index_t id = conduit::DataType::INT8_ID;
};

template <>
struct cpp2conduit<conduit::int16>
{
  using type = conduit::int16;
  static constexpr conduit::index_t id = conduit::DataType::INT16_ID;
};

template <>
struct cpp2conduit<conduit::int32>
{
  using type = conduit::int32;
  static constexpr conduit::index_t id = conduit::DataType::INT32_ID;
};

template <>
struct cpp2conduit<conduit::int64>
{
  using type = conduit::int64;
  static constexpr conduit::index_t id = conduit::DataType::INT64_ID;
};

template <>
struct cpp2conduit<conduit::uint8>
{
  using type = conduit::uint8;
  static constexpr conduit::index_t id = conduit::DataType::UINT8_ID;
};

template <>
struct cpp2conduit<conduit::uint16>
{
  using type = conduit::uint16;
  static constexpr conduit::index_t id = conduit::DataType::UINT16_ID;
};

template <>
struct cpp2conduit<conduit::uint32>
{
  using type = conduit::uint32;
  static constexpr conduit::index_t id = conduit::DataType::UINT32_ID;
};

template <>
struct cpp2conduit<conduit::uint64>
{
  using type = conduit::uint64;
  static constexpr conduit::index_t id = conduit::DataType::UINT64_ID;
};

template <>
struct cpp2conduit<conduit::float32>
{
  using type = conduit::float32;
  static constexpr conduit::index_t id = conduit::DataType::FLOAT32_ID;
};

template <>
struct cpp2conduit<conduit::float64>
{
  using type = conduit::float64;
  static constexpr conduit::index_t id = conduit::DataType::FLOAT64_ID;
};

//------------------------------------------------------------------------------
/**
 * \brief This class registers a Conduit allocator that can make Conduit allocate
 *        through Axom's allocate/deallocate functions using a specific allocator. This
 *        permits Conduit to allocate through Axom's UMPIRE logic.
 *
 */
template <typename ExecSpace>
class ConduitAllocateThroughAxom
{
public:
  static conduit::index_t getConduitAllocatorID()
  {
    static conduit::index_t conduitAllocatorID = -1;
    if(conduitAllocatorID == -1)
    {
      conduitAllocatorID =
        conduit::utils::register_allocator(internal_allocate, internal_free);
    }
    return conduitAllocatorID;
  }

private:
  static void *internal_allocate(size_t items, size_t item_size)
  {
    const auto axomAllocatorID = axom::execution_space<ExecSpace>::allocatorID();
    void *ptr = static_cast<void *>(
      axom::allocate<std::uint8_t>(items * item_size, axomAllocatorID));
    //std::cout << axom::execution_space<ExecSpace>::name() << ": Allocated for Conduit via axom: items=" << items << ", item_size=" << item_size << ", ptr=" << ptr << std::endl;
    return ptr;
  }

  static void internal_free(void *ptr)
  {
    //std::cout << axom::execution_space<ExecSpace>::name() << ": Dellocating for Conduit via axom: ptr=" << ptr << std::endl;
    axom::deallocate(ptr);
  }
};

//------------------------------------------------------------------------------
/**
 * \accelerated
 * \brief Make an unstructured representation of a structured topology.
 *
 * \tparam ExecSpace The execution space where the work will be done.
 *
 * \param topo The input topology to be turned into unstructured.
 * \param coordset The topology's coordset. It will be referenced as an external node in the output \a mesh.
 * \param topoName The name of the new topology to create.
 * \param mesh     The node that will contain the new topology and coordset.
 *
 * \note There are blueprint methods for this sort of thing but this one is accelerated.
 */
template <typename ExecSpace>
void to_unstructured(const conduit::Node &topo,
                     const conduit::Node &coordset,
                     const std::string &topoName,
                     conduit::Node &mesh)
{
  const std::string type = topo.fetch_existing("type").as_string();
  ConduitAllocateThroughAxom<ExecSpace> c2a;

  mesh["coordsets"][coordset.name()].set_external(coordset);
  conduit::Node &newtopo = mesh["topologies"][topoName];
  newtopo["coordset"] = coordset.name();

  if(type == "unstructured")
  {
    newtopo.set_external(topo);
  }
  else
  {
    newtopo["type"] = "unstructured";
    conduit::Node &n_newconn = newtopo["elements/connectivity"];
    conduit::Node &n_newsizes = newtopo["elements/sizes"];
    conduit::Node &n_newoffsets = newtopo["elements/offsets"];
    n_newconn.set_allocator(c2a.getConduitAllocatorID());
    n_newsizes.set_allocator(c2a.getConduitAllocatorID());
    n_newoffsets.set_allocator(c2a.getConduitAllocatorID());

    views::dispatch_structured_topologies(
      topo,
      coordset,
      [&](const std::string &shape, auto &topoView) {
        int ptsPerZone = 2;
        if(shape == "quad")
          ptsPerZone = 4;
        else if(shape == "hex")
          ptsPerZone = 8;

        newtopo["elements/shape"] = shape;

        // Allocate new mesh data.
        const auto nzones = topoView.numberOfZones();
        const auto connSize = nzones * ptsPerZone;
        n_newconn.set(conduit::DataType::index_t(connSize));
        n_newsizes.set(conduit::DataType::index_t(nzones));
        n_newoffsets.set(conduit::DataType::index_t(nzones));

        // Make views for the mesh data.
        axom::ArrayView<conduit::index_t> connView(
          static_cast<conduit::index_t *>(n_newconn.data_ptr()),
          connSize);
        axom::ArrayView<conduit::index_t> sizesView(
          static_cast<conduit::index_t *>(n_newsizes.data_ptr()),
          nzones);
        axom::ArrayView<conduit::index_t> offsetsView(
          static_cast<conduit::index_t *>(n_newoffsets.data_ptr()),
          nzones);

        // Fill in the new connectivity.
        topoView.template for_all_zones<ExecSpace>(
          AXOM_LAMBDA(auto zoneIndex, const auto &zone) {
            const auto start = zoneIndex * ptsPerZone;
            for(int i = 0; i < ptsPerZone; i++)
              connView[start + i] =
                static_cast<conduit::index_t>(zone.getIds()[i]);

            sizesView[zoneIndex] = ptsPerZone;
            offsetsView[zoneIndex] = start;
          });
      });
  }
}

/**
 * \brief Copies a Conduit tree in the \a src node to a new Conduit \a dest node,
 *        making sure to allocate array data in the appropriate memory space for
 *        the execution space.
 */
template <typename ExecSpace>
void copy(conduit::Node &dest, const conduit::Node &src)
{
  ConduitAllocateThroughAxom<ExecSpace> c2a;

  if(src.number_of_children() > 0)
  {
    for(conduit::index_t i = 0; i < src.number_of_children(); i++)
    {
      copy<ExecSpace>(dest[src[i].name()], src[i]);
    }
  }
  else
  {
    if(!src.dtype().is_string() && src.dtype().number_of_elements() > 1)
    {
      // Allocate the node's memory in the right place.
      dest.reset();
      dest.set_allocator(c2a.getConduitAllocatorID());
      dest.set(
        conduit::DataType(src.dtype().id(), src.dtype().number_of_elements()));

      // Copy the data to the destination node. Axom uses Umpire to manage that.
      if(src.is_compact())
        axom::copy(dest.data_ptr(), src.data_ptr(), src.dtype().bytes_compact());
      else
      {
        // NOTE: this assumes that src is on the host. Why would we have strided data on device?
        conduit::Node tmp;
        src.compact_to(tmp);
        axom::copy(dest.data_ptr(), tmp.data_ptr(), tmp.dtype().bytes_compact());
      }
    }
    else
    {
      // The node data fits in the node. It's on the host.
      dest.set(src);
    }
  }
}

/**
 * \brief Get the min/max values for the data in a Conduit node.
 *
 * \param[in] n The Conduit node whose data we're checking.
 *
 * \return A pair containing the min,max values in the node.
 */
template <typename ExecSpace>
std::pair<double, double> minmax(const conduit::Node &n)
{
  std::pair<double, double> retval {0., 0.};

  axom::mir::views::Node_to_ArrayView(n, [&](auto nview) {
    using value_type = typename decltype(nview)::value_type;
    using reduce_policy =
      typename axom::execution_space<ExecSpace>::reduce_policy;

    RAJA::ReduceMin<reduce_policy, value_type> vmin(
      std::numeric_limits<value_type>::max());
    RAJA::ReduceMax<reduce_policy, value_type> vmax(
      std::numeric_limits<value_type>::min());

    axom::for_all<ExecSpace>(
      nview.size(),
      AXOM_LAMBDA(auto index) {
        vmin.min(nview[index]);
        vmax.max(nview[index]);
      });

    retval = std::pair<double, double> {vmin.get(), vmax.get()};
  });

  return retval;
}

/**
 * \brief Save a Blueprint mesh to a legacy ASCII VTK file.
 *
 * \param node The node that contains the mesh data.
 * \param path The file path to save.
 *
 * \note This function currently handles only unstructured topos with explicit coordsets.
 */
void save_vtk(const conduit::Node &node, const std::string &path);

}  // end namespace blueprint
}  // end namespace utilities
}  // end namespace mir
}  // end namespace axom

#endif
