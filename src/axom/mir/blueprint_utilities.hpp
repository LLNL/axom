// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_BLUEPRINT_UTILITIES_HPP_
#define AXOM_MIR_BLUEPRINT_UTILITIES_HPP_

#include "axom/core/execution/execution_space.hpp"
#include "axom/core/Array.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/core/NumericLimits.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/mir/views/dispatch_structured_topology.hpp"
#include "axom/mir/views/NodeArrayView.hpp"

#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint.hpp>

#if defined(AXOM_USE_RAJA)
  #include <RAJA/RAJA.hpp>
#endif

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
{ };

template <>
struct cpp2conduit<conduit::int8>
{
  using type = conduit::int8;
  static constexpr conduit::index_t id = conduit::DataType::INT8_ID;
  static const char *name;
};

template <>
struct cpp2conduit<conduit::int16>
{
  using type = conduit::int16;
  static constexpr conduit::index_t id = conduit::DataType::INT16_ID;
  static const char *name;
};

template <>
struct cpp2conduit<conduit::int32>
{
  using type = conduit::int32;
  static constexpr conduit::index_t id = conduit::DataType::INT32_ID;
  static const char *name;
};

template <>
struct cpp2conduit<conduit::int64>
{
  using type = conduit::int64;
  static constexpr conduit::index_t id = conduit::DataType::INT64_ID;
  static const char *name;
};

template <>
struct cpp2conduit<conduit::uint8>
{
  using type = conduit::uint8;
  static constexpr conduit::index_t id = conduit::DataType::UINT8_ID;
  static const char *name;
};

template <>
struct cpp2conduit<conduit::uint16>
{
  using type = conduit::uint16;
  static constexpr conduit::index_t id = conduit::DataType::UINT16_ID;
  static const char *name;
};

template <>
struct cpp2conduit<conduit::uint32>
{
  using type = conduit::uint32;
  static constexpr conduit::index_t id = conduit::DataType::UINT32_ID;
  static const char *name;
};

template <>
struct cpp2conduit<conduit::uint64>
{
  using type = conduit::uint64;
  static constexpr conduit::index_t id = conduit::DataType::UINT64_ID;
  static const char *name;
};

template <>
struct cpp2conduit<conduit::float32>
{
  using type = conduit::float32;
  static constexpr conduit::index_t id = conduit::DataType::FLOAT32_ID;
  static const char *name;
};

template <>
struct cpp2conduit<conduit::float64>
{
  using type = conduit::float64;
  static constexpr conduit::index_t id = conduit::DataType::FLOAT64_ID;
  static const char *name;
};

//------------------------------------------------------------------------------

/**
 * \brief Make an axom::ArrayView from a Conduit node.
 *
 * \tparam T The type for the array view elements.
 *
 * \param n The conduit node for which we want an array view.
 *
 * \return An axom::ArrayView that wraps the data in the Conduit node.
 */
/// @{
template <typename T>
inline axom::ArrayView<T> make_array_view(conduit::Node &n)
{
  SLIC_ASSERT_MSG(
    cpp2conduit<T>::id == n.dtype().id(),
    axom::fmt::format("Cannot create ArrayView<{}> for Conduit {} data.",
                      cpp2conduit<T>::name,
                      n.dtype().name()));
  return axom::ArrayView<T>(static_cast<T *>(n.data_ptr()),
                            n.dtype().number_of_elements());
}

template <typename T>
inline axom::ArrayView<T> make_array_view(const conduit::Node &n)
{
  SLIC_ASSERT_MSG(
    cpp2conduit<T>::id == n.dtype().id(),
    axom::fmt::format("Cannot create ArrayView<{}> for Conduit {} data.",
                      cpp2conduit<T>::name,
                      n.dtype().name()));
  return axom::ArrayView<T>(static_cast<T *>(const_cast<void *>(n.data_ptr())),
                            n.dtype().number_of_elements());
}
/// @}

//------------------------------------------------------------------------------
/**
 * \brief This class registers a Conduit allocator that can make Conduit allocate
 *        through Axom's allocate/deallocate functions using a specific allocator.
 *        This permits Conduit to allocate through Axom's UMPIRE logic.
 *
 * \tparam ExecSpace The execution space.
 */
template <typename ExecSpace>
class ConduitAllocateThroughAxom
{
public:
  /**
   * \brief Get the Conduit allocator ID for this ExecSpace.
   *
   * \return The Conduit allocator ID for this ExecSpace.
   */
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
  /**
   * \brief A function we register with Conduit to allocate memory.
   *
   * \param items The number of items to allocate.
   * \param item_size The size of each item in bytes.
   *
   * \brief A block of newly allocated memory large enough for the requested items.
   */
  static void *internal_allocate(size_t items, size_t item_size)
  {
    const auto axomAllocatorID = axom::execution_space<ExecSpace>::allocatorID();
    void *ptr = static_cast<void *>(
      axom::allocate<std::uint8_t>(items * item_size, axomAllocatorID));
    //std::cout << axom::execution_space<ExecSpace>::name() << ": Allocated for Conduit via axom: items=" << items << ", item_size=" << item_size << ", ptr=" << ptr << std::endl;
    return ptr;
  }

  /**
   * \brief A deallocation function we register with Conduit.
   */
  static void internal_free(void *ptr)
  {
    //std::cout << axom::execution_space<ExecSpace>::name() << ": Dellocating for Conduit via axom: ptr=" << ptr << std::endl;
    axom::deallocate(ptr);
  }
};

//------------------------------------------------------------------------------
/**
 * \brief Fill an array with int values from a Conduit node.
 *
 * \tparam ArrayType The array type being filled.
 *
 * \param n The node that contains the data.
 * \param key The name of the node that contains the data in \a n.
 * \param[out] arr The array being filled.
 */
template <typename ArrayType>
bool fillFromNode(const conduit::Node &n, const std::string &key, ArrayType &arr)
{
  bool found = false;
  if((found = n.has_path(key)) == true)
  {
    const auto acc = n.fetch_existing(key).as_int_accessor();
    for(int i = 0; i < arr.size(); i++)
    {
      arr[i] = acc[i];
    }
  }
  return found;
}
//------------------------------------------------------------------------------

/**
 * \brief Returns the input index (no changes).
 */
struct DirectIndexing
{
  /**
   * \brief Return the input index (no changes).
   * \param index The input index.
   * \return The input index.
   */
  AXOM_HOST_DEVICE
  inline axom::IndexType operator[](axom::IndexType index) const
  {
    return index;
  }
};

//------------------------------------------------------------------------------
/**
 * \brief Help turn slice data zone indices into strided structured element field indices.
 * \tparam Indexing A StridedStructuredIndexing of some dimension.
 */
template <typename Indexing>
struct SSElementFieldIndexing
{
  /**
   * \brief Update the indexing offsets/strides from a Conduit node.
   * \param field The Conduit node for a field.
   */
  void update(const conduit::Node &field)
  {
    fillFromNode(field, "offsets", m_indexing.m_offsets);
    fillFromNode(field, "strides", m_indexing.m_strides);
  }

  /**
   * \brief Transforms the index from local to global through an indexing object.
   * \param index The local index
   * \return The global index for the field.
   */
  AXOM_HOST_DEVICE
  inline axom::IndexType operator[](axom::IndexType index) const
  {
    return m_indexing.LocalToGlobal(index);
  }

  Indexing m_indexing {};
};

//------------------------------------------------------------------------------
/**
 * \brief Help turn blend group node indices (global) into vertex field indices.
 * \tparam Indexing A StridedStructuredIndexing of some dimension.
 */
template <typename Indexing>
struct SSVertexFieldIndexing
{
  /**
   * \brief Update the indexing offsets/strides from a Conduit node.
   * \param field The Conduit node for a field.
   */
  void update(const conduit::Node &field)
  {
    fillFromNode(field, "offsets", m_fieldIndexing.m_offsets);
    fillFromNode(field, "strides", m_fieldIndexing.m_strides);
  }

  /**
   * \brief Transforms the index from local to global through an indexing object.
   * \param index The global index
   * \return The global index for the field.
   */
  AXOM_HOST_DEVICE
  inline axom::IndexType operator[](axom::IndexType index) const
  {
    // Make the global index into a global logical in the topo.
    const auto topoGlobalLogical = m_topoIndexing.GlobalToGlobal(index);
    // Make the global logical into a local logical in the topo.
    const auto topoLocalLogical = m_topoIndexing.GlobalToLocal(topoGlobalLogical);
    // Make the global logical index in the field.
    const auto fieldGlobalLogical =
      m_fieldIndexing.LocalToGlobal(topoLocalLogical);
    // Make the global index in the field.
    const auto fieldGlobalIndex =
      m_fieldIndexing.GlobalToGlobal(fieldGlobalLogical);
    return fieldGlobalIndex;
  }

  Indexing m_topoIndexing {};
  Indexing m_fieldIndexing {};
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

    axom::mir::views::dispatch_structured_topologies(
      topo,
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
        auto connView = make_array_view<conduit::index_t>(n_newconn);
        auto sizesView = make_array_view<conduit::index_t>(n_newsizes);
        auto offsetsView = make_array_view<conduit::index_t>(n_newoffsets);

        // Fill in the new connectivity.
        topoView.template for_all_zones<ExecSpace>(
          AXOM_LAMBDA(auto zoneIndex, const auto &zone) {
            const auto start = zoneIndex * ptsPerZone;
            for(int i = 0; i < ptsPerZone; i++)
            {
              connView[start + i] =
                static_cast<conduit::index_t>(zone.getIds()[i]);
            }
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
 *
 * \tparam The destination execution space (e.g. axom::SEQ_EXEC).
 *
 * \param dest The conduit node that will receive the copied data.
 * \param src The source data to be copied.
 */
template <typename ExecSpace>
void copy(conduit::Node &dest, const conduit::Node &src)
{
  ConduitAllocateThroughAxom<ExecSpace> c2a;
  dest.reset();
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
        // NOTE: This assumes that src is on the host.
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
template <typename ExecSpace, typename ReturnType>
std::pair<ReturnType, ReturnType> minmax(const conduit::Node &n)
{
  SLIC_ASSERT(n.dtype().number_of_elements() > 0);
  std::pair<ReturnType, ReturnType> retval;

  axom::mir::views::Node_to_ArrayView(n, [&](auto nview) {
    using value_type = typename decltype(nview)::value_type;
    using reduce_policy =
      typename axom::execution_space<ExecSpace>::reduce_policy;

    RAJA::ReduceMin<reduce_policy, value_type> vmin(
      axom::numeric_limits<value_type>::max());
    RAJA::ReduceMax<reduce_policy, value_type> vmax(
      axom::numeric_limits<value_type>::min());

    axom::for_all<ExecSpace>(
      nview.size(),
      AXOM_LAMBDA(auto index) {
        vmin.min(nview[index]);
        vmax.max(nview[index]);
      });

    retval =
      std::pair<ReturnType, ReturnType> {static_cast<ReturnType>(vmin.get()),
                                         static_cast<ReturnType>(vmax.get())};
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
