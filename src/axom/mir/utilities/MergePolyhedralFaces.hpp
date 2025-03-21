// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_MERGE_POLYHEDRAL_FACES_HPP_
#define AXOM_MIR_MERGE_POLYHEDRAL_FACES_HPP_

#include "axom/core.hpp"
#include "axom/mir/utilities/utilities.hpp"
#include "axom/mir/utilities/blueprint_utilities.hpp"

#include <conduit/conduit.hpp>

// RAJA
#if defined(AXOM_USE_RAJA)
  #include "RAJA/RAJA.hpp"
#endif

namespace axom
{
namespace mir
{
namespace utilities
{
namespace blueprint
{
/*!
 * \brief Take a Blueprint polyhedral mesh and merge like faces and rewrite the
 *        element and subelement connectivity.
 *
 * \tparam ExecSpace The execution space where the algorithm will execute.
 * \tparam ConnectivityType The type used in the topology connectivity arrays.
 */
template <typename ExecSpace, typename ConnectivityType>
class MergePolyhedralFaces
{
public:
  /*!
   * \brief Merge like faces using the face definitions and rewrite the element
   *        and subelement connectivity. The old definitions are replaced with
   *        new data.
   *
   * \param n_topology The topology to modify.
   */
  static void execute(conduit::Node &n_topo)
  {
    namespace bputils = axom::mir::utilities::blueprint;
    using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;

    SLIC_ASSERT(n_topo["elements/shape"].as_string() == "polyhedral");

    AXOM_ANNOTATE_SCOPE("MergePolyhedralFaces");
    bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;
    const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    // Get the data from the topology and make views.
    conduit::Node &n_elem_conn = n_topo["elements/connectivity"];
    //conduit::Node &n_elem_sizes = n_topo["elements/sizes"];
    //conduit::Node &n_elem_offsets = n_topo["elements/offsets"];
    conduit::Node &n_se_conn = n_topo["subelements/connectivity"];
    conduit::Node &n_se_sizes = n_topo["subelements/sizes"];
    conduit::Node &n_se_offsets = n_topo["subelements/offsets"];
    auto elem_conn = bputils::make_array_view<ConnectivityType>(n_elem_conn);
    //auto elem_sizes = bputils::make_array_view<ConnectivityType>(n_elem_sizes);
    //auto elem_offsets = bputils::make_array_view<ConnectivityType>(n_elem_offsets);
    auto se_conn = bputils::make_array_view<ConnectivityType>(n_se_conn);
    auto se_sizes = bputils::make_array_view<ConnectivityType>(n_se_sizes);
    auto se_offsets = bputils::make_array_view<ConnectivityType>(n_se_offsets);

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("maxnode");
    RAJA::ReduceMax<reduce_policy, ConnectivityType> reduceMaxNodeId(0);
    axom::for_all<ExecSpace>(se_conn.size(), AXOM_LAMBDA(axom::IndexType index)
    {
      reduceMaxNodeId.max(se_conn[index]);
    });
    const auto maxNodeId = reduceMaxNodeId.get();
    AXOM_ANNOTATE_END("maxnode");

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("naming");
    using NamingType = HashNaming<ConnectivityType, 5>;
    using KeyType = typename NamingType::KeyType;

    NamingType naming;
    naming.setMaxId(maxNodeId);
    auto namingView = naming.view();

    const axom::IndexType totalFaces = se_sizes.size();
    axom::Array<KeyType> faceNames(totalFaces, totalFaces, allocatorID);
    auto faceNamesView = faceNames.view();
    axom::for_all<ExecSpace>(totalFaces, AXOM_LAMBDA(axom::IndexType faceIndex)
    {
      // Get size and offset for current face.
      const auto numFaceIds = se_sizes[faceIndex];
      const auto faceIds = se_conn.data() + se_offsets[faceIndex];

      // Make a name for this face.
      faceNamesView[faceIndex] = namingView.makeName(faceIds, numFaceIds);
    });
    AXOM_ANNOTATE_END("naming");

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("unique");
    // The unique keys.
    axom::Array<std::uint64_t> uniqueKeys;
    // The index of the sorted key in the original keys. We can use this to
    // determine which face to use.
    axom::Array<axom::IndexType> selectedFaces;

    // Make faces unique.
    axom::mir::utilities::Unique<ExecSpace, std::uint64_t>::execute(faceNamesView, uniqueKeys, selectedFaces);
    const auto uniqueKeysView = uniqueKeys.view();
    const auto selectedFacesView = selectedFaces.view();
    AXOM_ANNOTATE_END("unique");

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("rewriting");
    conduit::Node n_new_se_sizes;
    n_new_se_sizes.set_allocator(c2a.getConduitAllocatorID());
    n_new_se_sizes.set(
      conduit::DataType(bputils::cpp2conduit<ConnectivityType>::id, selectedFaces.size()));
    auto new_se_sizes = bputils::make_array_view<ConnectivityType>(n_new_se_sizes);

    conduit::Node n_new_se_offsets;
    n_new_se_offsets.set_allocator(c2a.getConduitAllocatorID());
    n_new_se_offsets.set(
      conduit::DataType(bputils::cpp2conduit<ConnectivityType>::id, selectedFaces.size()));
    auto new_se_offsets = bputils::make_array_view<ConnectivityType>(n_new_se_offsets);

    // Copy the sizes of the selected faces into new_se_sizes and make new_se_offsets
    RAJA::ReduceSum<reduce_policy, axom::IndexType> reduceNewSizes(0);
    axom::for_all<ExecSpace>(selectedFaces.size(), AXOM_LAMBDA(axom::IndexType index)
    {
      const auto size =  se_sizes[selectedFacesView[index]];
      new_se_sizes[index] = size;
      reduceNewSizes += size;
    });
    const axom::IndexType newSEConnSize = reduceNewSizes.get();
    axom::exclusive_scan<ExecSpace>(new_se_sizes, new_se_offsets);

    // Allocate new_se_conn to contain the new face definitions.
    conduit::Node n_new_se_conn;
    n_new_se_conn.set_allocator(c2a.getConduitAllocatorID());
    n_new_se_conn.set(
      conduit::DataType(bputils::cpp2conduit<ConnectivityType>::id, newSEConnSize));
    auto new_se_conn = bputils::make_array_view<ConnectivityType>(n_new_se_conn);

    // Copy the selected faces into new_se_conn.
    axom::for_all<ExecSpace>(selectedFaces.size(), AXOM_LAMBDA(axom::IndexType index)
    {
      const auto numFaceIds = new_se_sizes[index];
      const auto destOffset = new_se_offsets[index];
      const auto srcOffset = se_offsets[selectedFacesView[index]];
      for(axom::IndexType i = 0; i < numFaceIds; i++)
      {
        new_se_conn[destOffset + i] = se_conn[srcOffset + i];
      }
    });

    // Now, rewrite the connectivity with the new face ids.
    axom::for_all<ExecSpace>(elem_conn.size(), AXOM_LAMBDA(axom::IndexType index)
    {
      const auto originalFaceId = elem_conn[index];

      // Get the "name" of the old face.
      const auto originalFaceKey = faceNamesView[originalFaceId];

      // Look for the index of the "name" in the new uniqueKeys.
      // That will be its face index in the new faces.
      const auto newId = axom::mir::utilities::bsearch(originalFaceKey, uniqueKeysView);
      SLIC_ASSERT(newId != -1);
      elem_conn[index] = static_cast<ConnectivityType>(newId);
    });

    // Move the "new" nodes into the Blueprint hierarchy.
    n_se_conn.move(n_new_se_conn);
    n_se_sizes.move(n_new_se_sizes);
    n_se_offsets.move(n_new_se_offsets);

    AXOM_ANNOTATE_END("rewriting");
  }
};

}  // end namespace blueprint
}  // end namespace utilities
}  // end namespace mir
}  // end namespace axom

#endif
