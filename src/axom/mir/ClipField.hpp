// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_CLIP_FIELD_HPP_
#define AXOM_MIR_CLIP_FIELD_HPP_

#include "axom/core.hpp"
#include "axom/mir/FieldBlender.hpp"
#include "axom/mir/CoordsetBlender.hpp"
#include "axom/mir/FieldSlicer.hpp"
#include "axom/mir/blueprint_utilities.hpp" // for cpp2conduit
#include "axom/mir/utilities.hpp" // for cpp2conduit
#include "axom/primal/Point.hpp"

#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint_mesh_utils.hpp>

namespace axom
{
namespace mir
{
namespace clipping
{
namespace details
{
template <typename KeyType, typename ValueType>
std::map<ValueType, KeyType> reverse_map(const std::map<KeyType, ValueType> &m)
{
  std::map<ValueType, KeyType> output;
  for(const auto& [key, value] : m)
  {
    output[value] = key;
  }
  return output;
}

std::map<std::string, int>
shapeMap_NameValue(const conduit::Node &n_shape_map)
{
  std::map<std::string, int> sm;
  for(conduit::index_t i = 0; i < n_shape_map.number_of_children(); i++)
  {
    sm[n_shape_map[i].name()] = n_shape_map[i].to_int();
  }
  return sm;
}

std::map<int, std::string>
shapeMap_ValueName(const conduit::Node &n_shape_map)
{
  std::map<int, std::string> sm;
  for(conduit::index_t i = 0; i < n_shape_map.number_of_children(); i++)
  {
    sm[n_shape_map[i].to_int()] = n_shape_map[i].name();
  }
  return sm;
}

std::map<std::string, int>
shapeMap_FromFlags(uint64 shapes)
{
  std::map<std::string, int> sm;

  if((shapes & LineShape<int>::id()) > 0)
    sm["line"] = LineShape<int>::id();

  if((shapes & TriShape<int>::id()) > 0)
    sm["tri"] = TriShape<int>::id();

  if((shapes & QuadShape<int>::id()) > 0)
    sm["quad"] = QuadShape<int>::id();

  if((shapes & TetShape<int>::id()) > 0)
    sm["tet"] = TetShape<int>::id();

  if((shapes & PyramidShape<int>::id()) > 0)
    sm["pyramid"] = PyramidShape<int>::id();

  if((shapes & WedgeShape<int>::id()) > 0)
    sm["wedge"] = WedgeShape<int>::id();

  if((shapes & HexShape<int>::id()) > 0)
    sm["hex"] = HexShape<int>::id();

  return sm;
}

AXOM_HOST_DEVICE
constexpr int getClipTableIndex(int dimension, int nnodes)
{
  return (dimension == 2) ? ((nnodes == 3) ? 0 : 1) : (nnodes - 2);
}
AXOM_HOST_DEVICE
constexpr bool color0Selected(int selection)
{
  return axom::utilities::bitIsSet(selection, 0);
}

AXOM_HOST_DEVICE
constexpr bool color1Selected(int selection)
{
  return axom::utilities::bitIsSet(selection, 1);
}

AXOM_HOST_DEVICE
constexpr bool generatedPointIsSelected(unsigned char color, int selection)
{
  return color == NOCOLOR ||
         (color0Selected(selection) && color == COLOR0) ||
         (color1Selected(selection) && color == COLOR1);
}

AXOM_HOST_DEVICE
constexpr bool shapeIsSelected(unsigned char color, int selection)
{
  return (color0Selected(selection) && color == COLOR0) ||
         (color1Selected(selection) && color == COLOR1);
}

AXOM_HOST_DEVICE
float computeT(float d0, float d1)
{
  const float delta = d1 - d0;
  const float abs_delta = (delta < 0) ? -delta : delta;
  const float t = (abs_delta != 0.) ? (-d0 / delta) : 0.;
  return t;
}

AXOM_HOST_DEVICE
template <typename ZoneType, typename ArrayViewType>
size_t clip_case(const ZoneType &zone, const ArrayViewType &view)
{
  size_t clipcase = 0;
  for(size_t i = 0; i < zone.numberOfNodes(); i++)
  {
    const auto id = zone.getId(i);
    if(view[id] > 0)
    {
      clipcase |= (1 << j);
    }
  }
  return clipcase;
}

//------------------------------------------------------------------------------
template <typename ZoneType>
constexpr size_t indexOfZoneType() { return 0; }

template <>
constexpr size_t indexOfZoneType<TriShape>() { return 0; }

template <>
constexpr size_t indexOfZoneType<QuadShape>() { return 1; }

template <>
constexpr size_t indexOfZoneType<TetShape>() { return 2; }

template <>
constexpr size_t indexOfZoneType<PyramidShape>() { return 3; }

template <>
constexpr size_t indexOfZoneType<WedgeShape>() { return 4; }

template <>
constexpr size_t indexOfZoneType<HexShape>() { return 5; }

} // end namespace details

/**
 * \accelerated
 * \brief This class clips a topology using a field and puts the new topology into a new Conduit node.
 */
template <typename ExecSpace, typename TopologyView, typename CoordsetView>
class ClipField
{
public:
/*
  void execute(const conduit::Node &inputMesh, const std::string &clipField, conduit::Node &outputMesh)
  {

 I was going to say that we could create the views here but there's a lot of dispatch code devoted to making typed array views, etc that would need to be done before we even know the types.

 I could include the dispatch_coordset method for example but then this class would have all the various instantiations of the algorithm rather than just the one we want.

    CoordsetView csview;
  }
 */

  void execute(const TopologyView &topoView, // I'd rather just pass the views to the method.
               const CoordsetView &coordsetView,
               const conduit::Node &n_inputMesh,
               const std::string &clipField,
               conduit::Node &n_outputMesh)
  {
    using KeyType = uint64;
    using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
    using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;
    const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    // Get the clip field.
    conduit::Node &n_clip_field = n_inputMesh["fields"].fetch_existing(clipField);

    // TODO: process options (or add some class members to hold attributes)
    int selection = -1; // All bits set.
    

    // Load clip table data and make views.
    ClipTableManager<ExecSpace> clipTables;
    clipTables.load(topoView.dimension());
    axom::StackArray<ClipTableManager<ExecSpace>::View, 6> clipTableViews;
    if(topoView.dimension() == 2)
    {
      clipTableViews[details::indexOfZoneType<TriShape>()] = clipTables[ST_TRI].view();
      clipTableViews[details::indexOfZoneType<QuadShape>()] = clipTables[ST_QUA].view();
    }
    if(topoView.dimension() == 2)
    {
      clipTableViews[details::indexOfZoneType<TetShape>()] = clipTables[ST_TET].view();
      clipTableViews[details::indexOfZoneType<PyramidShape>()] = clipTables[ST_PYR].view();
      clipTableViews[details::indexOfZoneType<WedgeShape>()] = clipTables[ST_WDG].view();
      clipTableViews[details::indexOfZoneType<HexShape>()] = clipTables[ST_HEX].view();
    }

    // ----------------------------------------------------------------------
    //
    // Stage 1: Iterate over elements and their respective clip cases to
    //          determine sizes of outputs.
    //
    // ----------------------------------------------------------------------
    RAJA::ReduceSum<reduce_policy, int> fragment_sum(0);
    RAJA::ReduceSum<reduce_policy, int> fragment_nids_sum(0);
    RAJA::ReduceSum<reduce_policy, int> blendGroups_sum(0);
    RAJA::ReduceSum<reduce_policy, int> blendGroupLen_sum(0);

    // Allocate some memory
    const auto nzones = topoView.numberOfZones();
    axom::Array<int> clipCases(nzones, nzones, allocatorID);      // The clip case for a zone.
    axom::Array<int64> pointsUsed(nzones, nzones, allocatorID);   // Which points are used over all selected fragments in a zone
    axom::Array<int> blendGroups(nzones, nzones, allocatorID);    // Number of blend groups in a zone.
    axom::Array<int> blendGroupsLen(nzones, nzones, allocatorID); // Length of the blend groups in a zone.
    axom::Array<int> fragments(nzones, nzones, allocatorID);      // The number of fragments (child zones) produced for a zone.
    axom::Array<int> fragmentsSize(nzones, nzones, allocatorID);  // The connectivity size for all selected fragments in a zone.

    auto clipCasesView = clipCases.view();
    auto pointsUsedView = pointsUsed.view();
    auto blendGroupsview = blendGroups.view();
    auto blendGroupsLenView = blendGroupsLen.view();
    auto fragmentsView = fragments.view();
    auto fragmentsSizeView = fragmentsSize.view();

    views::Node_to_ArrayView(n_clip_field, [&](auto clipFieldView)
    {
      topoView.template for_all_zones<ExecSpace>(AXOM_LAMBDA(auto zoneIndex, const auto &zone)
      {
        // Get the clip case for the current zone.
        const auto clipcase = details::clip_case(zone, clipFieldView);
        clipCasesView[zoneIndex] = clipcase;

        // Iterate over the shapes in this clip case to determine the number of blend groups.
        const auto clipTableIndex = details::getClipTableIndex(zone.dimension(), zone.numberOfNodes());
        const auto &ctView = clipTableViews[clipTableIndex];

        int thisBlendGroups = 0;    // The number of blend groups produced in this case.
        int thisBlendGroupLen = 0;  // The total length of the blend groups.
        int thisFragments = 0;      // The number of zone fragments produced in this case.
        int thisFragmentsNumIds = 0;// The number of points used to make all the fragment zones.
        int64 ptused = 0;           // A bitset indicating which ST_XX nodes are used.

        auto it = ctView.begin(clipcase);
        const auto end = ctView.end(clipcase);
        for(; it != end; it++)
        {
          // Get the current shape in the clip case.
          const auto caseData = *it;

          if(caseData[0] == ST_PNT)
          {
            if(details::generatedPointIsSelected(caseData[2], selection))
            {
              const size_t nIds = caseData[3];
              for(size_t ni = 4; ni < nIds; ni++)
              {
                const auto pid = caseData[ni];

                // Increase the blend size to include this center point.
                if(pid <= P7)
                {
                   // corner point.
                   thisBlendGroupLen++;
                }
                else if(pid >= EA && pid <= EL)
                {
                  // edge points are derived from 2 corner points. If
                  // those appear here then we're probably creating a
                  // face point. We can store the 2 corner points in place
                  // of the edge point (along with some blending coeff).
                  thisBlendGroupLen += 2;
                }
              }

              // This center or face point counts as a blend group.
              thisBlendGroups++;

              // Mark the point used.
              axom::utilities::setBit(ptused, N0 + caseData[1]);
            }
          }
          else
          {
            if(details::shapeIsSelected(caseData[1], selection))
            {
              thisFragments++;
              const int nIdsThisShape = caseData.size() - 2;
              thisFragmentsNumIds += nIdsThisShape;

              // Mark the points in this fragment used.
              for(size_t i = 2; i < caseData.size(); i++)
              {
                axom::utilities::setBit(ptused, caseData[i]);
              }
            }
          }
        }

        // Save the flags for the points that were used in this zone
        pointsUsedView[zoneIndex] = ptused;

        // Count which points in the original cell are used.
        for(unsigned char pid = P0; pid <= P7; pid++)
        {
          const int incr = axom::utilities::bitIsSet(ptused, pid) ? 1 : 0;
          thisBlendGroupLen += incr; // {p0}
          thisBlendGroups += incr;
        }

        // Count edges that are used.
        for(unsigned char pid = EA; pid <= EL; pid++)
        {
          const int incr = axom::utilities::bitIsSet(ptused, pid) ? 1 : 0;
          thisBlendGroupLen += 2 * incr; // {p0 p1}
          thisBlendGroups += incr;
        }

        // Save the results.
        blendGroupsView[zoneIndex] = thisBlendGroups;
        blendGroupLenView[zoneIndex] = thisBlendGroupLen;
        fragmentsView[zoneIndex] = thisFragments;
        fragmentsSizeView[zoneIndex] = thisFragmentsNumIds;

        // Sum up the sizes overall.
        fragment_sum += thisFragments;
        fragment_nids_sum += thisFragmentNumIds;
        blendGroups_sum += thisBlendGroups;
        blendGroupLen_sum += thisBlendGroupLen;
      });
    });


// TODO: change int to topoView::IndexType

    // ----------------------------------------------------------------------
    //
    // Stage 2: Do some scans to fill out blendOffset and blendGroupOffsets,
    //          which is where we fill in the real data.
    //
    // blendOffset : Starting offset for blending data like blendIds, blendCoeff.
    // blendGroupOffset : Starting offset for blendNames, blendGroupSizes.
    // fragmentOffsets : Where an element's fragments begin in the output.
    // ----------------------------------------------------------------------
    axom::Array<int> blendOffset(nzones, nzones, allocatorID);
    axom::Array<int> blendGroupOffsets(nzones, nzones, allocatorID);
    axom::Array<int> fragmentOffsets(nzones, nzones, allocatorID);
    axom::Array<int> fragmentSizeOffsets(nzones, nzones, allocatorID);

    auto blendOffsetView = blendOffset.view();
    auto blendGroupOffsetsView = blendGroupOffsets.view();
    auto fragmentOffsetsView = fragmentOffsets.view();

    // Make offsets via scan.
    RAJA::exclusive_scan<for_policy>(RAJA::make_span(blendGroupLenView, nzones),
                                     RAJA::make_span(blendOffsetView, nzones),
                                     RAJA::operators::plus<int>{});

    RAJA::exclusive_scan<for_policy>(RAJA::make_span(blendGroupsView, nzones),
                                     RAJA::make_span(blendGroupOffsetsView, nzones),
                                     RAJA::operators::plus<int>{});

    RAJA::exclusive_scan<for_policy>(RAJA::make_span(fragmentsView, nzones),
                                     RAJA::make_span(fragmentOffsetsView, nzones),
                                     RAJA::operators::plus<int>{});

    RAJA::exclusive_scan<for_policy>(RAJA::make_span(fragmentsSizeView, nzones),
                                     RAJA::make_span(fragmentSizeOffsetsView, nzones),
                                     RAJA::operators::plus<int>{});

    // ----------------------------------------------------------------------
    //
    // Stage 3: Iterate over the elements/cases again and fill in the blend
    //          groups that get produced: blendNames, blendGroupSizes,
    //          blendCoeff, blendIds. These are used to produce the new points.
    //
    //          NOTE: blendGroupStart is a scan of blendGroupSizes.
    //
    // ----------------------------------------------------------------------
    const auto blendGroupsSize = blendGroups_sum.get();
    const auto blendGroupLenSize = blendGroupLen_sum.get();

    axom::Array<KeyType> blendNames(blendGroupsSize, blendGroupsSize, allocatorID);
    axom::Array<int32> blendGroupSizes(blendGroupsSize, blendGroupsSize, allocatorID);
    axom::Array<int32> blendGroupStart(blendGroupsSize, blendGroupsSize, allocatorID);
    axom::Array<int32> blendIds(blendGroupLenSize, blendGroupLenSize, allocatorID);
    axom::Array<float> blendCoeff(blendGroupLenSize, blendGroupLenSize, allocatorID);

    auto blendNamesView = blendNames.view();
    auto blendGroupSizesView = blendGroupSizes.view();
    auto blendGroupStartView = blendGroupStart.view();
    auto blendIdsView = blendIds.view();
    auto blendCoeffView = blendCoeff.view();

    views::Node_to_ArrayView(n_clip_field, [&](auto clipFieldView)
    {
      topoView.template for_all_zones<ExecSpace>(AXOM_LAMBDA(auto zoneIndex, const auto &zone)
      {
        // Get the clip case for the current zone.
        const auto clipcase = clipCasesView[zoneIndex];

        // Iterate over the shapes in this clip case to determine the number of blend groups.
        const auto clipTableIndex = details::getClipTableIndex(zone.dimension(), zone.numberOfNodes());
        const auto &ctView = clipTableViews[clipTableIndex];

        const int64 ptused = pointsUsedView[zoneIndex];

        // Starting offset of where we store this element's blend groups.
        int32 bgStart = blendOffsetView[zoneIndex];
        int32 bgOffset = blendGroupOffsetsView[zoneIndex];

        auto it = ctView.begin(clipcase);
        const auto end = ctView.end(clipcase);
        for(; it != end; it++)
        {
          // Get the current shape in the clip case.
          const auto caseData = *it;

          if(caseData[0] == ST_PNT)
          {
            if(details::generatedPointIsSelected(caseData[2], selection))
            {
// TODO: put the ability to select on color back in... 0, 1, or both
              const int nIds = static_cast<int>(caseData[3]);
              const auto one_over_n = 1.f / static_cast<float>(nIds);
              const auto start = bgStart;

              for(int ni = 0; ni < nIds; ni++)
              {
                const auto ptid = caseData[4 + ni];

                // Add the point to the blend group.
                if(ptid <= P7)
                {
                  // corner point.
                  blendIdsView[bgStart] = zone.getId(ptid);
                  blendCoeffView[bgStart] = one_over_n;

                  bgStart++;
                }
                else if(ptid >= EA && ptid <= EL)
                {
                  // edge points are derived from 2 corner points. If
                  // those appear here then we're probably creating a
                  // face point. We can store the 2 corner points in place
                  // of the edge point (along with some blending coeff).
                  const auto edgeIndex = ptid - EA;
                  const auto edge = ZoneType::edges[edgeIndex];
                  const auto id0 = zone.getId(edge[0]);
                  const auto id1 = zone.getId(edge[1]);

                  // Figure out the blend for edge.
                  const float t = details::computeT(clipFieldView[id0], clipFieldView[id1]);
                
                  blendIdsView[bgStart]   = id0;
                  blendIdsView[bgStart+1] = id1;
                  blendCoeffView[bgStart] = one_over_n * (1. - t);
                  blendCoeffView[bgStart+1] = one_over_n * t;

                  bgStart += 2;
                }
              }

              // Store how many points make up this blend group. Note that the
              // size will not necessarily be equal to npts if edges were involved.
              int32 nblended = bgStart - blendGroupStartView[bgOffset];
              blendGroupSizesView[bgOffset] = nblended;

              // Store "name" of blend group.
              const auto blendName = make_name_n(blendIdsView, start, nblended);
              blendNamesView[bgOffset++] = blendName;
            }
          }
        }

        // Add blend group for each original point that was used.
        for(unsigned char pid = P0; pid <= P7; pid++)
        {
          if(axom::utilities::bitIsSet(ptused, pid))
          {
            // Store blend group info                
            blendIdsView[bgStart] = zone.getId(pid);
            blendCoeffView[bgStart] = 1.;

            // Store how many points make up this blend group.
            blendGroupSizesView[bgOffset] = 1;

            // Store where this blendGroup starts in the blendIds,blendCoeff.
            blendGroupStartView[bgOffset] = bgStart;

            // Store "name" of blend group.
            blendNamesView[bgOffset++] = make_name_1(zone.getId(pid));

            bgStart++;
          }
        }

        // Add blend group for each edge point that was used.
        for(unsigned char pid = EA; pid <= EL; pid++)
        {
          if(axom::utilities::bitIsSet(ptused, pid))
          {
            const auto edgeIndex = pid - EA;
            const auto edge = ZoneTypeEdges::edges[edgeIndex];
            const auto id0 = zone.getId(c[0]);
            const auto id1 = zone.getId(c[1]);

            // Figure out the blend for edge.
            const float t = details::computeT(clipFieldView[id0], clipfieldView[id1]);

            // Store blend group info                
            blendIdsView[bgStart]   = id0;
            blendIdsView[bgStart+1] = id1;
            blendCoeffView[bgStart] = (1. - t);
            blendCoeffView[bgStart+1] = t;

            // Store how many points make up this blend group.
            blendGroupSizesView[bgOffset] = 2;

            // Store where this blendGroup starts in the blendIds,blendCoeff.
            blendGroupStartView[bgOffset] = bgStart;

            // Store "name" of blend group.
            blendNamesView[bgOffset++] = make_name_2(id0, id1);

            bgStart += 2;
          }
        }
      });
    });

    // ----------------------------------------------------------------------
    //
    // Stage 4 - Make the blend groups unique based on their blendName.
    //
    // ----------------------------------------------------------------------
    // At this point, we have created the blend group data. We can now use the
    // blendNames to make unique blend groups. uNames contains a sorted list of
    // the unique blend group names while uIndices is their original index in
    // blendNames/blendGroupOffsets/blendGroupSizes.
    axom::Array<KeyType> uNames;
    axom::Array<axom::IndexType> uIndices;
    unique<ExecSpace>(blendNames, uNames, uIndices);

    auto uNamesView = uNames.view();
    auto uIndicesView = uIndices.view();

    // Bundle up blend data.
    axom::mir::utilities::blueprint::BlendData blend;
    blend.m_blendNamesView = uNames.view();
    blend.m_uniqueIdsView = uIndices.view();
    blend.m_blendGroupSizesView = blendGroupSizes.view();
    blend.m_blendGroupStartView = blendGroupStart.view();
    blend.m_blendIdsView = blendIds.view();
    blend.m_blendCoeffView = blendCoeff.view();


// I just had this thought like it would be nice to output the volumes... and a mapping of old2new zones
// I suppose with the mapping and the old/new meshes, I could compute the volume fractions.
// Should I pass in a outputTopo, outputCoordset, outputFields nodes?
// Then I don't have to care about the names - it's done outside this code.

    // ----------------------------------------------------------------------
    //
    // Stage 5 - Make new connectivity.
    //
    // ----------------------------------------------------------------------
    conduit::Node &n_topologies = outputMesh["topologies"];
    conduit::Node &n_topo = n_topologies[topoView.name()];

    const auto finalNumZones = fragment_sum.get();
    const auto finalConnSize = fragment_nids_sum.get();

    using ConnectivityType = typename topoView::IndexType;
    constexpr auto connTypeID = axom::mir::utilities::blueprint::cpp2conduit<ConnectivityType>::id;

    // Allocate connectivity.
    conduit::Node &n_conn = n_topo["elements/connectivity"];
    n_conn.set_allocator(allocatorID);
    n_conn.set(conduit::DataType(connTypeID, finalConnSize));
    auto connView = axom::ArrayView<ConnectivityType>(static_cast<ConnectivityType *>(n_conn.data_ptr()), finalConnSize);

    // Allocate shapes.
    conduit::Node &n_shapes = n_topo["elements/shapes"];
    n_shapes.set_allocator(allocatorID);
    n_shapes.set(conduit::DataType(connTypeID, finalNumZones));
    auto shapesView = axom::ArrayView<ConnectivityType>(static_cast<ConnectivityType *>(n_shapes.data_ptr()), finalNumZones);

    // Allocate sizes.
    conduit::Node &n_sizes = n_topo["elements/sizes"];
    n_sizes.set_allocator(allocatorID);
    n_sizes.set(conduit::DataType(connTypeID, finalNumZones));
    auto sizesView = axom::ArrayView<ConnectivityType>(static_cast<ConnectivityType *>(n_shapes.data_ptr()), finalNumZones);

    // Allocate offsets.
    conduit::Node &n_offsets = n_topo["elements/offsets"];
    n_offsets.set_allocator(allocatorID);
    n_offsets.set(conduit::DataType(connTypeID, finalNumZones));
    auto offsetsView = axom::ArrayView<ConnectivityType>(static_cast<ConnectivityType *>(n_offsets.data_ptr()), finalNumZones);


    const int32 uNames_len = uNames.size();

    RAJA::ReduceBitOr<reduce_policy, uint64> shapesUsed_reduce(0);

    // Here we get the node ids from the unique blend names so we might be able to skip the zone iteration since that is more expensive than a normal axom::for_all.
    // I'd have to store the clipTableIndex though.
#if 1
    AXOM_HOST_DEVICE
    constexpr int32 encode_case(int32 clipcase, int32 clipTableIndex, int32 nFragments)
    {
      return (clipcase & 0xff) |
             ((clipTableIndex & 0xff) << 8) |
             ((nFragments & 0xff) << 16);
    }
    AXOM_HOST_DEVICE
    constexpr void decode_case(int32 value, int32 &clipcase, int32 &clipTableIndex, int32 &nFragments)
    {
      clipcase = (value & 0xff);
      clipTableIndex = (value >> 8) & 0xff;
      nFragments = (value >> 16) & 0xff;
    }
    // Of course, I can't do an exclusive scan on the fragmentsView if it contains other stuff.

#endif
    topoView.template for_all_zones<ExecSpace>(AXOM_LAMBDA(auto zoneIndex, const auto &zone)
    {
      // If there are no fragments, return from lambda.
      if(fragmentsView[zoneIndex] == 0)
        return;

      // Seek to the start of the blend groups for this zone.
      const int32 bgStart = blendGroupOffsetsView[zoneIndex];

      // Go through the points in the order they would have been added as blend
      // groups, get their blendName, and then overall index of that blendName
      // in uNames, the unique list of new dof names. That will be their index
      // in the final points.
      const int64 ptused = pointsUsedView[zoneIndex];
      ConnectivityType point_2_new[N3 + 1];
      for(unsigned char pid = N0; pid <= N3; pid++)
      {
        if(axom::utilities::bitIsSet(ptused, pid))
        {
          const auto name = blendNamesView[bgStart++];
          point_2_new[pid] = bsearch(name, uNamesView);
        }
      }
      for(unsigned char pid = P0; pid <= P7; pid++)
      {
        if(axom::utilities::bitIsSet(ptused, pid))
        {
          const auto name = blendNamesView[bgStart++];
          point_2_new[pid] = bsearch(name, uNamesView);
        }
      }
      for(unsigned char pid = EA; pid <= EL; pid++)
      {
        if(axom::utilities::bitIsSet(ptused, pid))
        {
          const auto name = blendNamesView[bgStart++];
          point_2_new[pid] = bsearch(name, uNamesView);
        }
      }

      // This is where the output fragment connectivity start for this zone
      int outputIndex = fragmentSizeOffsetsView[zoneIndex];
      // This is where the output fragment sizes/shapes start for this zone.
      int sizeIndex = fragmentOffsetsView[zoneIndex];
      uint64 shapesUsed = 0;

      // Iterate over the selected fragments and emit connectivity for them.
      const auto clipcase = clipCasesView[zoneIndex];
      const auto clipTableIndex = details::getClipTableIndex(zone.dimension(), zone.numberOfNodes());
      const auto ctView = clipTableViews[clipTableIndex];
      it = ctView.begin(clipcase);
      for(; it != end; it++)
      {
        // Get the current shape in the clip case.
        const auto caseData = *it;

        if(caseData[0] != ST_PNT)
        {
          if(details::shapeIsSelected(caseData[1], selection))
          {
            // Output the nodes used in this zone.
            for(int i = 2; i < caseData.size(); i++)
              connView[outputIndex++] = point_2_new[caseData[i]];

            const auto nIdsInZone = caseData.size() - 2;
            sizesView[sizeIndex] = nIdsInZone;
            shapeView[sizeIndex] = zone.id();
            sizeIndex++;

            // Record which shape type was used. (zone.id() are already powers of 2)
            shapesUsed |= zone.id();
          }
        }
      }

      shapesUsed_reduce |= shapesUsed;
    });

    // Make offsets
    RAJA::exclusive_scan<for_policy>(RAJA::make_span(sizesView, nzones),
                                     RAJA::make_span(offsetsView, nzones),
                                     RAJA::operators::plus<int>{});

    // Add shape information to the connectivity.
    const auto shapesUsed = shapesUsed_reduce.get();
    const auto shapeMap = shapeMap_ValueName(shapesUsed);
    if(axom::utilities::countBits(shapesUsed) > 1)
    {
      n_topo["elements/shape"] = "mixed";
      conduit::Node &n_shape_map = n_topo["elements/shape_map"];
      for(const auto &[key, value] : shapeMap)
         n_shape_map[key] = value;
    }
    else
    {
      n_shapes.reset();
      n_topo["elements"].remove("shapes");

      n_topo["elements/shape"] = shapeMap.begin()->first;
    }

    //-----------------------------------------------------------------------
    // STAGE 6 - Make new coordset.
    //-----------------------------------------------------------------------
    
    make_fields(blend, inputMesh["fields"], outputMesh["fields"]);

    //-----------------------------------------------------------------------
    // STAGE 6 - Make new fields.
    //-----------------------------------------------------------------------
    if(inputMesh.has_path("fields"))
    {
       make_fields(blend, slice, topoName, inputMesh["fields"], outputMesh["fields"]);
    }

    //-----------------------------------------------------------------------
    // STAGE 6 - make originalElements (this will later be optional)
    //-----------------------------------------------------------------------
    conduit::Node &n_fields = outputMesh["fields"];
    if(n_fields.has_child("originalElements"))
    {
      // originalElements already exists. We need to map it forward.
      conduit::Node &n_orig = n_fields["originalElements"];
      conduit::Node &n_orig_values = n_orig["values"];
      views::IndexNode_to_ArrayView(n_orig_values, [&](auto origValuesView)
      {
        using value_type = typename origValuesView::value_type;
        conduit::Node n_values;
        n_values.set_allocator(allocatorID);
        n_values.set(conduit::DataType(n_orig_values.dtype().id(), finalNumZones));
        auto valuesView = axom::ArrayView<value_type>(static_cast<value_type *>(n_values.data_ptr()), finalNumZones);
        axom::for_all<ExecSpace>(nzones, AXOM_LAMBDA(auto zoneIndex)
        {
          int sizeIndex = fragmentOffsetsView[zoneIndex];
          int nFragments = fragmentsView[zoneIndex];
          for(int i = 0; i < nFragments; i++)
            valuesView[sizeIndex + i] = origValuesView[zoneIndex];
        });

        n_orig_values.move(n_values);
      });
    }
    else
    {
      // Make a new node and populate originalElement.
      conduit::Node &n_orig = n_fields["originalElements"];
      n_orig["association"] = "element";
      n_orig["topology"] = n_topo.name();
      conduit::Node &n_values = n_orig["values"];
      n_values.set_allocator(allocatorID);
      n_values.set(conduit::DataType(connTypeID, finalNumZones));
      auto valuesView = axom::ArrayView<ConnectivityType>(static_cast<ConnectivityType *>(n_values.data_ptr()), finalNumZones);
      axom::for_all<ExecSpace>(nzones, AXOM_LAMBDA(auto zoneIndex)
      {
        int sizeIndex = fragmentOffsetsView[zoneIndex];
        int nFragments = fragmentsView[zoneIndex];
        for(int i = 0; i < nFragments; i++)
          valuesView[sizeIndex + i] = zoneIndex;
      });
    }

//-----------------------------------------------------------------------------------
  } // end of execute

private:
  void make_coordset(const BlendData &blend, const CoordsetView &csview, conduit::Node &n_new_coordset) const
  {
    axom::mir::utilities::blueprint::CoordsetBlender<ExecSpace, CoordsetView> cb;
    b.execute(blend, csview, n_coordset, n_new_coordset);
  }

  void make_fields(const BlendData &blend, const SliceData &slice, const std::string &topoName, const conduit::Node &n_fields, conduit::Node &n_out_fields) const
  {
    for(conduit::index_t i = 0; i < n_fields.number_of_children(); i++)
    {
      const conduit::Node &n_field = n_fields[i];
      if(n_field["topology"].as_string() == topoName)
      {
        const std::string association = n_field["association"].as_string();
        if(association == "element")
        {
          axom::mir::utilities::blueprint::FieldSlicer<ExecSpace> s;
          s.execute(blend, n_field, n_out_fields[n_field.name()]);
        }
        else if(association == "vertex")
        {
          axom::mir::utilities::blueprint::FieldBlender<ExecSpace> b;
          b.execute(blend, n_field, n_out_fields[n_field.name()]);
        }
      }
    }
  }

  // NOTE: I probably want to make the clip tables be a class member so I can reuse the object without having to reload the clip tables each time.
};

} // end namespace clipping
} // end namespace mir
} // end namespace axom

#endif
