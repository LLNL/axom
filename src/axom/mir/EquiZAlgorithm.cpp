// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

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

void EquiZAlgorithm::execute(const conduit::Node &topo,
                             const conduit::Node &coordset,
                             const conduit::Node &matset,
                             const conduit::Node &options,
                             conduit::Node &new_topo,
                             conduit::Node &new_coordset,
                             conduit::Node &new_matset)
{
#if defined (AXOM_USE_RAJA) && defined (AXOM_USE_UMPIRE)
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
}

#if 0
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
struct zone2index { static constexpr size_t st_index = 0; };

template <>
struct zone2index<TriShape> { static constexpr size_t st_index = ST_TRI; };

template <>
struct zone2index<QuadShape> { static constexpr size_t st_index = ST_QUA; };

template <>
struct zone2index<TetShape> { static constexpr size_t st_index = ST_TET; };

template <>
struct zone2index<PyramidShape> { static constexpr size_t st_index = ST_PYR; };

template <>
struct zone2index<WedgeShape> { static constexpr size_t st_index = ST_WDG; };

template <>
struct zone2index<HexShape> { static constexpr size_t st_index = ST_HEX; };

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

//------------------------------------------------------------------------------
template <typename T>
struct cpp2conduit { static constexpr conduit::index_t type = conduit::DataType::EMPTY_ID; };

template <>
struct cpp2conduit<conduit::int8> { static constexpr conduit::index_t type = conduit::DataType::INT8_ID; };

template <>
struct cpp2conduit<conduit::int16> { static constexpr conduit::index_t type = conduit::DataType::INT16_ID; };

template <>
struct cpp2conduit<conduit::int32> { static constexpr conduit::index_t type = conduit::DataType::INT32_ID; };

template <>
struct cpp2conduit<conduit::int64> { static constexpr conduit::index_t type = conduit::DataType::INT64_ID; };

template <>
struct cpp2conduit<conduit::uint8> { static constexpr conduit::index_t type = conduit::DataType::UINT8_ID; };

template <>
struct cpp2conduit<conduit::uint16> { static constexpr conduit::index_t type = conduit::DataType::UINT16_ID; };

template <>
struct cpp2conduit<conduit::uint32> { static constexpr conduit::index_t type = conduit::DataType::UINT32_ID; };

template <>
struct cpp2conduit<conduit::uint64> { static constexpr conduit::index_t type = conduit::DataType::UINT64_ID; };

template <>
struct cpp2conduit<conduit::float32> { static constexpr conduit::index_t type = conduit::DataType::FLOAT32_ID; };

template <>
struct cpp2conduit<conduit::float64> { static constexpr conduit::index_t type = conduit::DataType::FLOAT64_ID; };


//------------------------------------------------------------------------------
AXOM_HOST_DEVICE
constexpr int getClipTableIndex(int dimension, int nnodes)
{
  return (dimension == 2) ? ((nnodes == 3) ? 0 : 1) : (nnodes - 2);
}

template <typename FlagType, typename BitType>
AXOM_HOST_DEVICE
constexpr bool bitIsSet(FlagType flags, BitType bit)
{
  return (flags & (1 << bit)) > 0;
}

template <typename FlagType, typename BitType>
AXOM_HOST_DEVICE
constexpr void setBit(FlagType &flags, BitType bit)
{
  flags |= (1 << bit);
}

template <typename IntType>
AXOM_HOST_DEVICE
int countBits(IntType value)
{
  constexpr int n = sizeof(IntType) * 8;
  int count = 0, mask = 1;
  for(int i = 0; i < n; i++)
    count += ((value & mask) > 0) ? 1 : 0;
  return count;
}

AXOM_HOST_DEVICE
constexpr bool color0Selected(int selection)
{
  return bitIsSet(selection, 0);
}

AXOM_HOST_DEVICE
constexpr bool color1Selected(int selection)
{
  return bitIsSet(selection, 1);
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

//-------------------------------------------------------------------------
template <typename T>
AXOM_HOST_DEVICE
int32 bsearch(T value, const axom::ArrayView<T> &view)
{
  int32 index = -1;
  int32 left = 0;
  int32 right = view.size() - 1;
  while(left <= right)
  {
    int32 m = (left + right) / 2;
    if(view[m] < value)
      left = m + 1;
    else if(view[m] > value)
      right = m - 1;
    else
    {
      index = m;
      break;
    }
  }

  return index;
}

//------------------------------------------------------------------------------
/// Based on a Jenkins hash, modified to include length and hash forwards and
/// backwards to make int64 rather than int32.
AXOM_HOST_DEVICE
uint64 hash_bytes(const uint8 *data, uint32 length)
{
  uint32 hash = 0;

  // Build the length into the hash.
  const auto ldata = reinterpret_cast<const uint8 *>(&length);
  for(int e = 0; e < 4; e++)
  {
    hash += ldata[e];
    hash += hash << 10;
    hash ^= hash >> 6;
  }

  uint32 hashr = hash;
  for(uint32 i = 0; i < length; i++)
  {
    hash += data[i];
    hash += hash << 10;
    hash ^= hash >> 6;

    hashr += data[length - 1 - i];
    hashr += hashr << 10;
    hashr ^= hashr >> 6;
  }
  hash += hash << 3;
  hash ^= hash >> 11;
  hash += hash << 15;

  hashr += hashr << 3;
  hashr ^= hashr >> 11;
  hashr += hashr << 15;

  return (static_cast<uint64>(hash) << 32) | hashr;
}

//------------------------------------------------------------------------------
template <typename ValueType>
uint64
AXOM_HOST_DEVICE
make_name_1(ValueType id)
{
  return hash_bytes(reinterpret_cast<uint8*>(&id), sizeof(ValueType));
};

//------------------------------------------------------------------------------
template <typename ValueType>
AXOM_HOST_DEVICE
uint64
make_name_2(ValueType id0, ValueType id1)
{
  ValueType data[2] = {id0, id1};
  if(id1 < id0)
  {
    data[0] = id1;
    data[1] = id0;
  }
  return hash_bytes(reinterpret_cast<uint8*>(data), 2 * sizeof(ValueType));
};

//------------------------------------------------------------------------------
template <typename ValueType, typename IndexType>
AXOM_HOST_DEVICE
uint32
sort_values(ValueType *v, IndexType n)
{
  for(IndexType i = 0; i < n-1; i++)
  {
    const IndexType m = n - i - 1;
    for(IndexType j = 0; j < m; j++)
    {
      if(v[j] > v[j+1])
      {
        axom::utilities::swap(v[j], v[j+1]);
      }
    }
  }
}

//-------------------------------------------------------------------------
template <typename ViewType>
uint64
AXOM_HOST_DEVICE
make_name_n(const ViewType &view, int32 start, int32 n)
{
  using value_type = typename ViewType::value_type;

  if(n == 2)
    return make_name_2(view[start], view[start + 1]);

  value_type v[14]={0,0,0,0,0,0,0,0,0,0,0,0,0,0}; // pick largest number of blends.
  for(int32 i = 0; i < n; i++)
  {
    v[i] = view[start + i];
  }
  sort_values(v, n);

  return hash_bytes(reinterpret_cast<uint8*>(v), n * sizeof(value_type));
}

//------------------------------------------------------------------------------
// I might want to do the real work here in a more targeted class. This way, a host code could choose to instantiate a single
template <typename ExecSpace, typename TopologyView, typename CoordsetView>
class ClipField
{
public:
  void execute(const TopologyView &topoView, // I'd rather just pass the views to the method.
               const CoordsetView &coordsetView,
               const conduit::Node &clipField,
               conduit::Node &outputMesh)
  {
    using int32 = conduit::int32;
    using uint32 = conduit::uint32;
    using int64 = conduit::int64;
    using uint64 = conduit::uint64;

    using KeyType = uint64;
    using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
    using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;
    const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    // Load clip table data and make views.
    ClipTableManager<ExecSpace> clipTables;
    clipTables.load(topoView.dimension());
    axom::StackArray<ClipTableManager<ExecSpace>::View, 6> clipTableViews;
    if(topoView.dimension() == 2)
    {
      clipTableViews[indexOfZoneType<TriShape>()] = clipTables[ST_TRI].view();
      clipTableViews[indexOfZoneType<QuadShape>()] = clipTables[ST_QUA].view();
    }
    if(topoView.dimension() == 2)
    {
      clipTableViews[indexOfZoneType<TetShape>()] = clipTables[ST_TET].view();
      clipTableViews[indexOfZoneType<PyramidShape>()] = clipTables[ST_PYR].view();
      clipTableViews[indexOfZoneType<WedgeShape>()] = clipTables[ST_WDG].view();
      clipTableViews[indexOfZoneType<HexShape>()] = clipTables[ST_HEX].view();
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
    axom::Array<int> blendGroups(nzones, nzones, allocatorID);    // Number of blend groups in a zone.
    axom::Array<int> blendGroupsLen(nzones, nzones, allocatorID); // Length of the blend groups in a zone.
    axom::Array<int> fragments(nzones, nzones, allocatorID);      // The number of fragments (child zones) produced for a zone.
    axom::Array<int> fragmentsSize(nzones, nzones, allocatorID);  // The total number of points produced for all fragments in a zone.
    auto clipCasesView = clipCases.view();
    auto blendGroupsview = blendGroups.view();
    auto blendGroupsLenView = blendGroupsLen.view();
    auto fragmentsView = fragments.view();
    auto fragmentsSizeView = fragmentsSize.view();

    views::Node_to_ArrayView(clipField, [&](auto clipFieldView)
    {
      topoView.template for_all_zones<ExecSpace>(AXOM_LAMBDA(auto zoneIndex, const auto &zone)
      {
        // Get the clip case for the current zone.
        const auto clipcase = clip_case(zone, clipFieldView);
        clipCasesView[zoneIndex] = clipcase;

        // Iterate over the shapes in this clip case to determine the number of blend groups.
        const auto clipTableIndex = getClipTableIndex(zone.dimension(), zone.numberOfNodes());
        const auto &ctView = clipTableViews[clipTableIndex];
        const size_t nShapes = ctView.shapesForCase(clipCases[zoneIndex]);

        int thisBlendGroups = 0;    // The number of blend groups produced in this case.
        int thisBlendGroupLen = 0;  // The total length of the blend groups.
        int thisFragments = 0;      // The number of zone fragments produced in this case.
        int thisFragmentsNumIds = 0;// The number of points used to make all the fragment zones.
        int64 ptused = 0;           // A bitset indicating which ST_XX nodes are used.

        for(size_t si = 0; si < nShapes; si++)
        {
          // Get the si'th shape in the clip case.
          const auto caseData = ctView.getShape(clipcase, si);

          if(caseData[0] == ST_PNT)
          {
            if(generatedPointIsSelected(caseData[2], selection))
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
            }
          }
          else
          {
            if(shapeShapeSelected(caseData[1], selection))
            {
              thisFragments++;
              const int nIdsThisShape = caseData.size() - 2;
              thisFragmentNumIds += nIdsThisShape;

              // Mark which points were used in this cell.
              for(size_t i = 2; i < caseData.size(); i++)
              {
                setBit(ptused, caseData[i]);
              }
            }
          }
        }

        // Count which points in the original cell are used.
        for(unsigned char pid = P0; pid <= P7; pid++)
        {
          const int incr = bitIsSet(ptused, pid) ? 1 : 0;
          thisBlendGroupLen += incr; // {p0}
          thisBlendGroups += incr;
        }

        // Count edges that are used.
        for(unsigned char pid = EA; pid <= EL; pid++)
        {
          const int incr = bitIsSet(ptused, pid) ? 1 : 0;
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

    views::Node_to_ArrayView(clipField, [&](auto clipFieldView)
    {
      topoView.template for_all_zones<ExecSpace>(AXOM_LAMBDA(auto zoneIndex, const auto &zone)
      {
        // Get the clip case for the current zone.
        const auto clipcase = clipCasesView[zoneIndex];

        // Iterate over the shapes in this clip case to determine the number of blend groups.
        const auto clipTableIndex = getClipTableIndex(zone.dimension(), zone.numberOfNodes());
        const auto &ctView = clipTableViews[clipTableIndex];
        const size_t nShapes = ctView.shapesForCase(clipCases[zoneIndex]);

        int64 ptused = 0;

        // Starting offset of where we store this element's blend groups.
        int32 bgStart = blendOffsetView[zoneIndex];
        int32 bgOffset = blendGroupOffsetsView[zoneIndex];

        for(size_t si = 0; si < nShapes; si++)
        {
          // Get the si'th shape in the clip case.
          const auto caseData = ctView.getShape(clipcase, si);

          if(caseData[0] == ST_PNT)
          {
            if(generatedPointIsSelected(caseData[2], selection))
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
                  const float t = computeT(clipFieldView[id0], clipFieldView[id1]);
                
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
          else
          {
            if(shapeShapeSelected(caseData[1], selection))
            {
              // Mark which points were used in this zone.
              for(size_t i = 2; i < caseData.size(); i++)
              {
                setBit(ptused, caseData[i]);
              }
            }
          }
        }

        // Add blend group for each original point that was used.
        for(unsigned char pid = P0; pid <= P7; pid++)
        {
          if(bitIsSet(ptused, pid))
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
          if(bitIsSet(ptused, pid))
          {
            const auto edgeIndex = pid - EA;
            const auto edge = ZoneTypeEdges::edges[edgeIndex];
            const auto id0 = zone.getId(c[0]);
            const auto id1 = zone.getId(c[1]);

            // Figure out the blend for edge.
            const float t = computeT(clipFieldView[id0], clipfieldView[id1]);

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
    const auto connTypeID = cpp2conduit<ConnectivityType>::type;

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

    topoView.template for_all_zones<ExecSpace>(AXOM_LAMBDA(auto zoneIndex, const auto &zone)
    {
      // If there are no fragments, return from lambda.
      if(fragmentsView[zoneIndex] == 0)
        return;

      const auto clipcase = clip_case(zone, clipFieldView);

      const auto clipTableIndex = getClipTableIndex(zone.dimension(), zone.numberOfNodes());
      const auto ctView = clipTableViews[clipTableIndex];
      const size_t nShapes = ctView.shapesForCase(clipCases[zoneIndex]);

      int64 ptused = 0;

      // Iterate over the tet fragments to see which points they use.
      // The points correspond to blend groups that we have made in
      // the final output.
      for(size_t si = 0; si < nShapes; si++)
      {
        // Get the si'th shape in the clip case.
        const auto caseData = ctView.getShape(clipcase, si);

        if(caseData[0] == ST_PNT)
        {
          // Count the point as used since a blend group would have
          // been emitted for it under these conditions. This makes
          // sure that we get the ordering for point_2_newdof right
          // when a set of shapes happens to not use the blended point.
          // That, of course, means that the lut needs to be fixed a bit.
          if(generatedPointIsSelected(caseData[2], selection))
          {
            setBit(ptused, N0 + caseData[1]);
          }
        }
        else
        {
          if(shapeShapeSelected(caseData[1], selection))
          {
            for(int i = 2; i < caseData.size(); i++)
              setBit(ptused, i);
          }
        }
      }

      // Seek to the start of the blend groups for this zone.
      const int32 bgStart = blendGroupOffsetsView[zoneIndex];

      // Go through the points in the order they would have been added as blend
      // groups, get their blendName, and then overall index of that blendName
      // in uNames, the unique list of new dof names. That will be their index
      // in the final points.
      ConnectivityType point_2_new[N3 + 1];
      for(unsigned char pid = N0; pid <= N3; pid++)
      {
        if(bitIsSet(ptused, pid))
        {
          const auto name = blendNamesView[bgStart++];
          point_2_new[pid] = bsearch(name, uNamesView);
        }
      }
      for(unsigned char pid = P0; pid <= P7; pid++)
      {
        if(bitIsSet(ptused, pid))
        {
          const auto name = blendNamesView[bgStart++];
          point_2_new[pid] = bsearch(name, uNamesView);
        }
      }
      for(unsigned char pid = EA; pid <= EL; pid++)
      {
        if(bitIsSet(ptused, pid))
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
      for(size_t si = 0; si < nShapes; si++)
      {
        // Get the si'th shape in the clip case.
        const auto caseData = ctView.getShape(clipcase, si);

        if(caseData[0] != ST_PNT)
        {
          if(shapeIsSelected(caseData[1], selection))
          {
            // Output the nodes used in this zone.
            for(int i = 2; i < caseData.size(); i++)
              connView[outputIndex++] = point_2_new[caseData[i]];

            const auto nIdsInZone = caseData.size() - 2;
            sizesView[sizeIndex] = nIdsInZone;
            shapeView[sizeIndex] = zone.id();
            setBit(shapesUsed, zone.id());
            sizeIndex++;
          }
        }
      }

      shapesUsed_reduce |= shapesUsed;
    });

    // Make offsets
    RAJA::exclusive_scan<for_policy>(RAJA::make_span(sizesView, nzones),
                                     RAJA::make_span(offsetsView, nzones),
                                     RAJA::operators::plus<int>{});

    // Add shape information
    const auto shapesUsed = shapesUsed_reduce.get();
    const auto shapeMap = shapeMap_ValueName(shapesUsed);
    if(countBits(shapesUsed) > 1)
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

    //-----------------------------------------------------------------------
    // STAGE 7 - Make new fields.
    //-----------------------------------------------------------------------

//-----------------------------------------------------------------------------------
  } // end of execute
};

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

template <typename ExecSpace, typename TopologyView, typename CoordsetView>
void
EquiZAlgorithm<ExecSpace, TopologyView, CoordsetView>::execute(
  const TopologyView &topoView,
  const CoordsetView &coordsetView,
  const conduit::Node &options)
{
}
//----------------------------------------------------------------------------------------
#endif

#if 0
/// Provide overloads of initialize_topology_view for every possible topology type.
template <typename IndexT>
void initialize_topology_view(const std::conduit::Node &topo, StructuredTopologyView<StridedStructuredIndexing<IndexT, 3>> &topoView)
{
  // Initialize the view from the Conduit node.
}

template <typename ExecSpace, typename TopologyView, typename CoordsetView, typename MatsetView>
void
EquiZAlgorithm<ExecSpace, TopologyView, CoordsetView, MatsetView>::execute(
  const TopologyView &topoView,
  const CoordsetView &coordsetView,
  const MatsetView   &matsetView,
  const conduit::Node &options)
{
  if(options.has_path("zones"))
  {
    const conduit::Node &n_zones = options.fetch_existing("zones");

/// NOTE: since each inner dispatch could be a lot of code, should I just make a zones array for the case where zones is not provided?

    // Operate on a list of zones.
    views::IndexNode_to_ArrayView(n_zones, [&](auto zonesView)
    {
      MaterialInformation matinfo = materials(matset);
      for(const auto &mat : matinfo)
      {
        const auto matID = mat.number;

        // Going this way, the relation builder should take in the topoView to do its work.
        axom::mir::utilities::NodeToZoneRelationBuilder<ExecSpace, TopologyView> nz;
        nz.execute(topoView);

        const auto relZonesView = nz.zones().view();
        const auto relSizesView = nz.sizes().view();
        const auto relOffsetsView = nz.offsets().view();

        // Create the clipping tables for the topo dimension.
        axom::mir::clipping::ClipTableManager<ExecSpace> clipManager;
        clipManager.load(topoView.dimension());

        // We need to get views for the various shape types.
        axom::StackArray<ClipTableView, ST_MAX>

        topoView. template for_selected_zones<ExecSpace>(zonesView, AXOM_LAMBDA(auto zoneIndex, const auto &zone)
        {          
                
        });
      }
    });
  }
  else
  {
    // Operate on all zones.
    
    topoView. template for_all_zones<ExecSpace>(AXOM_LAMBDA(auto zoneIndex, const auto &zone)
    {          

    });
  }
}


// I'm starting to want to just run EquiZ on a certain execution space with a certain input mesh rather than all of them...
//
// EquiZAlgorithm<cuda_exec, StructuredTopologyView<StridedStructuredIndexing<int, 3>>> mir;
// mir.execute(topo, coordset, matset, options, new_topo, new_coordset, new_matset);
//
// But what about fields??
//
// mir.execute(mesh, options, output);
// mir.execute(mesh, options, output["topologies/topo"], output["coordset/newcoords"], output["matsets/newmatset"]);
// mir.execute(mesh, options, output["topologies/topo"], output["coordset/newcoords"], output["matsets/newmatset"], output["fields"]);
#endif

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
#else
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
#if 1
                // We need to get views for the various shape types.
                using TableView = typename axom::mir::clipping::ClipTableManager<ExecSpace>::Table::View; 
                axom::StackArray<TableView, ST_MAX> tables;
                tables[ST_TET] = clipManager[ST_TET].view();
#endif

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

} // namespace mir
} // namespace axom
