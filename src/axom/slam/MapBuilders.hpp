// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file MapBuilders.hpp
 *
 * \brief Construct maps with policies deduced from a set, stride and buffer.
 */

#pragma once

#include "axom/core/ArrayView.hpp"
#include "axom/slic.hpp"

#include "axom/slam/Concepts.hpp"
#include "axom/slam/Map.hpp"
#include "axom/slam/policies/IndirectionPolicies.hpp"
#include "axom/slam/policies/StridePolicies.hpp"

namespace axom::slam
{
namespace detail
{
template <typename SetType, typename PosType>
axom::IndexType map_storage_size(const SetType* set, PosType stride)
{
  using Position = typename SetType::PositionType;
  return detail::checkedMapStorageSize(set ? set->size() : Position {}, stride);
}
}  // namespace detail

/// \name Map construction helpers
/// \brief Construct maps that borrow their sets and value buffers.
/// The set and backing allocation must outlive the map. ArrayView objects are copied.
/// Runtime strides must be non-Boolean integral values
/// and be positive and representable in the set's position type.
/// The returned map always uses that position type.
/// Compile-time strides must be positive.
/// \{

/*!
 * \brief Make a strided SLAM map backed by ArrayView storage.
 *
 * \param set pointer to the map's set
 * \param stride number of components per set element
 * \param data view of the value buffer
 *
 * \pre `data.size() == set->size() * stride`, or zero for a null set.
 */
template <typename SetType, typename T, typename StrideType>
  requires std::integral<StrideType> && detail::SetPositionConvertible<SetType, StrideType>
auto make_map(const SetType* set, StrideType stride, axom::ArrayView<T> data)
{
  using PosType = typename SetType::PositionType;
  using Indirection = policies::ArrayViewIndirection<PosType, T>;
  using Stride = policies::RuntimeStride<PosType>;
  using MapType = Map<T, SetType, Indirection, Stride>;
  const auto canonicalStride = detail::checkedMapStride<PosType>(stride);
  return MapType(set, data, canonicalStride);
}

/*!
 * \brief Make a stride-one SLAM map backed by ArrayView storage.
 *
 * \pre `data.size() == set->size()`, or zero for a null set.
 */
template <typename SetType, typename T, typename ExplicitPosType = void>
  requires detail::OptionalSetPositionSame<SetType, ExplicitPosType>
auto make_map(const SetType* set, axom::ArrayView<T> data)
{
  using PosType = typename SetType::PositionType;
  using Indirection = policies::ArrayViewIndirection<PosType, T>;
  using Stride = policies::StrideOne<PosType>;
  using MapType = Map<T, SetType, Indirection, Stride>;
  return MapType(set, data);
}

/*!
 * \brief Make a strided SLAM map backed by a raw pointer buffer.
 *
 * This overload wraps the buffer as an ArrayView with length `set->size() * stride`
 * and returns an ArrayView-backed map. A null set gives a zero-length view.
 * The caller must provide a buffer with at least that many elements.
 */
template <typename SetType, typename T, typename StrideType>
  requires std::integral<StrideType> && detail::SetPositionConvertible<SetType, StrideType>
auto make_map(const SetType* set, StrideType stride, T* data)
{
  using PosType = typename SetType::PositionType;
  const auto canonicalStride = detail::checkedMapStride<PosType>(stride);
  const auto n = detail::map_storage_size(set, canonicalStride);
  return make_map(set, canonicalStride, axom::ArrayView<T>(data, n));
}

/*!
 * \brief Make a stride-one SLAM map backed by a raw pointer buffer.
 *
 * \pre The buffer holds at least `set->size()` values, or zero for a null set.
 */
template <typename SetType, typename T, typename ExplicitPosType = void>
  requires detail::OptionalSetPositionSame<SetType, ExplicitPosType>
auto make_map(const SetType* set, T* data)
{
  using PosType = typename SetType::PositionType;
  const auto n = detail::map_storage_size(set, PosType {1});
  return make_map(set, axom::ArrayView<T>(data, n));
}

/*!
 * \brief Make a compile-time strided SLAM map backed by ArrayView storage.
 *
 * \tparam STRIDE number of values per set element
 *
 * \pre `data.size() == set->size() * STRIDE`, or zero for a null set.
 */
template <int STRIDE, typename SetType, typename T, typename ExplicitPosType = void>
  requires detail::PositiveStaticStrideFor<STRIDE, SetType> &&
  detail::OptionalSetPositionSame<SetType, ExplicitPosType>
auto make_map_ct(const SetType* set, axom::ArrayView<T> data)
{
  using PosType = typename SetType::PositionType;
  using Indirection = policies::ArrayViewIndirection<PosType, T>;
  using Stride = policies::CompileTimeStride<PosType, static_cast<PosType>(STRIDE)>;
  using MapType = Map<T, SetType, Indirection, Stride>;
  return MapType(set, data);
}

/*!
 * \brief Make a compile-time strided SLAM map backed by a raw pointer buffer.
 *
 * This overload wraps the buffer as an ArrayView with length `set->size() * STRIDE`
 * and returns an ArrayView-backed map. A null set gives a zero-length view.
 * The caller must provide a buffer with at least that many elements.
 */
template <int STRIDE, typename SetType, typename T, typename ExplicitPosType = void>
  requires detail::PositiveStaticStrideFor<STRIDE, SetType> &&
  detail::OptionalSetPositionSame<SetType, ExplicitPosType>
auto make_map_ct(const SetType* set, T* data)
{
  using PosType = typename SetType::PositionType;
  const PosType stride = static_cast<PosType>(STRIDE);
  const auto n = detail::map_storage_size(set, stride);
  return make_map_ct<STRIDE>(set, axom::ArrayView<T>(data, n));
}

/// \}

}  // namespace axom::slam
