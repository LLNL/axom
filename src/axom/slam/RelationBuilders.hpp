// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file RelationBuilders.hpp
 *
 * \brief Construct static relations, deducing set types and policies from their buffers.
 *
 * For example, bind existing begin offsets and to-set positions with:
 *
 * \code
 *   auto r = slam::make_variable_relation(&from, &to, offsets, indices);
 * \endcode
 *
 * The helpers borrow the sets and connectivity storage. They check buffer sizes
 * without reading connectivity values. Use StaticRelation::isValid() on accessible
 * storage to check the begin offsets and to-set positions.
 */

#pragma once

#include "axom/slam/Concepts.hpp"
#include "axom/slam/StaticRelation.hpp"
#include "axom/slam/policies/CardinalityPolicies.hpp"
#include "axom/slam/policies/IndirectionPolicies.hpp"

#include "axom/core/ArrayView.hpp"
#include "axom/slic.hpp"

#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

namespace axom::slam
{
namespace detail
{
template <typename FromSet, typename ToSet, typename ElemType>
concept RelationIndexBufferTypes =
  SetLike<FromSet> && SetLike<ToSet> && SetPositionSame<ToSet, ElemType>;

template <typename FromSet, typename ToSet>
using RelationFlatPosition =
  default_flat_position_t<typename FromSet::PositionType, typename ToSet::PositionType>;

template <typename FromSet, typename ToSet, typename ExplicitPosition>
using SelectedRelationFlatPosition =
  std::conditional_t<std::same_as<model_t<ExplicitPosition>, void>,
                     RelationFlatPosition<FromSet, ToSet>,
                     model_t<ExplicitPosition>>;

template <typename FromSet, typename ToSet, typename Position>
concept RelationFlatPositionFor = SetLike<FromSet> && SetLike<ToSet> &&
  PositionCanRepresent<Position, typename model_t<FromSet>::PositionType>;

template <typename FromSet, typename ToSet, typename ExplicitPosition>
concept OptionalRelationFlatPositionFor = SetLike<FromSet> && SetLike<ToSet> &&
  RelationFlatPositionFor<FromSet, ToSet, SelectedRelationFlatPosition<FromSet, ToSet, ExplicitPosition>>;

template <typename FromSet, typename ToSet, typename Value>
concept RelationFlatPositionConstructible =
  SetLike<FromSet> && SetLike<ToSet> && PositionValueLike<Value> &&
  std::constructible_from<RelationFlatPosition<FromSet, ToSet>, model_t<Value>>;

template <typename FromSet, typename ToSet, typename ExplicitPosition, typename Value>
concept SelectedRelationFlatPositionConstructible =
  OptionalRelationFlatPositionFor<FromSet, ToSet, ExplicitPosition> && PositionValueLike<Value> &&
  std::constructible_from<SelectedRelationFlatPosition<FromSet, ToSet, ExplicitPosition>, model_t<Value>>;

template <typename FromSet, typename ToSet, typename PosType, typename ElemType>
concept VariableRelationBufferTypes = RelationIndexBufferTypes<FromSet, ToSet, ElemType> &&
  RelationFlatPositionFor<FromSet, ToSet, PosType>;

template <typename PosType, typename Value>
inline PosType checked_relation_position(Value value, const char* description)
{
  using ValueType = std::remove_cvref_t<Value>;
  if constexpr(std::integral<PosType> && std::integral<ValueType>)
  {
    if constexpr(std::signed_integral<ValueType>)
    {
      SLIC_ERROR_IF(
        value < ValueType {},
        "SLAM relation " << description << " must be nonnegative; received " << value << ".");
    }
    SLIC_ERROR_IF(!std::in_range<PosType>(value),
                  "SLAM relation " << description << " (" << value
                                   << ") is not representable by the relation position type.");
  }
  return static_cast<PosType>(value);
}

template <typename PosType, typename Value>
inline PosType checked_relation_stride(Value value)
{
  using ValueType = std::remove_cvref_t<Value>;
  if constexpr(std::integral<ValueType>)
  {
    SLIC_ERROR_IF(
      value <= ValueType {},
      "slam::make_constant_relation -- runtime stride must be positive; received " << value << ".");
  }

  const PosType stride = checked_relation_position<PosType>(value, "runtime stride");
  SLIC_ERROR_IF(stride <= PosType {},
                "slam::make_constant_relation -- runtime stride must be positive.");
  return stride;
}

/*!
 * \brief Return the number of from-set elements in the relation position type.
 *
 * A null set has size zero. Non-null sizes must be nonnegative and representable.
 */
template <typename PosType, typename FromSet>
inline PosType relation_from_size(const FromSet* fromSet)
{
  return fromSet ? checked_relation_position<PosType>(fromSet->size(), "from-set size") : PosType {};
}

/*!
 * \brief Check that the begins array backing a variable-cardinality relation is correctly sized.
 *
 * A variable relation needs one begin offset per from-set element plus a final
 * offset. The buffer must contain exactly `fromSet->size() + 1` entries.
 * This check reads no buffer values, so the connectivity may reside on a device.
 * StaticRelation::isValid() checks the offset values against the index count.
 */
template <typename PosType, typename FromSet, typename SizeType>
inline PosType check_variable_relation_size(const FromSet* fromSet, SizeType beginsSize)
{
  const PosType fromSize = relation_from_size<PosType>(fromSet);
  const PosType canonicalBeginsSize =
    checked_relation_position<PosType>(beginsSize, "begins-buffer size");
  SLIC_ERROR_IF(fromSize == std::numeric_limits<PosType>::max(),
                "slam::make_variable_relation -- the from-set size cannot be represented together "
                "with the required terminal begin offset.");

  const PosType expected = fromSize + PosType {1};
  SLIC_ERROR_IF(canonicalBeginsSize != expected,
                "slam::make_variable_relation -- begins has "
                  << canonicalBeginsSize << " entries, but the from-set (size " << fromSize
                  << ") requires exactly " << expected
                  << " (one begin offset per element plus a terminal).");
  return canonicalBeginsSize;
}

/*!
 * \brief Check that the indices array backing a constant-cardinality
 * relation (with stride \a stride) is correctly sized.
 *
 * A constant-cardinality relation indexes through `pos * stride`, so \a indices must
 * contain exactly `fromSet->size() * stride` entries.
 * This check inspects only sizes, so it is safe for device-resident storage.
 */
template <typename FromSet, typename PosType>
inline void check_constant_relation_size(const FromSet* fromSet, PosType stride, PosType indicesSize)
{
  SLIC_ERROR_IF(
    stride <= PosType {},
    "slam::make_constant_relation -- runtime stride must be positive; received " << stride << ".");

  const PosType fromSize = relation_from_size<PosType>(fromSet);
  SLIC_ERROR_IF(fromSize != PosType {} && stride > std::numeric_limits<PosType>::max() / fromSize,
                "slam::make_constant_relation -- from-set size "
                  << fromSize << " and stride " << stride
                  << " overflow the relation position type when multiplied.");

  const PosType expected = fromSize * stride;
  SLIC_ERROR_IF(indicesSize != expected,
                "slam::make_constant_relation -- indices has "
                  << indicesSize << " entries, but the from-set (size " << fromSize
                  << ") with stride " << stride << " requires exactly " << expected << ".");
}
}  // namespace detail

/// \name Relation construction helpers
/// \brief Construct relations whose entries use the to-set position type.
/// Variable relations use their begins-buffer element type for flattened storage.
/// Constant relations use the common position type of the from-set and to-set,
/// unless an explicit flat type is supplied.
/// Runtime sizes and strides must be non-Boolean integral values representable
/// by the flat position type. Sizes are nonnegative and strides are positive.
/// Compile-time strides must be positive.
/// The sets and storage must outlive the relation. Array and vector objects are
/// borrowed by pointer. ArrayView objects are copied and borrow their allocations.
/// \{

/*!
 * \brief Make a static, variable-cardinality relation
 *  from \a fromSet to \a toSet, backed by std::vector storage for its begins and indices.
 *
 * The from/to set types are deduced from the pointers.
 * The begin offsets are expressed in the relation's flat position type.
 * The entries in \a indices are positions in the to-set. Both buffers must
 * outlive the relation. The relation uses STL-vector indirection.
 *
 * \param fromSet pointer to the from-set
 * \param toSet pointer to the to-set
 * \param begins begin offsets for the from-set positions, followed by a final offset
 * \param indices to-set positions in flat storage
 * \return a StaticRelation with VariableCardinality and STLVector indirection
 *
 * \pre begins.size() == fromSet->size() + 1
 */
template <typename FromSet, typename ToSet, typename PosType, typename ElemType>
  requires detail::VariableRelationBufferTypes<FromSet, ToSet, PosType, ElemType>
auto make_variable_relation(FromSet* fromSet,
                            ToSet* toSet,
                            std::vector<PosType>& begins,
                            std::vector<ElemType>& indices)
{
  using BeginsIndirection = policies::STLVectorIndirection<PosType, PosType>;
  using IndicesIndirection = policies::STLVectorIndirection<PosType, ElemType>;
  using Cardinality = policies::VariableCardinality<PosType, BeginsIndirection>;
  using RelationType =
    StaticRelation<PosType, ElemType, Cardinality, IndicesIndirection, FromSet, ToSet>;
  using Builder = typename RelationType::RelationBuilder;

  const PosType beginsSize = detail::check_variable_relation_size<PosType>(fromSet, begins.size());
  const PosType indicesSize =
    detail::checked_relation_position<PosType>(indices.size(), "indices-buffer size");
  return RelationType(
    Builder()
      .fromSet(fromSet)
      .toSet(toSet)
      .begins(typename Builder::BeginsSetBuilder().size(beginsSize).data(&begins))
      .indices(typename Builder::IndicesSetBuilder().size(indicesSize).data(&indices)));
}

/// \brief Reference overload for make_variable_relation (std::vector-backed).
template <typename FromSet, typename ToSet, typename PosType, typename ElemType>
  requires detail::VariableRelationBufferTypes<FromSet, ToSet, PosType, ElemType>
auto make_variable_relation(FromSet& fromSet,
                            ToSet& toSet,
                            std::vector<PosType>& begins,
                            std::vector<ElemType>& indices)
{
  return make_variable_relation(&fromSet, &toSet, begins, indices);
}

/*!
 * \brief Make a static, variable-cardinality relation backed by C array storage.
 *
 * \param fromSet pointer to the from-set
 * \param toSet pointer to the to-set
 * \param begins pointer to begin offsets
 * \param beginsSize number of begin offsets
 * \param indices pointer to to-set positions in flat storage
 * \param indicesSize number of indices
 * \pre beginsSize == fromSet->size() + 1, or one for a null from-set.
 */
template <typename FromSet, typename ToSet, typename PosType, typename ElemType>
  requires detail::VariableRelationBufferTypes<FromSet, ToSet, PosType, ElemType>
auto make_variable_relation(FromSet* fromSet,
                            ToSet* toSet,
                            PosType* begins,
                            PosType beginsSize,
                            ElemType* indices,
                            PosType indicesSize)
{
  using BeginsIndirection = policies::CArrayIndirection<PosType, PosType>;
  using IndicesIndirection = policies::CArrayIndirection<PosType, ElemType>;
  using Cardinality = policies::VariableCardinality<PosType, BeginsIndirection>;
  using RelationType =
    StaticRelation<PosType, ElemType, Cardinality, IndicesIndirection, FromSet, ToSet>;
  using Builder = typename RelationType::RelationBuilder;

  const PosType canonicalBeginsSize =
    detail::check_variable_relation_size<PosType>(fromSet, beginsSize);
  const PosType canonicalIndicesSize =
    detail::checked_relation_position<PosType>(indicesSize, "indices-buffer size");
  return RelationType(
    Builder()
      .fromSet(fromSet)
      .toSet(toSet)
      .begins(
        typename Builder::BeginsSetBuilder().size(canonicalBeginsSize).data(begins, canonicalBeginsSize))
      .indices(typename Builder::IndicesSetBuilder()
                 .size(canonicalIndicesSize)
                 .data(indices, canonicalIndicesSize)));
}

/// \brief Reference overload for make_variable_relation (C-array-backed).
template <typename FromSet, typename ToSet, typename PosType, typename ElemType>
  requires detail::VariableRelationBufferTypes<FromSet, ToSet, PosType, ElemType>
auto make_variable_relation(FromSet& fromSet,
                            ToSet& toSet,
                            PosType* begins,
                            PosType beginsSize,
                            ElemType* indices,
                            PosType indicesSize)
{
  return make_variable_relation(&fromSet, &toSet, begins, beginsSize, indices, indicesSize);
}

/*!
 * \brief Make a static, variable-cardinality relation backed by ArrayView storage.
 *
 * \param fromSet pointer to the from-set
 * \param toSet pointer to the to-set
 * \param begins view of begin offsets
 * \param indices view of to-set positions in flat storage
 * \pre begins.size() == fromSet->size() + 1, or one for a null from-set.
 */
template <typename FromSet, typename ToSet, typename PosType, typename ElemType>
  requires detail::VariableRelationBufferTypes<FromSet, ToSet, PosType, ElemType>
auto make_variable_relation(FromSet* fromSet,
                            ToSet* toSet,
                            axom::ArrayView<PosType> begins,
                            axom::ArrayView<ElemType> indices)
{
  using BeginsIndirection = policies::ArrayViewIndirection<PosType, PosType>;
  using IndicesIndirection = policies::ArrayViewIndirection<PosType, ElemType>;
  using Cardinality = policies::VariableCardinality<PosType, BeginsIndirection>;
  using RelationType =
    StaticRelation<PosType, ElemType, Cardinality, IndicesIndirection, FromSet, ToSet>;
  using Builder = typename RelationType::RelationBuilder;

  const PosType beginsSize = detail::check_variable_relation_size<PosType>(fromSet, begins.size());
  const PosType indicesSize =
    detail::checked_relation_position<PosType>(indices.size(), "indices-buffer size");
  return RelationType(
    Builder()
      .fromSet(fromSet)
      .toSet(toSet)
      .begins(typename Builder::BeginsSetBuilder().size(beginsSize).data(begins))
      .indices(typename Builder::IndicesSetBuilder().size(indicesSize).data(indices)));
}

/// \brief Reference overload for make_variable_relation (ArrayView-backed).
template <typename FromSet, typename ToSet, typename PosType, typename ElemType>
  requires detail::VariableRelationBufferTypes<FromSet, ToSet, PosType, ElemType>
auto make_variable_relation(FromSet& fromSet,
                            ToSet& toSet,
                            axom::ArrayView<PosType> begins,
                            axom::ArrayView<ElemType> indices)
{
  return make_variable_relation(&fromSet, &toSet, begins, indices);
}

/*!
 * \brief Make a static, variable-cardinality relation backed by axom::Array storage.
 *
 * \param fromSet pointer to the from-set
 * \param toSet pointer to the to-set
 * \param begins array of begin offsets
 * \param indices array of to-set positions in flat storage
 * \pre begins.size() == fromSet->size() + 1, or one for a null from-set.
 */
template <typename FromSet, typename ToSet, typename PosType, typename ElemType>
  requires detail::VariableRelationBufferTypes<FromSet, ToSet, PosType, ElemType>
auto make_variable_relation(FromSet* fromSet,
                            ToSet* toSet,
                            axom::Array<PosType>& begins,
                            axom::Array<ElemType>& indices)
{
  using BeginsIndirection = policies::ArrayIndirection<PosType, PosType>;
  using IndicesIndirection = policies::ArrayIndirection<PosType, ElemType>;
  using Cardinality = policies::VariableCardinality<PosType, BeginsIndirection>;
  using RelationType =
    StaticRelation<PosType, ElemType, Cardinality, IndicesIndirection, FromSet, ToSet>;
  using Builder = typename RelationType::RelationBuilder;

  const PosType beginsSize = detail::check_variable_relation_size<PosType>(fromSet, begins.size());
  const PosType indicesSize =
    detail::checked_relation_position<PosType>(indices.size(), "indices-buffer size");
  return RelationType(
    Builder()
      .fromSet(fromSet)
      .toSet(toSet)
      .begins(typename Builder::BeginsSetBuilder().size(beginsSize).data(&begins))
      .indices(typename Builder::IndicesSetBuilder().size(indicesSize).data(&indices)));
}

/// \brief Reference overload for make_variable_relation (axom::Array-backed).
template <typename FromSet, typename ToSet, typename PosType, typename ElemType>
  requires detail::VariableRelationBufferTypes<FromSet, ToSet, PosType, ElemType>
auto make_variable_relation(FromSet& fromSet,
                            ToSet& toSet,
                            axom::Array<PosType>& begins,
                            axom::Array<ElemType>& indices)
{
  return make_variable_relation(&fromSet, &toSet, begins, indices);
}

/*!
 * \brief Make a static, constant-cardinality relation with a runtime stride, backed by std::vector indices.
 *
 * \param fromSet pointer to the from-set
 * \param toSet pointer to the to-set
 * \param stride  number of to-set elements per from-set element
 * \param indices to-set positions in flat storage
 * \pre indices.size() == fromSet->size() * stride, or zero for a null from-set.
 */
template <typename FromSet, typename ToSet, typename StrideType, typename ElemType>
  requires detail::RelationIndexBufferTypes<FromSet, ToSet, ElemType> &&
  detail::RelationFlatPositionConstructible<FromSet, ToSet, StrideType>
auto make_constant_relation(FromSet* fromSet,
                            ToSet* toSet,
                            StrideType stride,
                            std::vector<ElemType>& indices)
{
  using PosType = detail::RelationFlatPosition<FromSet, ToSet>;
  using IndicesIndirection = policies::STLVectorIndirection<PosType, ElemType>;
  using CTy = policies::ConstantCardinality<PosType, policies::RuntimeStride<PosType>>;
  using RelationType = StaticRelation<PosType, ElemType, CTy, IndicesIndirection, FromSet, ToSet>;
  using Builder = typename RelationType::RelationBuilder;

  const PosType canonicalStride = detail::checked_relation_stride<PosType>(stride);
  const PosType canonicalSize =
    detail::checked_relation_position<PosType>(indices.size(), "indices-buffer size");
  auto begins_builder = typename Builder::BeginsSetBuilder().stride(canonicalStride);
  detail::check_constant_relation_size(fromSet, canonicalStride, canonicalSize);
  return RelationType(
    Builder()
      .fromSet(fromSet)
      .toSet(toSet)
      .begins(begins_builder)
      .indices(typename Builder::IndicesSetBuilder().size(canonicalSize).data(&indices)));
}

/// \brief Reference overload for make_constant_relation (std::vector-backed).
template <typename FromSet, typename ToSet, typename StrideType, typename ElemType>
  requires detail::RelationIndexBufferTypes<FromSet, ToSet, ElemType> &&
  detail::RelationFlatPositionConstructible<FromSet, ToSet, StrideType>
auto make_constant_relation(FromSet& fromSet,
                            ToSet& toSet,
                            StrideType stride,
                            std::vector<ElemType>& indices)
{
  return make_constant_relation(&fromSet, &toSet, stride, indices);
}

/// \brief Make a static, constant-cardinality relation with a runtime stride, backed by C array indices.
template <typename FromSet, typename ToSet, typename StrideType, typename ElemType, typename SizeType>
  requires detail::RelationIndexBufferTypes<FromSet, ToSet, ElemType> &&
  detail::RelationFlatPositionConstructible<FromSet, ToSet, StrideType> &&
  detail::RelationFlatPositionConstructible<FromSet, ToSet, SizeType>
auto make_constant_relation(FromSet* fromSet,
                            ToSet* toSet,
                            StrideType stride,
                            ElemType* indices,
                            SizeType indicesSize)
{
  using PosType = detail::RelationFlatPosition<FromSet, ToSet>;
  using IndicesIndirection = policies::CArrayIndirection<PosType, ElemType>;
  using CTy = policies::ConstantCardinality<PosType, policies::RuntimeStride<PosType>>;
  using RelationType = StaticRelation<PosType, ElemType, CTy, IndicesIndirection, FromSet, ToSet>;
  using Builder = typename RelationType::RelationBuilder;

  const PosType canonicalStride = detail::checked_relation_stride<PosType>(stride);
  const PosType canonicalSize =
    detail::checked_relation_position<PosType>(indicesSize, "indices-buffer size");
  auto begins_builder = typename Builder::BeginsSetBuilder().stride(canonicalStride);
  detail::check_constant_relation_size(fromSet, canonicalStride, canonicalSize);
  return RelationType(
    Builder()
      .fromSet(fromSet)
      .toSet(toSet)
      .begins(begins_builder)
      .indices(
        typename Builder::IndicesSetBuilder().size(canonicalSize).data(indices, canonicalSize)));
}

/// \brief Reference overload for make_constant_relation (C-array-backed).
template <typename FromSet, typename ToSet, typename StrideType, typename ElemType, typename SizeType>
  requires detail::RelationIndexBufferTypes<FromSet, ToSet, ElemType> &&
  detail::RelationFlatPositionConstructible<FromSet, ToSet, StrideType> &&
  detail::RelationFlatPositionConstructible<FromSet, ToSet, SizeType>
auto make_constant_relation(FromSet& fromSet,
                            ToSet& toSet,
                            StrideType stride,
                            ElemType* indices,
                            SizeType indicesSize)
{
  return make_constant_relation(&fromSet, &toSet, stride, indices, indicesSize);
}

/// \brief Make a static, constant-cardinality relation with a runtime stride, backed by ArrayView indices.
template <typename FromSet, typename ToSet, typename StrideType, typename ElemType>
  requires detail::RelationIndexBufferTypes<FromSet, ToSet, ElemType> &&
  detail::RelationFlatPositionConstructible<FromSet, ToSet, StrideType>
auto make_constant_relation(FromSet* fromSet,
                            ToSet* toSet,
                            StrideType stride,
                            axom::ArrayView<ElemType> indices)
{
  using PosType = detail::RelationFlatPosition<FromSet, ToSet>;
  using IndicesIndirection = policies::ArrayViewIndirection<PosType, ElemType>;
  using CTy = policies::ConstantCardinality<PosType, policies::RuntimeStride<PosType>>;
  using RelationType = StaticRelation<PosType, ElemType, CTy, IndicesIndirection, FromSet, ToSet>;
  using Builder = typename RelationType::RelationBuilder;

  const PosType canonicalStride = detail::checked_relation_stride<PosType>(stride);
  const PosType canonicalSize =
    detail::checked_relation_position<PosType>(indices.size(), "indices-buffer size");
  auto begins_builder = typename Builder::BeginsSetBuilder().stride(canonicalStride);
  detail::check_constant_relation_size(fromSet, canonicalStride, canonicalSize);
  return RelationType(
    Builder()
      .fromSet(fromSet)
      .toSet(toSet)
      .begins(begins_builder)
      .indices(typename Builder::IndicesSetBuilder().size(canonicalSize).data(indices)));
}

/// \brief Reference overload for make_constant_relation (ArrayView-backed).
template <typename FromSet, typename ToSet, typename StrideType, typename ElemType>
  requires detail::RelationIndexBufferTypes<FromSet, ToSet, ElemType> &&
  detail::RelationFlatPositionConstructible<FromSet, ToSet, StrideType>
auto make_constant_relation(FromSet& fromSet,
                            ToSet& toSet,
                            StrideType stride,
                            axom::ArrayView<ElemType> indices)
{
  return make_constant_relation(&fromSet, &toSet, stride, indices);
}

/// \brief Make a static, constant-cardinality relation with a runtime stride, backed by axom::Array indices.
template <typename FromSet, typename ToSet, typename StrideType, typename ElemType>
  requires detail::RelationIndexBufferTypes<FromSet, ToSet, ElemType> &&
  detail::RelationFlatPositionConstructible<FromSet, ToSet, StrideType>
auto make_constant_relation(FromSet* fromSet,
                            ToSet* toSet,
                            StrideType stride,
                            axom::Array<ElemType>& indices)
{
  using PosType = detail::RelationFlatPosition<FromSet, ToSet>;
  using IndicesIndirection = policies::ArrayIndirection<PosType, ElemType>;
  using CTy = policies::ConstantCardinality<PosType, policies::RuntimeStride<PosType>>;
  using RelationType = StaticRelation<PosType, ElemType, CTy, IndicesIndirection, FromSet, ToSet>;
  using Builder = typename RelationType::RelationBuilder;

  const PosType canonicalStride = detail::checked_relation_stride<PosType>(stride);
  const PosType canonicalSize =
    detail::checked_relation_position<PosType>(indices.size(), "indices-buffer size");
  auto begins_builder = typename Builder::BeginsSetBuilder().stride(canonicalStride);
  detail::check_constant_relation_size(fromSet, canonicalStride, canonicalSize);
  return RelationType(
    Builder()
      .fromSet(fromSet)
      .toSet(toSet)
      .begins(begins_builder)
      .indices(typename Builder::IndicesSetBuilder().size(canonicalSize).data(&indices)));
}

/// \brief Reference overload for make_constant_relation (axom::Array-backed).
template <typename FromSet, typename ToSet, typename StrideType, typename ElemType>
  requires detail::RelationIndexBufferTypes<FromSet, ToSet, ElemType> &&
  detail::RelationFlatPositionConstructible<FromSet, ToSet, StrideType>
auto make_constant_relation(FromSet& fromSet,
                            ToSet& toSet,
                            StrideType stride,
                            axom::Array<ElemType>& indices)
{
  return make_constant_relation(&fromSet, &toSet, stride, indices);
}

/*!
 * \brief Make a static, constant-cardinality relation with a compile-time stride, backed by C array indices.
 *
 * \tparam STRIDE number of to-set elements per from-set element
 */
template <int STRIDE, typename FromSet, typename ToSet, typename ExplicitPosType = void, typename ElemType, typename SizeType>
  requires detail::RelationIndexBufferTypes<FromSet, ToSet, ElemType> &&
  detail::PositiveStaticStrideForPosition<
             STRIDE,
             detail::SelectedRelationFlatPosition<FromSet, ToSet, ExplicitPosType>> &&
  detail::OptionalRelationFlatPositionFor<FromSet, ToSet, ExplicitPosType> &&
  detail::SelectedRelationFlatPositionConstructible<FromSet, ToSet, ExplicitPosType, SizeType>
auto make_constant_relation_ct(FromSet* fromSet, ToSet* toSet, ElemType* indices, SizeType indicesSize)
{
  using PosType = detail::SelectedRelationFlatPosition<FromSet, ToSet, ExplicitPosType>;
  using IndicesIndirection = policies::CArrayIndirection<PosType, ElemType>;
  using StridePolicy = policies::CompileTimeStride<PosType, static_cast<PosType>(STRIDE)>;
  using CTy = policies::ConstantCardinality<PosType, StridePolicy>;
  using RelationType = StaticRelation<PosType, ElemType, CTy, IndicesIndirection, FromSet, ToSet>;
  using Builder = typename RelationType::RelationBuilder;

  const PosType canonicalSize =
    detail::checked_relation_position<PosType>(indicesSize, "indices-buffer size");
  detail::check_constant_relation_size(fromSet, static_cast<PosType>(STRIDE), canonicalSize);
  return RelationType(Builder().fromSet(fromSet).toSet(toSet).indices(
    typename Builder::IndicesSetBuilder().size(canonicalSize).data(indices, canonicalSize)));
}

/// \brief Reference overload for make_constant_relation_ct (C-array-backed).
template <int STRIDE, typename FromSet, typename ToSet, typename ExplicitPosType = void, typename ElemType, typename SizeType>
  requires detail::RelationIndexBufferTypes<FromSet, ToSet, ElemType> &&
  detail::PositiveStaticStrideForPosition<
             STRIDE,
             detail::SelectedRelationFlatPosition<FromSet, ToSet, ExplicitPosType>> &&
  detail::OptionalRelationFlatPositionFor<FromSet, ToSet, ExplicitPosType> &&
  detail::SelectedRelationFlatPositionConstructible<FromSet, ToSet, ExplicitPosType, SizeType>
auto make_constant_relation_ct(FromSet& fromSet, ToSet& toSet, ElemType* indices, SizeType indicesSize)
{
  return make_constant_relation_ct<STRIDE, FromSet, ToSet, ExplicitPosType>(&fromSet,
                                                                            &toSet,
                                                                            indices,
                                                                            indicesSize);
}

/// \brief Make a static, constant-cardinality relation with a compile-time stride, backed by std::vector indices.
template <int STRIDE, typename FromSet, typename ToSet, typename ExplicitPosType = void, typename ElemType>
  requires detail::RelationIndexBufferTypes<FromSet, ToSet, ElemType> &&
  detail::PositiveStaticStrideForPosition<
             STRIDE,
             detail::SelectedRelationFlatPosition<FromSet, ToSet, ExplicitPosType>> &&
  detail::OptionalRelationFlatPositionFor<FromSet, ToSet, ExplicitPosType>
auto make_constant_relation_ct(FromSet* fromSet, ToSet* toSet, std::vector<ElemType>& indices)
{
  using PosType = detail::SelectedRelationFlatPosition<FromSet, ToSet, ExplicitPosType>;
  using IndicesIndirection = policies::STLVectorIndirection<PosType, ElemType>;
  using StridePolicy = policies::CompileTimeStride<PosType, static_cast<PosType>(STRIDE)>;
  using CTy = policies::ConstantCardinality<PosType, StridePolicy>;
  using RelationType = StaticRelation<PosType, ElemType, CTy, IndicesIndirection, FromSet, ToSet>;
  using Builder = typename RelationType::RelationBuilder;

  const PosType canonicalSize =
    detail::checked_relation_position<PosType>(indices.size(), "indices-buffer size");
  detail::check_constant_relation_size(fromSet, static_cast<PosType>(STRIDE), canonicalSize);
  return RelationType(Builder().fromSet(fromSet).toSet(toSet).indices(
    typename Builder::IndicesSetBuilder().size(canonicalSize).data(&indices)));
}

/// \brief Reference overload for make_constant_relation_ct (std::vector-backed).
template <int STRIDE, typename FromSet, typename ToSet, typename ExplicitPosType = void, typename ElemType>
  requires detail::RelationIndexBufferTypes<FromSet, ToSet, ElemType> &&
  detail::PositiveStaticStrideForPosition<
             STRIDE,
             detail::SelectedRelationFlatPosition<FromSet, ToSet, ExplicitPosType>> &&
  detail::OptionalRelationFlatPositionFor<FromSet, ToSet, ExplicitPosType>
auto make_constant_relation_ct(FromSet& fromSet, ToSet& toSet, std::vector<ElemType>& indices)
{
  return make_constant_relation_ct<STRIDE, FromSet, ToSet, ExplicitPosType>(&fromSet, &toSet, indices);
}

/// \brief Make a static, constant-cardinality relation with a compile-time stride, backed by ArrayView indices.
template <int STRIDE, typename FromSet, typename ToSet, typename ExplicitPosType = void, typename ElemType>
  requires detail::RelationIndexBufferTypes<FromSet, ToSet, ElemType> &&
  detail::PositiveStaticStrideForPosition<
             STRIDE,
             detail::SelectedRelationFlatPosition<FromSet, ToSet, ExplicitPosType>> &&
  detail::OptionalRelationFlatPositionFor<FromSet, ToSet, ExplicitPosType>
auto make_constant_relation_ct(FromSet* fromSet, ToSet* toSet, axom::ArrayView<ElemType> indices)
{
  using PosType = detail::SelectedRelationFlatPosition<FromSet, ToSet, ExplicitPosType>;
  using IndicesIndirection = policies::ArrayViewIndirection<PosType, ElemType>;
  using StridePolicy = policies::CompileTimeStride<PosType, static_cast<PosType>(STRIDE)>;
  using CTy = policies::ConstantCardinality<PosType, StridePolicy>;
  using RelationType = StaticRelation<PosType, ElemType, CTy, IndicesIndirection, FromSet, ToSet>;
  using Builder = typename RelationType::RelationBuilder;

  const PosType canonicalSize =
    detail::checked_relation_position<PosType>(indices.size(), "indices-buffer size");
  detail::check_constant_relation_size(fromSet, static_cast<PosType>(STRIDE), canonicalSize);
  return RelationType(Builder().fromSet(fromSet).toSet(toSet).indices(
    typename Builder::IndicesSetBuilder().size(canonicalSize).data(indices)));
}

/// \brief Reference overload for make_constant_relation_ct (ArrayView-backed).
template <int STRIDE, typename FromSet, typename ToSet, typename ExplicitPosType = void, typename ElemType>
  requires detail::RelationIndexBufferTypes<FromSet, ToSet, ElemType> &&
  detail::PositiveStaticStrideForPosition<
             STRIDE,
             detail::SelectedRelationFlatPosition<FromSet, ToSet, ExplicitPosType>> &&
  detail::OptionalRelationFlatPositionFor<FromSet, ToSet, ExplicitPosType>
auto make_constant_relation_ct(FromSet& fromSet, ToSet& toSet, axom::ArrayView<ElemType> indices)
{
  return make_constant_relation_ct<STRIDE, FromSet, ToSet, ExplicitPosType>(&fromSet, &toSet, indices);
}

/// \brief Make a static, constant-cardinality relation with a compile-time stride, backed by axom::Array indices.
template <int STRIDE, typename FromSet, typename ToSet, typename ExplicitPosType = void, typename ElemType>
  requires detail::RelationIndexBufferTypes<FromSet, ToSet, ElemType> &&
  detail::PositiveStaticStrideForPosition<
             STRIDE,
             detail::SelectedRelationFlatPosition<FromSet, ToSet, ExplicitPosType>> &&
  detail::OptionalRelationFlatPositionFor<FromSet, ToSet, ExplicitPosType>
auto make_constant_relation_ct(FromSet* fromSet, ToSet* toSet, axom::Array<ElemType>& indices)
{
  using PosType = detail::SelectedRelationFlatPosition<FromSet, ToSet, ExplicitPosType>;
  using IndicesIndirection = policies::ArrayIndirection<PosType, ElemType>;
  using StridePolicy = policies::CompileTimeStride<PosType, static_cast<PosType>(STRIDE)>;
  using CTy = policies::ConstantCardinality<PosType, StridePolicy>;
  using RelationType = StaticRelation<PosType, ElemType, CTy, IndicesIndirection, FromSet, ToSet>;
  using Builder = typename RelationType::RelationBuilder;

  const PosType canonicalSize =
    detail::checked_relation_position<PosType>(indices.size(), "indices-buffer size");
  detail::check_constant_relation_size(fromSet, static_cast<PosType>(STRIDE), canonicalSize);
  return RelationType(Builder().fromSet(fromSet).toSet(toSet).indices(
    typename Builder::IndicesSetBuilder().size(canonicalSize).data(&indices)));
}

/// \brief Reference overload for make_constant_relation_ct (axom::Array-backed).
template <int STRIDE, typename FromSet, typename ToSet, typename ExplicitPosType = void, typename ElemType>
  requires detail::RelationIndexBufferTypes<FromSet, ToSet, ElemType> &&
  detail::PositiveStaticStrideForPosition<
             STRIDE,
             detail::SelectedRelationFlatPosition<FromSet, ToSet, ExplicitPosType>> &&
  detail::OptionalRelationFlatPositionFor<FromSet, ToSet, ExplicitPosType>
auto make_constant_relation_ct(FromSet& fromSet, ToSet& toSet, axom::Array<ElemType>& indices)
{
  return make_constant_relation_ct<STRIDE, FromSet, ToSet, ExplicitPosType>(&fromSet, &toSet, indices);
}

/// \}

}  // end namespace axom::slam
