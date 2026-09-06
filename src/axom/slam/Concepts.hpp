// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file Concepts.hpp
 * \brief Semantic concepts for SLAM sets, relations, and maps.
 *
 * Public object concepts ignore top-level cv/ref qualification,
 * which are normalized internally in the detail namespace.
 * Standard range integration is provided separately by Ranges.hpp.
 */

#pragma once

#include <concepts>
#include <limits>
#include <type_traits>
#include <utility>

namespace axom::slam
{
namespace detail
{
template <typename T>
using model_t = std::remove_cvref_t<T>;

template <typename T>
using position_t = typename model_t<T>::PositionType;

template <typename T>
using element_t = typename model_t<T>::ElementType;

/// Integer arguments may be unsigned, but cannot be Boolean.
template <typename T>
concept PositionValueLike = std::integral<model_t<T>> && !std::same_as<model_t<T>, bool>;

/// The operations needed to traverse values, without prescribing iterator aliases.
template <typename Iterator, typename Sentinel, typename Value>
concept IteratesAs = requires(Iterator it, Sentinel end) {
  { *it } -> std::convertible_to<Value>;
  ++it;
  { it != end } -> std::convertible_to<bool>;
};

template <typename R, typename Value>
concept IterableValues = requires(const model_t<R>& range) {
  requires IteratesAs<decltype(range.begin()), decltype(range.end()), Value>;
};

template <typename R, typename Value>
concept SizedValues = IterableValues<R, Value> && requires(const model_t<R>& range) {
  { range.size() } -> PositionValueLike;
};

/// Scalar references may add constness, but cannot discard element qualification.
template <typename Reference, typename Data>
concept MapReferenceFor = std::is_lvalue_reference_v<Reference> &&
  (std::same_as<std::remove_reference_t<Reference>, Data> ||
   std::same_as<std::remove_reference_t<Reference>, const Data>);
}  // namespace detail

/// \brief A signed integral position type supported by the current SLAM arithmetic.
/// Tagged positions need a separate point/difference/size design, not an opt-in here.
template <typename T>
concept PositionLike = std::signed_integral<detail::model_t<T>>;

/// \brief A sized collection with positional access to its elements.
/// size() is nonnegative; empty() agrees with size() == 0. Access requires a valid position.
template <typename T>
concept SetLike = PositionLike<detail::position_t<T>> &&
  requires(const detail::model_t<T>& set, detail::position_t<T> pos) {
    typename detail::element_t<T>;
    { set.size() } -> std::same_as<detail::position_t<T>>;
    { set.empty() } -> std::convertible_to<bool>;
    { set.at(pos) } -> std::convertible_to<detail::element_t<T>>;
  };

/// \brief A SetLike type whose const traversal visits the elements in positional order.
template <typename T>
concept IterableSetLike = SetLike<T> && detail::IterableValues<T, detail::element_t<T>>;

/// \brief A container that can check its own internal consistency.
/// Validation is a separate capability, not a requirement of every set or relation.
template <typename T>
concept Validatable = requires(const detail::model_t<T>& container) {
  { container.isValid(false) } -> std::convertible_to<bool>;
};

/// \brief A set of coordinate pairs with a row of second-set positions per first-set position.
/// Concatenating rows in first-set order agrees with at(flat). The flat size
/// equals the sum of the row sizes. Search and flat projections are conveniences.
template <typename T>
concept BivariateSetLike = SetLike<T> && SetLike<typename detail::model_t<T>::FirstSetType> &&
  SetLike<typename detail::model_t<T>::SecondSetType> &&
  requires(detail::element_t<T> coordinate) {
    requires std::same_as<std::remove_cvref_t<decltype(coordinate.first)>,
                          typename detail::model_t<T>::FirstSetType::PositionType>;
    requires std::same_as<std::remove_cvref_t<decltype(coordinate.second)>,
                          typename detail::model_t<T>::SecondSetType::PositionType>;
  } &&
  requires(const detail::model_t<T>& set,
           typename detail::model_t<T>::FirstSetType::PositionType first) {
    { set.getFirstSet() } -> std::same_as<const typename detail::model_t<T>::FirstSetType*>;
    { set.getSecondSet() } -> std::same_as<const typename detail::model_t<T>::SecondSetType*>;
    requires detail::SizedValues<decltype(set.getElements(first)),
                                 typename detail::model_t<T>::SecondSetType::PositionType>;
  };

/// \brief A relation with a row of to-set positions for each from-set position.
/// Row entries must identify valid to-set positions. The sets may themselves
/// contain coordinates. No flat storage, validation method, or row alias is required.
template <typename T>
concept RelationLike = SetLike<typename detail::model_t<T>::FromSetType> &&
  SetLike<typename detail::model_t<T>::ToSetType> &&
  requires(const detail::model_t<T>& relation,
           typename detail::model_t<T>::FromSetType::PositionType from) {
    { relation.fromSet() } -> std::same_as<const typename detail::model_t<T>::FromSetType*>;
    { relation.toSet() } -> std::same_as<const typename detail::model_t<T>::ToSetType*>;
    requires detail::SizedValues<decltype(relation[from]),
                                 typename detail::model_t<T>::ToSetType::PositionType>;
  };

/// \brief Values associated with entries, addressed by entry position and local component.
/// index(pos) identifies the associated set element. flatValue(pos, component)
/// accesses one scalar component, including for tensor-valued maps. "flat" here
/// selects an entry of a bivariate set, not a global component-storage position.
/// Bound maps have positive numComp(); a default unbound SubMap may be empty with zero components.
template <typename T>
concept MapLike = PositionLike<detail::position_t<T>> &&
  requires(detail::model_t<T>& map,
           const detail::model_t<T>& constMap,
           detail::position_t<T> pos,
           detail::position_t<T> component) {
    typename detail::model_t<T>::DataType;
    { constMap.size() } -> std::same_as<detail::position_t<T>>;
    { constMap.numComp() } -> std::same_as<detail::position_t<T>>;
    constMap.index(pos);
    requires std::is_object_v<std::remove_cvref_t<decltype(constMap.index(pos))>>;
    {
      map.flatValue(pos, component)
    } -> detail::MapReferenceFor<typename detail::model_t<T>::DataType>;
    {
      constMap.flatValue(pos, component)
    } -> detail::MapReferenceFor<typename detail::model_t<T>::DataType>;
  };

/// \brief A MapLike type explicitly bound to all positions of exactly S.
/// MappedSetType names that binding and set() returns it. size() agrees with
/// set()->size(), and index(pos) agrees with set()->at(pos). SubMap has no such
/// binding: its set() selects parent positions, and index() projects their elements.
template <typename M, typename S>
concept MapOver = MapLike<M> && SetLike<S> &&
  std::same_as<typename detail::model_t<M>::MappedSetType, detail::model_t<S>> &&
  std::same_as<detail::position_t<M>, detail::position_t<S>> &&
  requires(const detail::model_t<M>& map, detail::position_t<M> pos) {
    { map.set() } -> std::same_as<const detail::model_t<S>*>;
    { map.index(pos) } -> std::convertible_to<detail::element_t<S>>;
  };

/// \brief A non-reference type with a trivially copyable C++ representation.
/// This does not certify device-callable operations or the accessibility/lifetime
/// of referenced objects. Device capture is a separate view-conversion contract.
template <typename T>
concept TriviallyCopyableRepresentation =
  !std::is_reference_v<T> && std::is_trivially_copyable_v<std::remove_cv_t<T>>;

namespace detail
{
/// Type-level representability of nonnegative positions and sizes.
template <typename Position, typename RepresentedPosition>
concept PositionCanRepresent = PositionLike<Position> && PositionLike<RepresentedPosition> &&
  (std::numeric_limits<model_t<Position>>::digits >=
   std::numeric_limits<model_t<RepresentedPosition>>::digits);

template <int Stride, typename Position>
concept PositiveStaticStrideForPosition = PositionLike<Position> && (Stride > 0) &&
  (Stride <= std::numeric_limits<model_t<Position>>::max());

template <int Stride, typename Set>
concept PositiveStaticStrideFor =
  SetLike<Set> && PositiveStaticStrideForPosition<Stride, position_t<Set>>;

template <typename Set, typename Value>
concept SetPositionConvertible =
  SetLike<Set> && PositionValueLike<Value> && std::convertible_to<model_t<Value>, position_t<Set>>;

template <typename Set, typename Position>
concept SetPositionSame = SetLike<Set> && std::same_as<model_t<Position>, position_t<Set>>;

template <typename Set, typename Position>
concept OptionalSetPositionSame =
  SetLike<Set> && (std::same_as<model_t<Position>, void> || SetPositionSame<Set, Position>);

/// Selected positions used by SubMap and the current BivariateMap adapter.
template <typename R, typename Position>
concept FlatRangeOver = IterableSetLike<R> && std::same_as<element_t<R>, model_t<Position>> &&
  requires(const model_t<R>& range, position_t<R> pos) {
    { range[pos] } -> std::convertible_to<model_t<Position>>;
  };

/// Operations consumed by the current BivariateMap implementation, not the set abstraction.
template <typename T>
concept BivariateMapSet = BivariateSetLike<T> && IterableSetLike<T> && Validatable<T> &&
  requires(const model_t<T>& set,
           position_t<T> flat,
           typename model_t<T>::FirstSetType::PositionType first,
           typename model_t<T>::SecondSetType::PositionType second) {
    typename model_t<T>::SubsetType;
    { model_t<T>::INVALID_POS } -> std::convertible_to<position_t<T>>;
    { set.elementRangeSet(first) } -> FlatRangeOver<position_t<T>>;
    {
      set.findElementIndex(first, second)
    } -> std::same_as<decltype(set.getElements(first).size())>;
    { set.findElementFlatIndex(first, second) } -> std::same_as<position_t<T>>;
    { set.flatToFirstIndex(flat) } -> std::same_as<typename model_t<T>::FirstSetType::PositionType>;
    {
      set.flatToSecondIndex(flat)
    } -> std::same_as<typename model_t<T>::SecondSetType::PositionType>;
    { set.firstSetSize() } -> std::same_as<typename model_t<T>::FirstSetType::PositionType>;
    { set.secondSetSize() } -> std::same_as<typename model_t<T>::SecondSetType::PositionType>;
    { set.size(first) } -> std::same_as<position_t<T>>;
  };

/// Flat storage consumed by RelationSet, not a public relation requirement.
template <typename T>
concept RelationSetSource = RelationLike<model_t<T>> && Validatable<model_t<T>> &&
  requires {
    typename model_t<T>::FlatPositionType;
    typename model_t<T>::RelationSubset;
  } && PositionLike<typename model_t<T>::FlatPositionType> &&
  requires(const model_t<T>& relation,
           typename model_t<T>::FromSetType::PositionType fromPosition,
           typename model_t<T>::FlatPositionType flatPosition) {
    { relation[fromPosition] } -> std::convertible_to<typename model_t<T>::RelationSubset>;
    { relation.size(fromPosition) } -> PositionValueLike;
    { relation.offset(fromPosition) } -> std::convertible_to<typename model_t<T>::FlatPositionType>;
    {
      static_cast<typename model_t<T>::FromSetType::PositionType>(relation.firstIndex(flatPosition))
    } -> std::same_as<typename model_t<T>::FromSetType::PositionType>;
    {
      relation.relationData().size()
    } -> std::convertible_to<typename model_t<T>::FlatPositionType>;
    {
      relation.relationData()[flatPosition]
    } -> std::convertible_to<typename model_t<T>::ToSetType::PositionType>;
  };

template <typename T>
concept HasIndirectionAssociatedTypes = requires {
  typename T::IndirectionResult;
  typename T::ConstIndirectionResult;
  typename T::IndirectionBufferType;
  typename T::IndirectionPtrType;
};

template <typename T>
concept HasTypedIndirectionAssociatedTypes = requires {
  typename T::PositionType;
  typename T::ElementType;
};

template <typename T>
using policy_default_t = std::remove_cv_t<decltype(T::DEFAULT_VALUE)>;

/// \brief A scalar runtime or compile-time value policy.
template <typename T>
concept ValuePolicy = requires(const T& policy) {
  typename T::TagType;
  typename T::IntType;
  { policy.value() } -> std::same_as<typename T::IntType>;
  { policy.isValid(false) } -> std::convertible_to<bool>;
};

/// \brief A policy that reports a size and whether that size is empty.
template <typename T>
concept SizePolicy = requires(const T& policy) {
  T::DEFAULT_VALUE;
  { policy.size() } -> std::same_as<policy_default_t<T>>;
  { policy.empty() } -> std::convertible_to<bool>;
  { policy.isValid(false) } -> std::convertible_to<bool>;
};

/// \brief A scalar value policy that reports an offset.
template <typename T>
concept OffsetPolicy = ValuePolicy<T> && requires(const T& policy) {
  T::DEFAULT_VALUE;
  { policy.offset() } -> std::same_as<typename T::IntType>;
};

/*!
 * \brief The common capability shared by scalar and multi-dimensional stride policies.
 *
 * Use OrderedSetStridePolicyFor or MapStridePolicyFor when checking whether a stride
 * can actually be substituted into one of those owners.
 */
template <typename T>
concept StridePolicy = requires(const T& policy) {
  typename T::IndexType;
  typename T::ShapeType;
  T::NumDims;
  { T::DefaultSize() } -> std::same_as<typename T::ShapeType>;
  { policy.stride() } -> std::same_as<typename T::IndexType>;
  { policy.shape() } -> std::same_as<typename T::ShapeType>;
};

//---- substitutability --------------------------------------------------
//
// The four '*PolicyFor' concepts below check whether the policy can be substituted
// into an owner indexed by Position. These name the individual clauses,
// so a failed constraint identifies the problem.

/// \brief The policy's default value is exactly \a Position, and it is constructible from one.
template <typename T, typename Position>
concept PolicyDefaultedOver = PositionLike<Position> && requires { T::DEFAULT_VALUE; } &&
  std::same_as<policy_default_t<T>, Position> && std::constructible_from<T, Position>;

/// \brief A ValuePolicy whose scalar value type is exactly \a Position.
template <typename T, typename Position>
concept ScalarValuePolicyOver =
  ValuePolicy<T> && PolicyDefaultedOver<T, Position> && std::same_as<typename T::IntType, Position>;

/// \brief A StridePolicy carrying a single scalar stride measured in \a Position.
template <typename T, typename Position>
concept ScalarStridePolicyOver = StridePolicy<T> && PositionLike<Position> && (T::NumDims == 1) &&
  std::same_as<typename T::IndexType, Position> && std::same_as<typename T::ShapeType, Position>;

/// \brief A stride policy that reports one stride per dimension.
template <typename T>
concept ExposesPerDimensionStrides = requires(const T& policy) {
  { policy.strides() } -> std::same_as<typename T::ShapeType>;
};

//---- substitutability into a specific owner ----------------------------------

/// \brief A SizePolicy usable by a set whose position type is \a Position.
template <typename T, typename Position>
concept SetSizePolicyFor = SizePolicy<T> && PolicyDefaultedOver<T, Position>;

/// \brief An OffsetPolicy usable by an OrderedSet whose position type is \a Position.
template <typename T, typename Position>
concept OrderedSetOffsetPolicyFor = OffsetPolicy<T> && ScalarValuePolicyOver<T, Position>;

/*!
 * \brief A scalar stride policy usable by OrderedSet with \a Position.
 *
 * OrderedSet constructs its stride from a position and validates it as a scalar value policy.
 * Multi-dimensional map strides do not satisfy this refinement.
 */
template <typename T, typename Position>
concept OrderedSetStridePolicyFor =
  ScalarValuePolicyOver<T, Position> && ScalarStridePolicyOver<T, Position>;

/*!
 * \brief A scalar or multi-dimensional stride policy usable by Map with \a Position.
 *
 * A map is more permissive than an ordered set: its stride index type only has
 * to convert to the map's position type, and it may have more than one dimension.
 */
template <typename T, typename Position>
concept MapStridePolicyFor = StridePolicy<T> && PositionLike<Position> &&
  PositionLike<typename T::IndexType> && std::convertible_to<typename T::IndexType, Position> &&
  (T::NumDims > 0) && std::constructible_from<T, typename T::ShapeType> &&
  ((T::NumDims == 1) || ExposesPerDimensionStrides<T>);

//------------------------------------------------------------------------------
// Indirection policies
//------------------------------------------------------------------------------

//---- base capability ---------------------------------------------------------

/*!
 * \brief The common storage/indirection-policy capability.
 *
 * Use IndirectionPolicyFor when the calling position type is available and the
 * indirection operation itself should also be checked.
 * Use OrderedSetIndirectionPolicyFor or MapIndirectionPolicyFor
 * when checking substitutability into those owners.
 */
template <typename T>
concept IndirectionPolicy = HasIndirectionAssociatedTypes<T> && requires(const T& policy) {
  T::DeviceAccessible;
  { policy.hasIndirection() } -> std::convertible_to<bool>;
};

/// \brief An IndirectionPolicy callable with Position.
template <typename T, typename Position>
concept IndirectionPolicyFor =
  IndirectionPolicy<T> && requires(T& policy, const T& constPolicy, Position pos) {
    { policy.indirection(pos) } -> std::convertible_to<typename T::IndirectionResult>;
    { constPolicy.indirection(pos) } -> std::convertible_to<typename T::ConstIndirectionResult>;
    requires(!std::is_reference_v<typename T::IndirectionResult> ||
             std::same_as<decltype(policy.indirection(pos)), typename T::IndirectionResult>);
    requires(!std::is_reference_v<typename T::ConstIndirectionResult> ||
             std::same_as<decltype(constPolicy.indirection(pos)), typename T::ConstIndirectionResult>);
  };

//---- substitutability atoms --------------------------------------------------

/*!
 * \brief The policy indirects to exactly \a Data.
 *
 * Compares cv but not ref: a policy over `const double` 
 * is a different policy from one over `double`.
 */
template <typename T, typename Data>
concept IndirectsExactly = std::same_as<typename T::ElementType, std::remove_reference_t<Data>> &&
  (std::same_as<std::remove_reference_t<typename T::IndirectionResult>, std::remove_reference_t<Data>> ||
   std::same_as<std::remove_reference_t<typename T::IndirectionResult>,
                const std::remove_reference_t<Data>>) &&
  (std::same_as<std::remove_reference_t<typename T::ConstIndirectionResult>, std::remove_reference_t<Data>> ||
   std::same_as<std::remove_reference_t<typename T::ConstIndirectionResult>,
                const std::remove_reference_t<Data>>);

/// \brief Both access paths return stable lvalue references, as Map's element access requires.
template <typename T>
concept YieldsStableReferences = std::is_lvalue_reference_v<typename T::IndirectionResult> &&
  std::is_lvalue_reference_v<typename T::ConstIndirectionResult>;

/*!
 * \brief The policy names the pointer types its static accessors return.
 *
 * Both are fixed by the corresponding result type, so this states a consistency
 * requirement rather than a free choice. Map exposes them through `data_ptr()`.
 */
template <typename T>
concept HasResultPointerAliases =
  requires {
    typename T::ResultPtr;
    typename T::ConstResultPtr;
  } &&
  std::same_as<typename T::ResultPtr,
               std::add_pointer_t<std::remove_reference_t<typename T::IndirectionResult>>> &&
  std::same_as<typename T::ConstResultPtr,
               std::add_pointer_t<std::remove_reference_t<typename T::ConstIndirectionResult>>>;

/*!
 * \brief Static, pointer-returning access to the policy's buffer.
 *
 * \note Both the positioned and the whole-buffer form are required:
 *    Map's element access calls the former and `Map::data_ptr()` calls the latter.
 */
template <typename T, typename Position>
concept HasStaticBufferAccess = requires(typename T::IndirectionBufferType& buffer,
                                         const typename T::IndirectionBufferType& constBuffer,
                                         Position pos) {
  { T::getIndirection(buffer, pos) } -> std::same_as<typename T::ResultPtr>;
  { T::getConstIndirection(constBuffer, pos) } -> std::same_as<typename T::ConstResultPtr>;
  { T::getIndirection(buffer) } -> std::same_as<typename T::ResultPtr>;
  { T::getConstIndirection(constBuffer) } -> std::same_as<typename T::ConstResultPtr>;
};

/// \brief The buffer reports its extent, and can be resized when the policy owns it.
template <typename T, typename Position>
concept HasSizedBuffer = requires {
  T::IsMutableBuffer;
  std::integral_constant<bool, T::IsMutableBuffer> {};
} && requires(const typename T::IndirectionBufferType& buffer) {
  { buffer.size() } -> std::convertible_to<Position>;
  { buffer.empty() } -> std::convertible_to<bool>;
} && (!T::IsMutableBuffer || requires(typename T::IndirectionBufferType& buffer, Position size) {
                           buffer.resize(size);
                         });

/// \brief Default-constructible, and bindable to an existing buffer, as OrderedSet requires.
template <typename T>
concept BindableIndirection =
  std::default_initializable<T> && std::constructible_from<T, typename T::IndirectionPtrType>;

/// \brief Validates a (size, offset, stride) triple against the buffer it indirects through.
template <typename T, typename Position>
concept ValidatesSetRange =
  requires(const T& policy, Position size, Position offset, Position stride) {
    { policy.isValid(size, offset, stride, false) } -> std::convertible_to<bool>;
  };

//---- substitutability into a specific owner ----------------------------------

/*!
 * \brief An indirection policy usable by OrderedSet over Position and Element.
 *
 * \note \a Element keeps its cv-qualification: a policy over `const double` is a
 *       different policy from one over `double`.
 */
template <typename T, typename Position, typename Element>
concept OrderedSetIndirectionPolicyFor =
  IndirectionPolicyFor<T, Position> && PositionLike<Position> &&
  HasTypedIndirectionAssociatedTypes<T> && std::same_as<typename T::PositionType, Position> &&
  IndirectsExactly<T, Element> && BindableIndirection<T> && ValidatesSetRange<T, Position>;

/*!
 * \brief An indirection policy providing Map's buffer and static access API.
 *
 * Both access paths must return stable lvalue references, and their pointer aliases must point
 * to the same cv-qualified value types. Const access may retain shallow view semantics.
 *
 * \note \a Data keeps its cv-qualification, as for OrderedSetIndirectionPolicyFor.
 */
template <typename T, typename Position, typename Data>
concept MapIndirectionPolicyFor = IndirectionPolicy<T> && PositionLike<Position> &&
  HasTypedIndirectionAssociatedTypes<T> && std::same_as<typename T::PositionType, Position> &&
  IndirectsExactly<T, Data> && YieldsStableReferences<T> && HasResultPointerAliases<T> &&
  HasStaticBufferAccess<T, Position> && HasSizedBuffer<T, Position>;

/// \brief A MapIndirectionPolicyFor that can allocate and initialize its buffer.
template <typename T, typename Position, typename Data>
concept AllocatingMapIndirectionPolicyFor = MapIndirectionPolicyFor<T, Position, Data> &&
  requires(Position size, const std::remove_cvref_t<Data>& value, int allocatorId) {
    { T::create(size, value, allocatorId) } -> std::same_as<typename T::IndirectionBufferType>;
  };

}  // namespace detail

// Public policy extension protocols.
// The owner compositions and diagnostic clauses above remain implementation details.

/// Reports size, emptiness, and validity, with a default size value.
template <typename T>
concept SizePolicy = detail::SizePolicy<detail::model_t<T>>;

/// Reports scalar component count and shape for scalar or multidimensional strides.
template <typename T>
concept StridePolicy = detail::StridePolicy<detail::model_t<T>>;

/// Reports a scalar offset, its default value, and validity.
template <typename T>
concept OffsetPolicy = detail::OffsetPolicy<detail::model_t<T>>;

/// Callable, bindable indirection for an OrderedSet.
template <typename T, typename Position, typename Element>
concept OrderedSetIndirectionPolicyFor =
  detail::OrderedSetIndirectionPolicyFor<detail::model_t<T>, detail::model_t<Position>, Element>;

/// Sized-buffer and stable-reference access for a Map.
template <typename T, typename Position, typename Data>
concept MapIndirectionPolicyFor =
  detail::MapIndirectionPolicyFor<detail::model_t<T>, detail::model_t<Position>, Data>;

/// Map storage that can also allocate and initialize its buffer.
template <typename T, typename Position, typename Data>
concept AllocatingMapIndirectionPolicyFor =
  detail::AllocatingMapIndirectionPolicyFor<detail::model_t<T>, detail::model_t<Position>, Data>;

}  // namespace axom::slam
