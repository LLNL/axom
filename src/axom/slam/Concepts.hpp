// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file Concepts.hpp
 *
 * \brief C++20 concepts for SLAM containers, policies, and index properties.
 *
 * The concepts in this header describe public operations and relationships between associated types.
 * They do not require concrete SLAM implementations or include the standard ranges library.
 *
 * \note We normalize most types (by removing cvref) at the public boundary, 
 *   so `SetLike<const MySet&>` and `SetLike<MySet>` agree.
 *   We do not normalize a container's \a Element and a map's \a Data 
 *   since their cv-qualification is part of the contract, 
 *   e.g. an indirection policy over `const double` is different than one over `double`.
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
/// \brief The bare value type a concept is checked against.
template <typename T>
using model_t = std::remove_cvref_t<T>;
}  // namespace detail

/*!
 * \brief Opt-in customization for a non-integral SLAM position type.
 *
 * A specialized type is still responsible for providing the arithmetic
 * and ordering required by the SLAM APIs in which it is used.
 */
template <typename T>
inline constexpr bool enable_position_like = false;

namespace detail::model
{
//------------------------------------------------------------------------------
// Associated-type detectors. Each names the member typedefs that distinguish
// one SLAM concept from another; the semantic requirements are layered on top.
//------------------------------------------------------------------------------

template <typename T>
concept HasSetAssociatedTypes = requires {
  typename T::PositionType;
  typename T::ElementType;
};

template <typename T>
concept HasBivariateSetAssociatedTypes = requires {
  typename T::FirstSetType;
  typename T::SecondSetType;
};

template <typename T>
concept HasRelationAssociatedTypes = requires {
  typename T::FromSetType;
  typename T::ToSetType;
};

template <typename T>
concept HasMapAssociatedTypes = requires {
  typename T::DataType;
  typename T::SetPosition;
  typename T::SetElement;
  typename T::ValueType;
  typename T::ConstValueType;
};

template <typename T>
concept HasUnivariateMapAssociatedTypes = requires { typename T::SetType; };

template <typename T>
concept HasBivariateMapAssociatedTypes = requires { typename T::BivariateSetType; };

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

//------------------------------------------------------------------------------
// Index properties
//------------------------------------------------------------------------------

/// \brief A built-in signed integral or explicitly opted-in SLAM position type.
template <typename T>
concept PositionLike = std::signed_integral<T> || (!std::integral<T> && enable_position_like<T>);

/*!
 * \brief A non-Boolean integral or an opted-in value convertible to a position.
 *
 * \note Unlike the other concepts here, this one normalizes its own argument.
 *       It is applied to deduced expression types in return-type-constraint position,
 *       where the value category of the expression is incidental.
 */
template <typename T>
concept PositionValueLike =
  (std::integral<std::remove_cvref_t<T>> && !std::same_as<std::remove_cvref_t<T>, bool>) ||
  enable_position_like<std::remove_cvref_t<T>>;

template <typename Position, typename RepresentedPosition>
consteval bool positionTypeCanRepresent()
{
  if constexpr(std::integral<Position> && std::integral<RepresentedPosition>)
  {
    // Positions and sizes are nonnegative, so compare the number of value bits.
    return std::numeric_limits<Position>::digits >= std::numeric_limits<RepresentedPosition>::digits;
  }
  else
  {
    return std::constructible_from<Position, RepresentedPosition>;
  }
}

template <typename Position, typename RepresentedPosition>
concept PositionCanRepresent = PositionLike<Position> && PositionLike<RepresentedPosition> &&
  positionTypeCanRepresent<Position, RepresentedPosition>();

template <int Stride, typename Position>
consteval bool positiveStaticStrideRepresentable()
{
  if constexpr(std::integral<Position>)
  {
    return std::in_range<Position>(Stride);
  }
  else
  {
    // Opted-in position types state their own construction contract.
    return std::constructible_from<Position, int>;
  }
}

template <int Stride, typename Position>
concept PositiveStaticStrideForPosition =
  PositionLike<Position> && (Stride > 0) && positiveStaticStrideRepresentable<Stride, Position>();

//------------------------------------------------------------------------------
// Sets
//------------------------------------------------------------------------------

/*!
 * \brief A set: a sized, possibly empty sequence of elements addressed by position.
 *
 * This is the base of the set hierarchy. BivariateSetLike refines it, so generic code
 * that only needs size and positional access can be written once against SetLike.
 * Use UnivariateSetLike where a bivariate set must be excluded.
 */
template <typename T>
concept SetLike = HasSetAssociatedTypes<T> && PositionLike<typename T::PositionType> &&
  requires(const T& set, typename T::PositionType pos) {
    { set.size() } -> std::same_as<typename T::PositionType>;
    { set.empty() } -> std::convertible_to<bool>;
    { set.at(pos) } -> std::convertible_to<typename T::ElementType>;
  };

/// \brief A container that can check its own internal consistency.
template <typename T>
concept Validatable = requires(const T& container) {
  { container.isValid(false) } -> std::convertible_to<bool>;
};

/*!
 * \brief A set that is not a bivariate set.
 *
 * The input sets to a relation, of a product set and the domain of a univariate Map
 * are all univariate.
 */
template <typename T>
concept UnivariateSetLike = SetLike<T> && !HasBivariateSetAssociatedTypes<T>;

/*!
 * \brief A SetLike type with const iteration over its elements.
 *
 * This concept states SLAM's ordered-set surface without importing the standard ranges taxonomy.
 * Standard iterator and range categories are tested separately.
 */
template <typename T>
concept OrderedSetLike = SetLike<T> && requires(const T& set) {
  set.begin();
  { set.end() } -> std::same_as<decltype(set.begin())>;
  { *set.begin() } -> std::convertible_to<typename T::ElementType>;
};

template <int Stride, typename Set>
concept PositiveStaticStrideFor =
  UnivariateSetLike<Set> && PositiveStaticStrideForPosition<Stride, typename Set::PositionType>;

/*!
 * \brief A set whose elements are indexed by positions from two component sets.
 *
 * PositionType indexes the flattened sequence. ElementType is a coordinate
 * whose `first` and `second` members have the exact endpoint position types.
 * Row-local positions are inferred from row operations rather than imposed as
 * another required associated type.
 */
template <typename T>
concept BivariateSetLike = SetLike<T> && HasBivariateSetAssociatedTypes<T> &&
  UnivariateSetLike<typename T::FirstSetType> && UnivariateSetLike<typename T::SecondSetType> &&
  requires(typename T::ElementType coordinate) {
    requires std::same_as<std::remove_cvref_t<decltype(coordinate.first)>,
                          typename T::FirstSetType::PositionType>;
    requires std::same_as<std::remove_cvref_t<decltype(coordinate.second)>,
                          typename T::SecondSetType::PositionType>;
  } &&
  requires(const T& set,
           typename T::PositionType flatPosition,
           typename T::FirstSetType::PositionType firstPosition,
           typename T::SecondSetType::PositionType secondPosition) {
    { set.getFirstSet() } -> std::same_as<const typename T::FirstSetType*>;
    { set.getSecondSet() } -> std::same_as<const typename T::SecondSetType*>;
    { set.getElements(firstPosition).size() } -> PositionValueLike;
    set.getElements(firstPosition).begin();
    {
      set.getElements(firstPosition).end()
    } -> std::same_as<decltype(set.getElements(firstPosition).begin())>;
    {
      *set.getElements(firstPosition).begin()
    } -> std::convertible_to<typename T::SecondSetType::PositionType>;
    {
      set.findElementIndex(firstPosition, secondPosition)
    } -> std::same_as<decltype(set.getElements(firstPosition).size())>;
    {
      set.findElementFlatIndex(firstPosition, secondPosition)
    } -> std::same_as<typename T::PositionType>;
    { set.flatToFirstIndex(flatPosition) } -> std::same_as<typename T::FirstSetType::PositionType>;
    {
      set.flatToSecondIndex(flatPosition)
    } -> std::same_as<typename T::SecondSetType::PositionType>;
    { set.firstSetSize() } -> std::same_as<typename T::FirstSetType::PositionType>;
    { set.secondSetSize() } -> std::same_as<typename T::SecondSetType::PositionType>;
    { set.size(firstPosition) } -> std::same_as<typename T::PositionType>;
  };

/// \brief An ordered set of flat positions (the shape returned by elementRangeSet())
template <typename R, typename Position>
concept FlatRangeOver = OrderedSetLike<std::remove_cvref_t<R>> &&
  std::same_as<typename std::remove_cvref_t<R>::ElementType, Position> &&
  requires(const std::remove_cvref_t<R>& range,
           typename std::remove_cvref_t<R>::PositionType pos) {
    { range[pos] } -> std::convertible_to<Position>;
  };

/*!
 * \brief A BivariateSetLike type that a BivariateMap can bind field data over.
 *
 * BivariateSetLike models only the coordinate structure.
 * This adds what would be needed to bind a map.
 */
template <typename T>
concept BivariateMapDomain = BivariateSetLike<T> && OrderedSetLike<T> && Validatable<T> &&
  requires(const T& set, typename T::FirstSetType::PositionType firstPosition) {
    { T::INVALID_POS } -> std::convertible_to<typename T::PositionType>;
    { set.elementRangeSet(firstPosition) } -> FlatRangeOver<typename T::PositionType>;
  };

//------------------------------------------------------------------------------
// Relations
//------------------------------------------------------------------------------

/*!
 * \brief A relation that exposes its two sets and a const iterable row for a
 * position in the from-set.
 *
 * A relation row is selected by the FromSetType's PositionType and contains
 * positions in the ToSetType. Endpoint elements are obtained by projecting
 * those positions through the endpoint sets. A flattened representation is an
 * implementation capability, not a requirement of the relation abstraction.
 */
template <typename T>
concept RelationLike = HasRelationAssociatedTypes<T> &&
  UnivariateSetLike<typename T::FromSetType> && UnivariateSetLike<typename T::ToSetType> &&
  requires(const T& relation, typename T::FromSetType::PositionType fromPosition) {
    { relation.fromSet() } -> std::same_as<const typename T::FromSetType*>;
    { relation.toSet() } -> std::same_as<const typename T::ToSetType*>;
    { relation[fromPosition].size() } -> PositionValueLike;
    relation[fromPosition].begin();
    { relation[fromPosition].end() } -> std::same_as<decltype(relation[fromPosition].begin())>;
    { *relation[fromPosition].begin() } -> std::convertible_to<typename T::ToSetType::PositionType>;
  };

/*!
 * \brief A RelationLike type exposing the flattened storage used by RelationSet.
 *
 * RelationLike models row access. This refinement adds an additional static-storage
 * operations required to adapt a relation into a bivariate RelationSet.
 */
template <typename T>
concept FlatRelationLike = RelationLike<T> && Validatable<T> &&
  requires {
    typename T::FlatPositionType;
    typename T::RelationSubset;
  } && PositionLike<typename T::FlatPositionType> &&
  requires(const T& relation,
           typename T::FromSetType::PositionType fromPosition,
           typename T::FlatPositionType flatPosition) {
    { relation[fromPosition] } -> std::convertible_to<typename T::RelationSubset>;
    { relation.size(fromPosition) } -> PositionValueLike;
    { relation.offset(fromPosition) } -> std::convertible_to<typename T::FlatPositionType>;
    {
      static_cast<typename T::FromSetType::PositionType>(relation.firstIndex(flatPosition))
    } -> std::same_as<typename T::FromSetType::PositionType>;
    { relation.relationData().size() } -> std::convertible_to<typename T::FlatPositionType>;
    {
      relation.relationData()[flatPosition]
    } -> std::convertible_to<typename T::ToSetType::PositionType>;
  };

//------------------------------------------------------------------------------
// Maps
//------------------------------------------------------------------------------

/// \brief \a Value names the same underlying data type as \a Data, ignoring qualifiers.
template <typename Value, typename Data>
concept MapValueFor = std::same_as<std::remove_cvref_t<Value>, std::remove_cvref_t<Data>>;

/// \brief The size and element-access surface shared by univariate and bivariate maps.
template <typename T>
concept CommonMapModel =
  HasMapAssociatedTypes<T> && MapValueFor<typename T::ValueType, typename T::DataType> &&
  MapValueFor<typename T::ConstValueType, typename T::DataType> &&
  requires(T& map, const T& constMap, typename T::SetPosition pos) {
    { constMap.size() } -> std::same_as<typename T::SetPosition>;
    { map[pos] } -> std::same_as<typename T::ValueType>;
    { constMap[pos] } -> std::same_as<typename T::ConstValueType>;
  };

/// \brief A map whose domain is a univariate SetType.
template <typename T>
concept UnivariateMapLike = CommonMapModel<T> && HasUnivariateMapAssociatedTypes<T> &&
  UnivariateSetLike<typename T::SetType> &&
  std::same_as<typename T::SetPosition, typename T::SetType::PositionType> &&
  std::same_as<typename T::SetElement, typename T::SetType::ElementType> && requires(const T& map) {
    { map.set() } -> std::same_as<const typename T::SetType*>;
  };

/// \brief A map whose domain is a BivariateSetType.
template <typename T>
concept BivariateMapLike = CommonMapModel<T> && HasBivariateMapAssociatedTypes<T> &&
  BivariateSetLike<typename T::BivariateSetType> &&
  std::same_as<typename T::SetPosition, typename T::BivariateSetType::PositionType> &&
  std::same_as<typename T::SetElement, typename T::BivariateSetType::ElementType> &&
  requires(const T& map) {
    { map.set() } -> std::same_as<const typename T::BivariateSetType*>;
  };

/// \brief A univariate or bivariate SLAM map.
template <typename T>
concept MapLike = UnivariateMapLike<T> || BivariateMapLike<T>;

/// \brief A map whose semantic domain is exactly S.
template <typename M, typename S>
concept MapOver = (UnivariateMapLike<M> && std::same_as<typename M::SetType, S>) ||
  (BivariateMapLike<M> && std::same_as<typename M::BivariateSetType, S>);

//------------------------------------------------------------------------------
// Scalar value policies
//------------------------------------------------------------------------------

//---- base capabilities -------------------------------------------------------

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
  MapValueFor<typename T::IndirectionResult, Data> &&
  MapValueFor<typename T::ConstIndirectionResult, Data>;

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
concept HasResultPointerAliases = requires {
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
} &&
  requires(const typename T::IndirectionBufferType& buffer) {
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

/*!
 * \brief A MapLike type that support SubMap
 *
 * MapLike states the abstraction: a sized domain you can index.
 * A SubMap needs more than that, because it re-uses its storage policies.
 * It inherits the super-map's stride policy, derives its own reference type
 * from the super-map's indirection policy, and forwards range iteration.
 *
 * It relates to MapLike like FlatRelationLike relates to RelationLike, 
 * and like BivariateMapDomain relates to BivariateSetLike.
 */
template <typename T>
concept SubMappable = MapLike<T> && requires {
  typename T::StridePolicyType;
  typename T::IndirectionPolicy;
  typename T::range_iterator;
  typename T::const_range_iterator;
} && MapStridePolicyFor<typename T::StridePolicyType, typename T::SetPosition> &&
  MapIndirectionPolicyFor<typename T::IndirectionPolicy, typename T::SetPosition, typename T::DataType> &&
  requires(const T& map) {
    { map.shape() } -> std::same_as<typename T::StridePolicyType::ShapeType>;
  };

}  // namespace detail::model

//------------------------------------------------------------------------------
// Public concepts.
//
// Each alias normalizes the types that name a model (the container, policy,
// or position being checked) and forwards Element/Data unchanged,
// since their cv-qualification is part of the contract.
//------------------------------------------------------------------------------

/// \brief A built-in signed integral or explicitly opted-in SLAM position type.
/// \tparam T the candidate position type
template <typename T>
concept PositionLike = detail::model::PositionLike<detail::model_t<T>>;

/// \brief A univariate set with position-based size and element access.
/// \tparam T the candidate set type
/// \note BivariateSetLike refines this. Use UnivariateSetLike to exclude bivariate sets.
template <typename T>
concept SetLike = detail::model::SetLike<detail::model_t<T>>;

/// \brief A set that is not a bivariate set.
/// \tparam T the candidate set type
/// \note This is the one place the univariate/bivariate exclusion is stated.
template <typename T>
concept UnivariateSetLike = detail::model::UnivariateSetLike<detail::model_t<T>>;

/// \brief A SetLike type with const iteration over its elements.
/// \tparam T the candidate set type
template <typename T>
concept OrderedSetLike = detail::model::OrderedSetLike<detail::model_t<T>>;

/// \brief A set whose elements are coordinate pairs of positions from two component sets.
/// \tparam T the candidate bivariate set type
template <typename T>
concept BivariateSetLike = detail::model::BivariateSetLike<detail::model_t<T>>;

/// \brief A container that can check its own internal consistency.
/// \tparam T the candidate container type
template <typename T>
concept Validatable = detail::model::Validatable<detail::model_t<T>>;

/// \brief A BivariateSetLike type that a BivariateMap can bind field data over.
/// \tparam T the candidate bivariate set type
template <typename T>
concept BivariateMapDomain = detail::model::BivariateMapDomain<detail::model_t<T>>;

/// \brief A relation exposing its two sets and a const iterable row per from-set position.
/// \tparam T the candidate relation type
template <typename T>
concept RelationLike = detail::model::RelationLike<detail::model_t<T>>;

/// \brief A RelationLike type exposing the flattened storage used by RelationSet.
/// \tparam T the candidate relation type
template <typename T>
concept FlatRelationLike = detail::model::FlatRelationLike<detail::model_t<T>>;

/// \brief A map whose domain is a univariate SetType.
/// \tparam T the candidate map type
template <typename T>
concept UnivariateMapLike = detail::model::UnivariateMapLike<detail::model_t<T>>;

/// \brief A map whose domain is a BivariateSetType.
/// \tparam T the candidate map type
template <typename T>
concept BivariateMapLike = detail::model::BivariateMapLike<detail::model_t<T>>;

/// \brief A univariate or bivariate SLAM map.
/// \tparam T the candidate map type
template <typename T>
concept MapLike = detail::model::MapLike<detail::model_t<T>>;

/// \brief A subscriptable ordered set whose elements are flat positions in a larger space.
/// \tparam R the candidate range type
/// \tparam Position the flat position type its elements index
template <typename R, typename Position>
concept FlatRangeOver = detail::model::FlatRangeOver<R, detail::model_t<Position>>;

/// \brief A MapLike type that a SubMap can be taken of.
/// \tparam T the candidate super-map type
template <typename T>
concept SubMappable = detail::model::SubMappable<detail::model_t<T>>;

/// \brief A map whose semantic domain is exactly \a S.
/// \tparam M the candidate map type
/// \tparam S the expected domain: a SetLike for univariate maps, a BivariateSetLike otherwise
template <typename M, typename S>
concept MapOver = detail::model::MapOver<detail::model_t<M>, detail::model_t<S>>;

/// \brief A scalar runtime or compile-time value policy.
/// \tparam T the candidate policy type
template <typename T>
concept ValuePolicy = detail::model::ValuePolicy<detail::model_t<T>>;

/// \brief A policy that reports a size and whether that size is empty.
/// \tparam T the candidate policy type
template <typename T>
concept SizePolicy = detail::model::SizePolicy<detail::model_t<T>>;

/// \brief A SizePolicy usable by a set whose position type is \a Position.
/// \tparam T the candidate policy type
/// \tparam Position the owning set's position type
template <typename T, typename Position>
concept SetSizePolicyFor =
  detail::model::SetSizePolicyFor<detail::model_t<T>, detail::model_t<Position>>;

/// \brief The capability shared by scalar and multi-dimensional stride policies.
/// \tparam T the candidate policy type
template <typename T>
concept StridePolicy = detail::model::StridePolicy<detail::model_t<T>>;

/// \brief A scalar stride policy usable by OrderedSet with \a Position.
/// \tparam T the candidate policy type
/// \tparam Position the owning set's position type
template <typename T, typename Position>
concept OrderedSetStridePolicyFor =
  detail::model::OrderedSetStridePolicyFor<detail::model_t<T>, detail::model_t<Position>>;

/// \brief A scalar or multi-dimensional stride policy usable by Map with \a Position.
/// \tparam T the candidate policy type
/// \tparam Position the owning map's position type
template <typename T, typename Position>
concept MapStridePolicyFor =
  detail::model::MapStridePolicyFor<detail::model_t<T>, detail::model_t<Position>>;

/// \brief A scalar value policy that reports an offset.
/// \tparam T the candidate policy type
template <typename T>
concept OffsetPolicy = detail::model::OffsetPolicy<detail::model_t<T>>;

/// \brief An OffsetPolicy usable by an OrderedSet whose position type is \a Position.
/// \tparam T the candidate policy type
/// \tparam Position the owning set's position type
template <typename T, typename Position>
concept OrderedSetOffsetPolicyFor =
  detail::model::OrderedSetOffsetPolicyFor<detail::model_t<T>, detail::model_t<Position>>;

/// \brief The common storage/indirection-policy capability.
/// \tparam T the candidate policy type
template <typename T>
concept IndirectionPolicy = detail::model::IndirectionPolicy<detail::model_t<T>>;

/// \brief An IndirectionPolicy callable with \a Position.
/// \tparam T the candidate policy type
/// \tparam Position the position type the indirection is invoked with
template <typename T, typename Position>
concept IndirectionPolicyFor =
  detail::model::IndirectionPolicyFor<detail::model_t<T>, detail::model_t<Position>>;

/// \brief An indirection policy usable by OrderedSet over \a Position and \a Element.
/// \tparam T the candidate policy type
/// \tparam Position the owning set's position type
/// \tparam Element the owning set's element type; its cv-qualification is significant
template <typename T, typename Position, typename Element>
concept OrderedSetIndirectionPolicyFor =
  detail::model::OrderedSetIndirectionPolicyFor<detail::model_t<T>, detail::model_t<Position>, Element>;

/// \brief An indirection policy providing Map's buffer and static access API.
/// \tparam T the candidate policy type
/// \tparam Position the owning map's position type
/// \tparam Data the owning map's data type; its cv-qualification is significant
template <typename T, typename Position, typename Data>
concept MapIndirectionPolicyFor =
  detail::model::MapIndirectionPolicyFor<detail::model_t<T>, detail::model_t<Position>, Data>;

/// \brief A MapIndirectionPolicyFor that can allocate and initialize its buffer.
/// \tparam T the candidate policy type
/// \tparam Position the owning map's position type
/// \tparam Data the owning map's data type; its cv-qualification is significant
template <typename T, typename Position, typename Data>
concept AllocatingMapIndirectionPolicyFor =
  detail::model::AllocatingMapIndirectionPolicyFor<detail::model_t<T>, detail::model_t<Position>, Data>;

/// \brief A non-reference type that can be copied byte-for-byte into device code.
/// \tparam T the candidate type
template <typename T>
concept DeviceCapturable =
  !std::is_reference_v<T> && std::is_trivially_copyable_v<std::remove_cv_t<T>>;

namespace detail
{
//------------------------------------------------------------------------------
// Glue used by the construction helpers in SetBuilders, MapBuilders and
// RelationBuilders. These normalize like the public concepts above.
//------------------------------------------------------------------------------

/// \brief A non-Boolean integral or an opted-in value convertible to a position.
template <typename T>
concept PositionValueLike = model::PositionValueLike<T>;

/// \brief \a Position has at least as many value bits as \a RepresentedPosition.
template <typename Position, typename RepresentedPosition>
concept PositionCanRepresent =
  model::PositionCanRepresent<model_t<Position>, model_t<RepresentedPosition>>;

/// \brief \a Value is a position value convertible to \a Set's position type.
template <typename Set, typename Value>
concept SetPositionConvertible = UnivariateSetLike<Set> && PositionValueLike<Value> &&
  std::convertible_to<model_t<Value>, typename model_t<Set>::PositionType>;

/// \brief \a Position is exactly \a Set's position type.
template <typename Set, typename Position>
concept SetPositionSame =
  UnivariateSetLike<Set> && std::same_as<model_t<Position>, typename model_t<Set>::PositionType>;

/// \brief \a Position is \c void (meaning "unspecified") or exactly \a Set's position type.
template <typename Set, typename Position>
concept OptionalSetPositionSame = UnivariateSetLike<Set> &&
  (std::same_as<model_t<Position>, void> || SetPositionSame<Set, Position>);

/// \brief \a Stride is positive and representable by \a Position.
template <int Stride, typename Position>
concept PositiveStaticStrideForPosition =
  model::PositiveStaticStrideForPosition<Stride, model_t<Position>>;

/// \brief \a Stride is positive and representable by \a Set's position type.
template <int Stride, typename Set>
concept PositiveStaticStrideFor = model::PositiveStaticStrideFor<Stride, model_t<Set>>;
}  // namespace detail

}  // namespace axom::slam
