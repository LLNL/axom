// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file slam_concepts.cpp
 *
 * \brief Positive and negative compile-time tests for SLAM's C++20 concepts.
 */

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/slam/Aliases.hpp"
#include "axom/slam/BivariateMap.hpp"
#include "axom/slam/Concepts.hpp"
#include "axom/slam/DynamicConstantRelation.hpp"
#include "axom/slam/DynamicMap.hpp"
#include "axom/slam/DynamicVariableRelation.hpp"
#include "axom/slam/Map.hpp"
#include "axom/slam/ProductSet.hpp"
#include "axom/slam/RangeSet.hpp"
#include "axom/slam/RelationSet.hpp"
#include "axom/slam/RelationBuilders.hpp"
#include "axom/slam/Traits.hpp"
#include "axom/slam/policies/CardinalityPolicies.hpp"
#include "axom/slam/policies/IndirectionPolicies.hpp"
#include "axom/slam/policies/OffsetPolicies.hpp"
#include "axom/slam/policies/SizePolicies.hpp"
#include "axom/slam/policies/StridePolicies.hpp"

#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <vector>

namespace slam_concept_test
{
namespace slam = axom::slam;
namespace policies = axom::slam::policies;

using Position = slam::DefaultPositionType;
using Element = slam::DefaultElementType;
using Range = slam::RangeSet<Position, Element>;
using ConcreteRange = typename Range::ConcreteSet;
using Product = typename slam::ProductSet<ConcreteRange, ConcreteRange>::ConcreteSet;
using ArrayIndirection = policies::ArrayIndirection<Position, double>;
using ViewIndirection = policies::ArrayViewIndirection<Position, double>;
using ConstViewIndirection = policies::ArrayViewIndirection<Position, const double>;
using UnaryMap = slam::Map<double, ConcreteRange, ViewIndirection>;
using BinaryMap = slam::BivariateMap<double, Product, ViewIndirection>;
using DynamicSet = slam::DynamicSet<Position, Element>;
using DynamicMap = slam::DynamicMap<DynamicSet, double>;
using WrongDataIndirection = policies::ArrayViewIndirection<Position, int>;
using RawPointerMapIndirection = policies::CArrayIndirection<Position, double>;
using VectorMapIndirection = policies::STLVectorIndirection<Position, double>;
using VariableRelation = slam::VariableRelationView<ConcreteRange, ConcreteRange>;
using DynamicVariableRelation = slam::DynamicVariableRelation<ConcreteRange, ConcreteRange>;
using NarrowRange = slam::RangeSet<std::int32_t, std::int32_t>;
using WideRange = slam::RangeSet<std::int64_t, std::int64_t>;
using FirstHandleRange = slam::RangeSet<std::int32_t, double>;
using SecondHandleRange = slam::RangeSet<std::int32_t, float>;
using DistinctHandleProduct =
  typename slam::ProductSet<FirstHandleRange, SecondHandleRange>::ConcreteSet;
using HeterogeneousProduct = typename slam::ProductSet<NarrowRange, WideRange>::ConcreteSet;
using HeterogeneousNullBivariateSet = slam::NullBivariateSet<NarrowRange, WideRange>;
using HeterogeneousMapIndirection = policies::ArrayViewIndirection<std::int64_t, double>;
using HeterogeneousBinaryMap =
  slam::BivariateMap<double, HeterogeneousProduct, HeterogeneousMapIndirection>;
using HeterogeneousVariableRelation = slam::VariableRelationView<WideRange, NarrowRange>;
using HeterogeneousDynamicVariableRelation = slam::DynamicVariableRelation<NarrowRange, WideRange>;
using DynamicConstantCardinality =
  policies::ConstantCardinality<Position, policies::CompileTimeStride<Position, 3>>;
using DynamicConstantRelation =
  slam::DynamicConstantRelation<Position, Element, DynamicConstantCardinality>;
using DistinctElementDynamicConstantRelation =
  slam::DynamicConstantRelation<Position, std::int64_t, DynamicConstantCardinality>;

template <typename FirstSet, typename SecondSet>
concept CanFormProductSet = requires { sizeof(slam::ProductSet<FirstSet, SecondSet>); };

template <typename FirstSet, typename SecondSet, typename FlatPosition>
concept CanFormProductSetWithFlatPosition = requires {
  sizeof(slam::ProductSet<FirstSet, SecondSet, policies::ConcreteInterface, FlatPosition>);
};

template <typename RelationType>
concept CanFormRelationSet = slam::detail::RelationSetSource<RelationType> &&
  requires { sizeof(slam::RelationSet<RelationType>); };

struct StrongPosition
{
  std::int64_t value;
};

struct ExplicitFirstPosition
{
  explicit ExplicitFirstPosition(std::int32_t value = 0) : value(value) { }
  std::int32_t value;
};

struct ExplicitSecondPosition
{
  explicit ExplicitSecondPosition(std::int64_t value = 0) : value(value) { }
  std::int64_t value;
};

struct ExplicitFlatPosition
{
  explicit ExplicitFlatPosition(std::int64_t value = 0) : value(value) { }
  explicit ExplicitFlatPosition(ExplicitFirstPosition position) : value(position.value) { }
  explicit ExplicitFlatPosition(ExplicitSecondPosition position) : value(position.value) { }
  std::int64_t value;
};

struct ExplicitRowPosition
{
  explicit ExplicitRowPosition(std::int32_t value = 0) : value(value) { }
  std::int32_t value;
};

template <typename PositionType_, typename ElementType_>
struct MinimalSet
{
  using PositionType = PositionType_;
  using ElementType = ElementType_;

  PositionType size() const;
  bool empty() const;
  ElementType at(PositionType) const;
};

using ExplicitFirstSet = MinimalSet<ExplicitFirstPosition, double>;
using ExplicitSecondSet = MinimalSet<ExplicitSecondPosition, float>;

struct ExplicitPositionRow
{
  ExplicitRowPosition size() const;
  const ExplicitSecondPosition* begin() const;
  const ExplicitSecondPosition* end() const;
};

struct ExplicitPositionBivariateSet
{
  using FirstSetType = ExplicitFirstSet;
  using SecondSetType = ExplicitSecondSet;
  using PositionType = ExplicitFlatPosition;
  using ElementType = std::pair<ExplicitFirstPosition, ExplicitSecondPosition>;

  PositionType size() const;
  typename FirstSetType::PositionType firstSetSize() const;
  typename SecondSetType::PositionType secondSetSize() const;
  PositionType size(typename FirstSetType::PositionType) const;
  bool empty() const;
  ElementType at(PositionType) const;
  const FirstSetType* getFirstSet() const;
  const SecondSetType* getSecondSet() const;
  ExplicitPositionRow getElements(ExplicitFirstPosition) const;
  ExplicitRowPosition findElementIndex(ExplicitFirstPosition, ExplicitSecondPosition) const;
  PositionType findElementFlatIndex(ExplicitFirstPosition, ExplicitSecondPosition) const;
  ExplicitFirstPosition flatToFirstIndex(PositionType) const;
  ExplicitSecondPosition flatToSecondIndex(PositionType) const;
};

// Associated types alone must not satisfy a semantic concept.
struct TypedefOnlySet
{
  using PositionType = Position;
  using ElementType = Element;
};

struct WrongSizeSet
{
  using PositionType = Position;
  using ElementType = Element;

  double size() const;
  bool empty() const;
  Element at(Position) const;
};

struct MutableAccessOnlySet
{
  using PositionType = Position;
  using ElementType = Element;

  Position size() const;
  bool empty() const;
  Element at(Position);
};

struct FloatingPositionSet
{
  using PositionType = double;
  using ElementType = Element;

  double size() const;
  bool empty() const;
  Element at(double) const;
};

struct TypedefOnlyBivariateSet
{
  using FirstSetType = ConcreteRange;
  using SecondSetType = ConcreteRange;
  using PositionType = Position;
  using ElementType =
    std::pair<typename FirstSetType::PositionType, typename SecondSetType::PositionType>;
};

// Coordinate members use the positions of the component sets.
struct MinimalCoordinate
{
  NarrowRange::PositionType first;
  WideRange::PositionType second;
};

struct MinimalBivariateSet
{
  using FirstSetType = NarrowRange;
  using SecondSetType = WideRange;
  using PositionType = std::int64_t;
  using ElementType = MinimalCoordinate;

  PositionType size() const { return 3; }
  bool empty() const { return false; }
  ElementType at(PositionType pos) const { return {pos == 0 ? 0 : 1, pos}; }
  const FirstSetType* getFirstSet() const { return &first; }
  const SecondSetType* getSecondSet() const { return &second; }
  SecondSetType getElements(FirstSetType::PositionType pos) const
  {
    return pos == 0 ? SecondSetType(0, 1) : SecondSetType(1, 3);
  }

  FirstSetType first {2};
  SecondSetType second {3};
};

// A bivariate set is a set of coordinate pairs, so it models SetLike too.
struct BivariateSetLikeIsASet
{
  using FirstSetType = NarrowRange;
  using SecondSetType = WideRange;
  using PositionType = std::int64_t;
  using ElementType = MinimalCoordinate;

  PositionType size() const;
  typename FirstSetType::PositionType firstSetSize() const;
  typename SecondSetType::PositionType secondSetSize() const;
  PositionType size(typename FirstSetType::PositionType) const;
  bool empty() const;
  ElementType at(PositionType) const;
  const FirstSetType* getFirstSet() const;
  const SecondSetType* getSecondSet() const;
  SecondSetType getElements(FirstSetType::PositionType) const;
  PositionType findElementIndex(FirstSetType::PositionType, SecondSetType::PositionType) const;
  PositionType findElementFlatIndex(FirstSetType::PositionType, SecondSetType::PositionType) const;
  FirstSetType::PositionType flatToFirstIndex(PositionType) const;
  SecondSetType::PositionType flatToSecondIndex(PositionType) const;
};

// Negative control for the refinement: BivariateSetLike requires empty()
struct BivariateSetMissingEmpty
{
  using FirstSetType = NarrowRange;
  using SecondSetType = WideRange;
  using PositionType = std::int64_t;
  using ElementType = MinimalCoordinate;

  PositionType size() const;
  typename FirstSetType::PositionType firstSetSize() const;
  typename SecondSetType::PositionType secondSetSize() const;
  PositionType size(typename FirstSetType::PositionType) const;
  ElementType at(PositionType) const;
  const FirstSetType* getFirstSet() const;
  const SecondSetType* getSecondSet() const;
  SecondSetType getElements(FirstSetType::PositionType) const;
  PositionType findElementIndex(FirstSetType::PositionType, SecondSetType::PositionType) const;
  PositionType findElementFlatIndex(FirstSetType::PositionType, SecondSetType::PositionType) const;
  FirstSetType::PositionType flatToFirstIndex(PositionType) const;
  SecondSetType::PositionType flatToSecondIndex(PositionType) const;
};

struct HeterogeneousPositionBivariateSet
{
  using PositionType = typename NarrowRange::PositionType;
  using ElementType = typename WideRange::PositionType;
  using FirstSetType = NarrowRange;
  using SecondSetType = WideRange;

  PositionType size() const;
  typename FirstSetType::PositionType firstSetSize() const;
  typename SecondSetType::PositionType secondSetSize() const;
  PositionType size(typename FirstSetType::PositionType) const;
  bool empty() const;
  ElementType at(PositionType) const;
  const FirstSetType* getFirstSet() const;
  const SecondSetType* getSecondSet() const;
  PositionType findElementFlatIndex(PositionType, PositionType) const;
  PositionType flatToFirstIndex(PositionType) const;
  PositionType flatToSecondIndex(PositionType) const;
};

struct WrongElementBivariateSet
{
  using FirstSetType = NarrowRange;
  using SecondSetType = NarrowRange;
  using FirstPositionType = typename FirstSetType::PositionType;
  using SecondPositionType = typename SecondSetType::PositionType;
  using PositionType = std::int32_t;
  using ElementType = double;

  PositionType size() const;
  bool empty() const;
  PositionType size(FirstPositionType) const;
  ElementType at(PositionType) const;
  const FirstSetType* getFirstSet() const;
  const SecondSetType* getSecondSet() const;
  FirstPositionType firstSetSize() const;
  SecondPositionType secondSetSize() const;
  SecondSetType getElements(FirstPositionType) const;
  PositionType findElementIndex(FirstPositionType, SecondPositionType) const;
  PositionType findElementFlatIndex(FirstPositionType, SecondPositionType) const;
  FirstPositionType flatToFirstIndex(PositionType) const;
  SecondPositionType flatToSecondIndex(PositionType) const;
};

struct WrongCoordinateBivariateSet
{
  using FirstSetType = NarrowRange;
  using SecondSetType = WideRange;
  using FirstPositionType = typename FirstSetType::PositionType;
  using SecondPositionType = typename SecondSetType::PositionType;
  using PositionType = std::int64_t;
  using ElementType = std::pair<SecondPositionType, SecondPositionType>;

  PositionType size() const;
  bool empty() const;
  PositionType size(FirstPositionType) const;
  ElementType at(PositionType) const;
  const FirstSetType* getFirstSet() const;
  const SecondSetType* getSecondSet() const;
  FirstPositionType firstSetSize() const;
  SecondPositionType secondSetSize() const;
  SecondSetType getElements(FirstPositionType) const;
  PositionType findElementIndex(FirstPositionType, SecondPositionType) const;
  PositionType findElementFlatIndex(FirstPositionType, SecondPositionType) const;
  FirstPositionType flatToFirstIndex(PositionType) const;
  SecondPositionType flatToSecondIndex(PositionType) const;
};

struct TypedefOnlyRelation
{
  using FromSetType = ConcreteRange;
  using ToSetType = ConcreteRange;
};

struct MinimalRelationRow
{
  std::size_t size() const { return count; }
  const WideRange::PositionType* begin() const { return data; }
  const WideRange::PositionType* end() const { return data + count; }

  const WideRange::PositionType* data;
  std::size_t count;
};

struct MinimalRelation
{
  using FromSetType = NarrowRange;
  using ToSetType = WideRange;

  const FromSetType* fromSet() const { return &from; }
  const ToSetType* toSet() const { return &to; }
  MinimalRelationRow operator[](FromSetType::PositionType pos) const
  {
    return pos == 0 ? MinimalRelationRow {indices, 1} : MinimalRelationRow {indices + 1, 2};
  }

  FromSetType from {2};
  ToSetType to {3};
  WideRange::PositionType indices[3] {2, 0, 1};
};

struct NotAPosition
{ };

struct WrongRelationRow
{
  std::size_t size() const;
  NotAPosition* begin() const;
  NotAPosition* end() const;
};

struct WrongRelationEntry
{
  using FromSetType = ConcreteRange;
  using ToSetType = WideRange;

  const FromSetType* fromSet() const;
  const ToSetType* toSet() const;
  WrongRelationRow operator[](FromSetType::PositionType) const;
};

struct ConvertibleRelationRow
{
  using RelationSubset = typename VariableRelation::RelationSubset;

  VariableRelation::FlatPositionType size() const;
  typename RelationSubset::const_iterator begin() const;
  typename RelationSubset::const_iterator end() const;
  operator RelationSubset() const;
};

struct ProxyRowRelation : VariableRelation
{
  ConvertibleRelationRow operator[](FromSetType::PositionType) const;
};

struct ExplicitFirstIndexResult
{
  using FromPositionType = typename VariableRelation::FromPositionType;
  explicit operator FromPositionType() const;
};

struct ExplicitFirstIndexRelation : VariableRelation
{
  ExplicitFirstIndexResult firstIndex(FlatPositionType) const;
};

struct TypedefOnlyMap
{
  using DataType = double;
  using SetType = ConcreteRange;
  using PositionType = Position;
  using SetElement = Element;
  using ValueType = double&;
  using ConstValueType = const double&;
};

// A map needs no SLAM storage aliases or bound set object.
struct MinimalBivariateMap
{
  using DataType = double;
  using PositionType = std::int64_t;

  PositionType size() const { return 2; }
  PositionType numComp() const { return 3; }
  MinimalCoordinate index(PositionType pos) const { return {1, 2 * pos}; }
  double& flatValue(PositionType pos, PositionType component)
  {
    return values[3 * pos + component];
  }
  const double& flatValue(PositionType pos, PositionType component) const
  {
    return values[3 * pos + component];
  }

  double values[6] {};
};

struct WrongDomainMap : MinimalBivariateMap
{
  using MappedSetType = MinimalBivariateSet;
  const TypedefOnlySet* set() const;
};

struct WrongPositionMap : MinimalBivariateMap
{
  using PositionType = short;
};

struct WrongValueTypeMap : MinimalBivariateMap
{
  int& flatValue(PositionType, PositionType);
  const int& flatValue(PositionType, PositionType) const;
};

struct WrongMutableAccessMap : MinimalBivariateMap
{
  double flatValue(PositionType, PositionType);
  const double& flatValue(PositionType, PositionType) const;
};

struct DropsDataConstnessMap : MinimalBivariateMap
{
  using DataType = const double;
};

struct MissingComponentCountMap : MinimalBivariateMap
{
  void numComp() const = delete;
};

struct MissingIndexMap : MinimalBivariateMap
{
  void index(PositionType) const;
};

struct MissingConstAccessMap : MinimalBivariateMap
{
  using MinimalBivariateMap::flatValue;
  const double& flatValue(PositionType, PositionType) const = delete;
};

// Unused aliases must not change the classification of a working map.
struct ExtraAliasesMap : MinimalBivariateMap
{
  using SetType = TypedefOnlySet;
  using BivariateSetType = void;
  using ValueType = int;
};

struct TypedefOnlyValuePolicy
{
  struct TagType;
  using IntType = int;
};

struct WrongValuePolicy
{
  struct TagType;
  using IntType = int;

  double value() const;
  bool isValid(bool) const;
};

struct TypedefOnlyIndirection
{
  using IndirectionResult = double&;
  using ConstIndirectionResult = const double&;
  using IndirectionBufferType = double*;
  using IndirectionPtrType = double*;

  static constexpr bool DeviceAccessible = true;
};

struct PrvalueMapIndirection : ViewIndirection
{
  using IndirectionResult = double;
  using ConstIndirectionResult = const double&;
  using ResultPtr = double*;
  using ConstResultPtr = const double*;
};

struct MismatchedConstPointerIndirection : ViewIndirection
{
  using IndirectionResult = double&;
  using ConstIndirectionResult = double&;
  using ResultPtr = double*;
  using ConstResultPtr = const double*;
};

// Access still returns mutable data despite the declared const element type.
struct DropsReferentConst : ViewIndirection
{
  using ElementType = const double;
};

struct AccessDropsReferentConst : ConstViewIndirection
{
  double& indirection(Position);
  double& indirection(Position) const;
};

struct NonDefaultBuffer
{
  NonDefaultBuffer() = delete;
  explicit NonDefaultBuffer(Position size) : m_size(size) { }

  Position size() const { return m_size; }
  bool empty() const { return m_size == 0; }
  void resize(Position size) { m_size = size; }

  Position m_size;
};

struct NonDefaultBufferIndirection
{
  using PositionType = Position;
  using ElementType = double;
  using IndirectionResult = double&;
  using ConstIndirectionResult = const double&;
  using IndirectionBufferType = NonDefaultBuffer;
  using IndirectionPtrType = IndirectionBufferType*;
  using IndirectionRefType = IndirectionBufferType&;
  using IndirectionConstRefType = const IndirectionBufferType&;
  using ResultPtr = double*;
  using ConstResultPtr = const double*;

  static constexpr bool DeviceAccessible = false;
  static constexpr bool IsMutableBuffer = true;

  bool hasIndirection() const { return true; }

  static ResultPtr getIndirection(IndirectionBufferType&, Position = 0) { return nullptr; }
  static ConstResultPtr getConstIndirection(const IndirectionBufferType&, Position = 0)
  {
    return nullptr;
  }
  static IndirectionBufferType create(Position size, const double&, int)
  {
    return IndirectionBufferType(size);
  }
};

// Supplies everything the Map family consumes, but omits Indirection[Const]RefType
struct LeanMapIndirection
{
  using PositionType = Position;
  using ElementType = double;
  using IndirectionResult = double&;
  using ConstIndirectionResult = const double&;
  using IndirectionBufferType = std::vector<double>;
  using IndirectionPtrType = IndirectionBufferType*;
  using ResultPtr = double*;
  using ConstResultPtr = const double*;

  static constexpr bool DeviceAccessible = false;
  static constexpr bool IsMutableBuffer = true;

  bool hasIndirection() const { return true; }

  static ResultPtr getIndirection(IndirectionBufferType&, Position = 0) { return nullptr; }
  static ConstResultPtr getConstIndirection(const IndirectionBufferType&, Position = 0)
  {
    return nullptr;
  }
};

// Positioned access only. Map::data_ptr() needs the whole-buffer form.
struct PositionedAccessOnlyIndirection : LeanMapIndirection
{
  static ResultPtr getIndirection(IndirectionBufferType&, Position) { return nullptr; }
  static ConstResultPtr getConstIndirection(const IndirectionBufferType&, Position)
  {
    return nullptr;
  }
};

struct TrivialCapture
{
  int value;
};

struct NonTrivialCapture
{
  NonTrivialCapture(const NonTrivialCapture&) { }
};

}  // namespace slam_concept_test

namespace slam_concept_test
{
// Sets and bivariate sets
static_assert(slam::SetLike<slam::Set<>>);
static_assert(slam::SetLike<Range>);
static_assert(slam::SetLike<const Range&>);
static_assert(slam::IterableSetLike<Range>);
static_assert(slam::IterableSetLike<const Range&>);
static_assert(!slam::BivariateSetLike<Range>);
static_assert(slam::BivariateSetLike<Product>);
static_assert(slam::BivariateSetLike<const Product&>);
static_assert(std::same_as<typename Product::ElementType, std::pair<Position, Position>>);
static_assert(slam::BivariateSetLike<DistinctHandleProduct>);
static_assert(
  std::same_as<typename DistinctHandleProduct::ElementType, std::pair<std::int32_t, std::int32_t>>);
static_assert(slam::BivariateSetLike<HeterogeneousProduct>);
static_assert(slam::BivariateSetLike<HeterogeneousNullBivariateSet>);
static_assert(slam::BivariateSetLike<MinimalBivariateSet>);
static_assert(
  !std::same_as<typename MinimalBivariateSet::ElementType, std::pair<std::int32_t, std::int64_t>>);
static_assert(!slam::BivariateSetLike<ExplicitPositionBivariateSet>);
static_assert(!std::convertible_to<ExplicitFirstPosition, ExplicitFlatPosition>);
static_assert(!std::convertible_to<ExplicitSecondPosition, ExplicitFlatPosition>);
static_assert(!std::convertible_to<ExplicitRowPosition, ExplicitFlatPosition>);
static_assert(std::same_as<typename HeterogeneousProduct::FirstPositionType, std::int32_t>);
static_assert(std::same_as<typename HeterogeneousProduct::SecondPositionType, std::int64_t>);
static_assert(std::same_as<typename HeterogeneousProduct::PositionType, std::int64_t>);
static_assert(
  std::same_as<typename HeterogeneousProduct::ElementType, std::pair<std::int32_t, std::int64_t>>);
// BivariateSetLike refines SetLike: a bivariate set is a set of coordinates.
static_assert(slam::SetLike<Product>);
static_assert(slam::SetLike<ConcreteRange>);
static_assert(slam::BivariateSetLike<BivariateSetLikeIsASet>);
static_assert(slam::SetLike<BivariateSetLikeIsASet>);
static_assert(slam::SetLike<MinimalBivariateSet>);
static_assert(!slam::BivariateSetLike<BivariateSetMissingEmpty>);
static_assert(!slam::SetLike<BivariateSetMissingEmpty>);

// detail::BivariateMapSet is the contract BivariateMap relies on.
static_assert(slam::detail::BivariateMapSet<Product>);
static_assert(slam::detail::BivariateMapSet<HeterogeneousProduct>);
static_assert(!slam::detail::BivariateMapSet<MinimalBivariateSet>,
              "coordinate structure alone is not enough to bind a BivariateMap");
static_assert(slam::Validatable<Product>);
static_assert(!slam::Validatable<MinimalBivariateSet>);

// detail::RelationSetSource gained isValid(), which RelationSet::isValid() forwards to.
static_assert(slam::Validatable<VariableRelation>);

template <typename S>
  requires slam::SetLike<S>
std::integral_constant<int, 1> selectByConstraint();
template <typename S>
  requires slam::BivariateSetLike<S>
std::integral_constant<int, 2> selectByConstraint();

static_assert(decltype(selectByConstraint<ConcreteRange>())::value == 1);
static_assert(decltype(selectByConstraint<Product>())::value == 2);
static_assert(decltype(selectByConstraint<MinimalBivariateSet>())::value == 2);
static_assert(!slam::BivariateSetLike<ConcreteRange>);
static_assert(!slam::SetLike<TypedefOnlySet>);
static_assert(!slam::SetLike<WrongSizeSet>);
static_assert(!slam::SetLike<MutableAccessOnlySet>);
static_assert(!slam::SetLike<FloatingPositionSet>);
static_assert(!slam::BivariateSetLike<TypedefOnlyBivariateSet>);
static_assert(!slam::BivariateSetLike<HeterogeneousPositionBivariateSet>);
static_assert(!slam::BivariateSetLike<WrongElementBivariateSet>);
static_assert(!slam::BivariateSetLike<WrongCoordinateBivariateSet>);
static_assert(CanFormProductSet<NarrowRange, NarrowRange>);
static_assert(CanFormProductSet<NarrowRange, WideRange>);
static_assert(CanFormProductSetWithFlatPosition<NarrowRange, WideRange, std::int64_t>);
static_assert(!CanFormProductSetWithFlatPosition<NarrowRange, WideRange, std::int32_t>);
static_assert(!CanFormProductSetWithFlatPosition<NarrowRange, NarrowRange, std::uint64_t>);
static_assert(!CanFormProductSet<TypedefOnlySet, NarrowRange>);
static_assert(!slam::SetLike<int>);

// Relations
static_assert(slam::RelationLike<VariableRelation>);
static_assert(slam::RelationLike<const VariableRelation&>);
static_assert(slam::RelationLike<HeterogeneousVariableRelation>);
static_assert(slam::RelationLike<DynamicVariableRelation>);
static_assert(slam::RelationLike<HeterogeneousDynamicVariableRelation>);
static_assert(std::same_as<typename HeterogeneousDynamicVariableRelation::FromPositionType,
                           typename NarrowRange::PositionType>);
static_assert(std::same_as<typename HeterogeneousDynamicVariableRelation::ToPositionType,
                           typename WideRange::PositionType>);
static_assert(CanFormRelationSet<VariableRelation>);
static_assert(CanFormRelationSet<HeterogeneousVariableRelation>);
static_assert(std::constructible_from<slam::RelationSet<VariableRelation>, VariableRelation*>,
              "an accepted flat relation must instantiate its RelationSet consumer");
static_assert(
  std::constructible_from<slam::RelationSet<HeterogeneousVariableRelation>, HeterogeneousVariableRelation*>,
  "heterogeneous from-set and to-set positions must instantiate their RelationSet consumer");
static_assert(slam::detail::RelationSetSource<VariableRelation>);
static_assert(slam::detail::RelationSetSource<const VariableRelation&>);
static_assert(slam::detail::RelationSetSource<ProxyRowRelation>);
static_assert(slam::detail::RelationSetSource<ExplicitFirstIndexRelation>);
static_assert(slam::RelationLike<DynamicConstantRelation>);
static_assert(slam::RelationLike<DistinctElementDynamicConstantRelation>);
static_assert(slam::RelationLike<const DynamicConstantRelation&>);
static_assert(slam::is_relation_like_v<DynamicConstantRelation>);
static_assert(slam::RelationLike<MinimalRelation>);
static_assert(!slam::detail::RelationSetSource<MinimalRelation>);
static_assert(!slam::detail::RelationSetSource<DynamicVariableRelation>);
static_assert(!slam::detail::RelationSetSource<DynamicConstantRelation>);
static_assert(!CanFormRelationSet<MinimalRelation>);
static_assert(!CanFormRelationSet<DynamicVariableRelation>);
static_assert(!CanFormRelationSet<DynamicConstantRelation>);
static_assert(!slam::RelationLike<TypedefOnlyRelation>);
static_assert(!slam::RelationLike<WrongRelationEntry>);

// Maps
// BivariateMap::SetType is its flat backing set.
// MappedSetType names the bivariate set that MapOver checks.
static_assert(slam::MapLike<UnaryMap>);
static_assert(slam::MapLike<BinaryMap>);
static_assert(slam::MapLike<HeterogeneousBinaryMap>);
static_assert(slam::MapLike<MinimalBivariateMap>);
static_assert(std::same_as<typename HeterogeneousBinaryMap::FirstPositionType, std::int32_t>);
static_assert(std::same_as<typename HeterogeneousBinaryMap::SecondPositionType, std::int64_t>);
static_assert(std::same_as<decltype(HeterogeneousBinaryMap::iterator {}.firstIndex()), std::int32_t>);
static_assert(std::same_as<decltype(HeterogeneousBinaryMap::iterator {}.secondIndex()), std::int64_t>);
static_assert(
  std::same_as<decltype(HeterogeneousBinaryMap::range_iterator {}.firstIndex()), std::int32_t>);
static_assert(
  std::same_as<decltype(HeterogeneousBinaryMap::range_iterator {}.secondIndex()), std::int64_t>);
static_assert(
  std::same_as<typename HeterogeneousBinaryMap::SetElement, std::pair<std::int32_t, std::int64_t>>);
static_assert(slam::MapLike<const UnaryMap&>);
static_assert(slam::MapLike<DynamicMap>);
static_assert(slam::MapLike<const DynamicMap&>);
static_assert(slam::MapOver<UnaryMap, ConcreteRange>);
static_assert(slam::MapOver<BinaryMap, Product>);
static_assert(slam::MapOver<const BinaryMap&, const Product&>);
static_assert(slam::MapOver<DynamicMap, DynamicSet>);
static_assert(!slam::MapOver<BinaryMap, typename BinaryMap::SetType>);
static_assert(!slam::MapLike<TypedefOnlyMap>);
static_assert(slam::MapLike<WrongDomainMap>);
static_assert(!slam::MapOver<WrongDomainMap, MinimalBivariateSet>);
static_assert(slam::MapLike<ExtraAliasesMap>);
static_assert(!slam::MapOver<MinimalBivariateMap, MinimalBivariateSet>);
static_assert(!slam::MapLike<DropsDataConstnessMap>);
static_assert(!slam::MapLike<MissingComponentCountMap>);
static_assert(!slam::MapLike<MissingIndexMap>);
static_assert(!slam::MapLike<MissingConstAccessMap>);
static_assert(slam::MapLike<const MinimalBivariateMap&>);
static_assert(!slam::MapLike<WrongPositionMap>);
static_assert(!slam::MapLike<WrongValueTypeMap>);
static_assert(!slam::MapLike<WrongMutableAccessMap>);
static_assert(!slam::MapLike<int>);
static_assert(!slam::MapOver<int, ConcreteRange>);
using UnarySubMap = slam::SubMap<UnaryMap, ConcreteRange>;
using ConstUnarySubMap = slam::SubMap<const UnaryMap, ConcreteRange>;
using BinarySubMap = typename BinaryMap::SubMapType;
using ConstBinarySubMap = typename BinaryMap::ConstSubMapType;
static_assert(slam::MapLike<UnarySubMap>);
static_assert(slam::MapLike<ConstUnarySubMap>);
static_assert(!slam::MapOver<UnarySubMap, ConcreteRange>);
static_assert(!slam::MapOver<ConstUnarySubMap, ConcreteRange>);
static_assert(slam::MapLike<BinarySubMap>);
static_assert(slam::MapLike<ConstBinarySubMap>);
static_assert(!slam::MapOver<BinarySubMap, typename BinaryMap::SetType>);

// SubMap consumes parent operations, not parent storage policies.
static_assert(slam::detail::SubMapSource<UnaryMap>);
static_assert(slam::detail::SubMapSource<BinaryMap>);
static_assert(slam::detail::SubMapSource<const UnaryMap>);
static_assert(!slam::detail::SubMapSource<ConcreteRange>);
static_assert(!slam::detail::SubMapSource<int>);
// A SubMap can itself serve as another SubMap's parent.
static_assert(slam::detail::SubMapSource<BinarySubMap>);
static_assert(
  slam::detail::FlatRangeOver<typename BinarySubMap::IndexSetType, typename BinarySubMap::PositionType>);
using NestedSubMap = slam::SubMap<BinarySubMap, typename BinarySubMap::IndexSetType>;
static_assert(slam::MapLike<NestedSubMap>);
static_assert(slam::detail::SubMapSource<NestedSubMap>, "and it composes to any depth");
static_assert(std::same_as<typename NestedSubMap::ProjectedElement, Product::ElementType>);
static_assert(
  std::same_as<decltype(std::declval<const BinaryMap&>().index(Position {})), Product::ElementType>);

// MapLike is agnostic about whether const access is deep or shallow.
// and SLAM has both use-cases.
using DeepMap = slam::Map<double, ConcreteRange, ArrayIndirection>;
using ShallowMap = slam::Map<double, ConcreteRange, ViewIndirection>;
static_assert(slam::MapLike<DeepMap> && slam::MapLike<ShallowMap> && slam::MapLike<BinarySubMap>);
static_assert(std::is_same_v<typename DeepMap::ConstValueType, const double&>,
              "an owning indirection is deep-const");
static_assert(std::is_same_v<typename ShallowMap::ConstValueType, double&>,
              "the same Map over a view indirection is already shallow-const");
static_assert(std::is_same_v<typename BinarySubMap::ValueType,
                             typename BinarySubMap::ConstValueType>,
              "a SubMap is a view: constness rides on SuperMapType, not on the object");
// The two axes are independent. The yielded reference is decided by whether SuperMapType is const,
// while the super-map's indirection policy determines how deep that const goes.
using DeepBinaryMap = slam::BivariateMap<double, Product, ArrayIndirection>;
static_assert(std::is_same_v<
                typename std::remove_const_t<typename DeepBinaryMap::ConstSubMapType>::ConstValueType,
                const double&>,
              "over a deep-const super-map, a SubMap does yield const references");
static_assert(std::is_same_v<
                typename std::remove_const_t<typename BinaryMap::ConstSubMapType>::ConstValueType,
                double&>,
              "over a view-backed super-map it stays shallow, as that policy dictates");

// detail::FlatRangeOver is the index set a SubMap is built over, and the row range reported by a bivariate set.
// Subscript lives here rather than in IterableSetLike because a bivariate set is IterableSetLike
// and has no operator[].
static_assert(
  slam::detail::FlatRangeOver<typename BinaryMap::SetType, typename BinaryMap::PositionType>);
static_assert(slam::IterableSetLike<Product>);
static_assert(!slam::detail::FlatRangeOver<Product, typename Product::PositionType>,
              "a bivariate set is ordered but not subscriptable");
static_assert(!slam::MapOver<ConstBinarySubMap, typename BinaryMap::SetType>);
static_assert(std::same_as<typename BinarySubMap::IndexSetType, typename BinaryMap::SetType>);
static_assert(std::same_as<typename BinarySubMap::ProjectedElement, typename Product::ElementType>);
static_assert(std::same_as<decltype(std::declval<const BinarySubMap&>().index(Position {})),
                           typename Product::ElementType>);

// Policies
using Size = policies::CompileTimeSize<int, 5>;
using EmptySize = policies::ZeroSize<int>;
using RuntimeSize = policies::RuntimeSize<int>;
using WrongRuntimeSize = policies::RuntimeSize<double>;
using WrongCompileTimeSize = policies::CompileTimeSize<std::int64_t, 5>;
using WrongEmptySize = policies::ZeroSize<std::int64_t>;
using ScalarStride = policies::CompileTimeStride<int, 3>;
using MatrixStride = policies::MultiDimStride<int, 2>;
using Offset = policies::CompileTimeOffset<int, 4>;
using RuntimeOffset = policies::RuntimeOffset<int>;
using EmptyOffset = policies::ZeroOffset<int>;
using WrongRuntimeOffset = policies::RuntimeOffset<std::int64_t>;
using WrongCompileTimeOffset = policies::CompileTimeOffset<std::int64_t, 4>;
using WrongEmptyOffset = policies::ZeroOffset<std::int64_t>;
using NoIndirection = policies::NoIndirection<Position, Element>;
using OwningIndirection = policies::ArrayIndirection<Position, double>;
using OwningMap = slam::Map<double, ConcreteRange, OwningIndirection>;

static_assert(slam::detail::ValuePolicy<Size>);
static_assert(slam::SizePolicy<Size>);
static_assert(slam::SizePolicy<EmptySize>);
static_assert(!slam::detail::ValuePolicy<EmptySize>);
static_assert(slam::detail::SetSizePolicyFor<RuntimeSize, int>);
static_assert(slam::detail::SetSizePolicyFor<Size, int>);
static_assert(slam::detail::SetSizePolicyFor<EmptySize, int>);
static_assert(slam::SizePolicy<WrongRuntimeSize>);
static_assert(!slam::detail::SetSizePolicyFor<WrongRuntimeSize, int>);
static_assert(!slam::detail::SetSizePolicyFor<WrongCompileTimeSize, int>);
static_assert(!slam::detail::SetSizePolicyFor<WrongEmptySize, int>);
static_assert(slam::StridePolicy<ScalarStride>);
static_assert(slam::StridePolicy<MatrixStride>);
static_assert(slam::detail::OrderedSetStridePolicyFor<ScalarStride, int>);
static_assert(!slam::detail::OrderedSetStridePolicyFor<MatrixStride, int>);
static_assert(slam::detail::MapStridePolicyFor<ScalarStride, int>);
static_assert(slam::detail::MapStridePolicyFor<MatrixStride, int>);
// Maps support a stride index that converts to their (possibly wider) position type
// OrderedSet's scalar value policy must use that exact position type.
#if !defined(AXOM_NO_INT64_T)
// Use explicit types for wide/narrow pairs so they do not depend on AXOM_USE_64BIT_INDEXTYPE.
using WidePosition = std::int64_t;
using NarrowStride = policies::CompileTimeStride<std::int32_t, 3>;
static_assert(!std::is_same_v<NarrowStride::IndexType, WidePosition>,
              "the wide and narrow types must actually differ for this to test anything");
static_assert(slam::detail::MapStridePolicyFor<NarrowStride, WidePosition>);
static_assert(!slam::detail::OrderedSetStridePolicyFor<NarrowStride, WidePosition>);
#endif
// The substitutability checks behind the *For concepts, asserted directly so a
// failure names the clause rather than an unnamed conjunction.
namespace model = slam::detail;
static_assert(model::PolicyDefaultedOver<RuntimeSize, int>);
static_assert(!model::PolicyDefaultedOver<WrongRuntimeSize, int>);
static_assert(model::ScalarValuePolicyOver<ScalarStride, int>);
static_assert(!model::ScalarValuePolicyOver<EmptySize, int>, "a size policy carries no IntType");
static_assert(model::ScalarStridePolicyOver<ScalarStride, int>);
static_assert(!model::ScalarStridePolicyOver<MatrixStride, int>, "multi-dim stride is not scalar");
static_assert(!model::ExposesPerDimensionStrides<ScalarStride>);
static_assert(model::ExposesPerDimensionStrides<MatrixStride>);

static_assert(slam::OffsetPolicy<Offset>);
static_assert(slam::detail::OrderedSetOffsetPolicyFor<RuntimeOffset, int>);
static_assert(slam::detail::OrderedSetOffsetPolicyFor<Offset, int>);
static_assert(slam::detail::OrderedSetOffsetPolicyFor<EmptyOffset, int>);
static_assert(!slam::detail::OrderedSetOffsetPolicyFor<WrongRuntimeOffset, int>);
static_assert(!slam::detail::OrderedSetOffsetPolicyFor<WrongCompileTimeOffset, int>);
static_assert(!slam::detail::OrderedSetOffsetPolicyFor<WrongEmptyOffset, int>);
static_assert(!slam::detail::ValuePolicy<TypedefOnlyValuePolicy>);
static_assert(!slam::detail::ValuePolicy<WrongValuePolicy>);
static_assert(!slam::SizePolicy<int>);
static_assert(!slam::StridePolicy<int>);
static_assert(!slam::OffsetPolicy<int>);
static_assert(slam::detail::IndirectionPolicy<ViewIndirection>);
static_assert(slam::detail::IndirectionPolicyFor<ViewIndirection, Position>);
static_assert(slam::OrderedSetIndirectionPolicyFor<ViewIndirection, Position, double>);
// The set and map indirection concepts agree on cv-qualification:
// the policy's ElementType must match the container's element/data type exactly.
static_assert(slam::OrderedSetIndirectionPolicyFor<ConstViewIndirection, Position, const double>);
static_assert(!slam::OrderedSetIndirectionPolicyFor<ViewIndirection, Position, const double>);
static_assert(!slam::OrderedSetIndirectionPolicyFor<ConstViewIndirection, Position, double>);
static_assert(slam::MapIndirectionPolicyFor<ViewIndirection, Position, double>);
static_assert(slam::MapIndirectionPolicyFor<ConstViewIndirection, Position, const double>);
static_assert(!slam::MapIndirectionPolicyFor<ViewIndirection, Position, const double>);
static_assert(!slam::MapIndirectionPolicyFor<ConstViewIndirection, Position, double>);
static_assert(!slam::MapIndirectionPolicyFor<DropsReferentConst, Position, const double>);
static_assert(!slam::OrderedSetIndirectionPolicyFor<DropsReferentConst, Position, const double>);
static_assert(!slam::OrderedSetIndirectionPolicyFor<AccessDropsReferentConst, Position, const double>);
static_assert(!slam::AllocatingMapIndirectionPolicyFor<ViewIndirection, Position, double>);
static_assert(!slam::MapIndirectionPolicyFor<RawPointerMapIndirection, Position, double>);
static_assert(!slam::MapIndirectionPolicyFor<PrvalueMapIndirection, Position, double>);
static_assert(!slam::MapIndirectionPolicyFor<MismatchedConstPointerIndirection, Position, double>);
static_assert(slam::MapIndirectionPolicyFor<VectorMapIndirection, Position, double>);
static_assert(slam::AllocatingMapIndirectionPolicyFor<VectorMapIndirection, Position, double>);
static_assert(slam::detail::IndirectionPolicy<NoIndirection>);
static_assert(slam::detail::IndirectionPolicyFor<NoIndirection, Position>);
static_assert(slam::OrderedSetIndirectionPolicyFor<NoIndirection, Position, Element>);
static_assert(!slam::MapIndirectionPolicyFor<NoIndirection, Position, Element>);
static_assert(slam::MapIndirectionPolicyFor<OwningIndirection, Position, double>);
static_assert(slam::AllocatingMapIndirectionPolicyFor<OwningIndirection, Position, double>);
static_assert(slam::AllocatingMapIndirectionPolicyFor<NonDefaultBufferIndirection, Position, double>);
using NonDefaultBufferMap = slam::Map<double, ConcreteRange, NonDefaultBufferIndirection>;
static_assert(std::constructible_from<NonDefaultBufferMap, const ConcreteRange*>);
static_assert(std::default_initializable<OwningMap>);
static_assert(!std::default_initializable<UnaryMap>);
static_assert(std::constructible_from<UnaryMap, const ConcreteRange*, typename UnaryMap::OrderedMap>);
static_assert(!slam::MapIndirectionPolicyFor<WrongDataIndirection, Position, double>);

// The Map indirection contract that Map family reads.
// Indirection[Const]RefType are StaticRelation's aliases
static_assert(slam::MapIndirectionPolicyFor<LeanMapIndirection, Position, double>);
// ... but the whole-buffer accessor behind Map::data_ptr() is required.
static_assert(!slam::MapIndirectionPolicyFor<PositionedAccessOnlyIndirection, Position, double>);
static_assert(model::HasStaticBufferAccess<LeanMapIndirection, Position>);
static_assert(!model::HasStaticBufferAccess<PositionedAccessOnlyIndirection, Position>);
// Map's data_ptr() (private) is declared in terms of the policy's ResultPtr alias
static_assert(std::same_as<typename ViewIndirection::ResultPtr, double*>);
static_assert(std::same_as<typename ViewIndirection::ConstResultPtr, double*>,
              "ArrayView indirection has shallow constness");

// The atoms behind the two indirection *For concepts.
static_assert(model::IndirectsExactly<ViewIndirection, double>);
static_assert(!model::IndirectsExactly<ViewIndirection, const double>);
static_assert(model::IndirectsExactly<ConstViewIndirection, const double>);
static_assert(model::YieldsStableReferences<ViewIndirection>);
static_assert(!model::YieldsStableReferences<PrvalueMapIndirection>);
static_assert(model::HasResultPointerAliases<ViewIndirection>);
static_assert(!model::HasResultPointerAliases<MismatchedConstPointerIndirection>);
static_assert(model::HasSizedBuffer<VectorMapIndirection, Position>);
static_assert(model::BindableIndirection<ViewIndirection>);
static_assert(model::ValidatesSetRange<ViewIndirection, Position>);
static_assert(!slam::detail::IndirectionPolicy<TypedefOnlyIndirection>);
static_assert(!slam::detail::IndirectionPolicy<int>);

// Position types and representation properties
static_assert(slam::PositionLike<int>);
static_assert(!slam::PositionLike<unsigned>);
static_assert(!slam::PositionLike<bool>);
static_assert(!slam::PositionLike<double>);
static_assert(!slam::PositionLike<StrongPosition>);
static_assert(!slam::PositionLike<ExplicitFirstPosition>);
static_assert(!slam::PositionLike<ExplicitSecondPosition>);
static_assert(!slam::PositionLike<ExplicitFlatPosition>);
static_assert(!slam::PositionLike<ExplicitRowPosition>);
static_assert(slam::PositionLike<const std::int64_t&>);
static_assert(slam::TriviallyCopyableRepresentation<TrivialCapture>);
static_assert(!slam::TriviallyCopyableRepresentation<NonTrivialCapture>);
static_assert(!slam::TriviallyCopyableRepresentation<TrivialCapture&>);
// These host-object views are trivially copyable. That says nothing about device use.
struct DerivedSubMap : UnarySubMap
{ };
static_assert(std::is_trivially_copyable_v<UnarySubMap>);
static_assert(std::is_trivially_copyable_v<DerivedSubMap>);
static_assert(std::is_trivially_copyable_v<NestedSubMap>);
static_assert(slam::TriviallyCopyableRepresentation<UnarySubMap>);
static_assert(slam::TriviallyCopyableRepresentation<const UnarySubMap>);
static_assert(slam::TriviallyCopyableRepresentation<DerivedSubMap>);
static_assert(slam::TriviallyCopyableRepresentation<NestedSubMap>);

// Compatibility trait spellings remain exact Boolean wrappers around the concepts.
static_assert(slam::is_set_like_v<Range> == slam::SetLike<Range>);
static_assert(slam::is_bivariate_set_like_v<Product> == slam::BivariateSetLike<Product>);
static_assert(slam::is_relation_like_v<VariableRelation> == slam::RelationLike<VariableRelation>);
static_assert(slam::is_map_like_v<BinaryMap> == slam::MapLike<BinaryMap>);
static_assert(slam::is_map_over_v<BinaryMap, Product> == slam::MapOver<BinaryMap, Product>);
static_assert(slam::is_size_policy_v<EmptySize> == slam::SizePolicy<EmptySize>);
static_assert(slam::is_stride_policy_v<MatrixStride> == slam::StridePolicy<MatrixStride>);
static_assert(slam::is_position_like_v<StrongPosition> == slam::PositionLike<StrongPosition>);
static_assert(!slam::is_set_like_v<int>);
static_assert(!slam::is_map_over_v<int, ConcreteRange>);

using LegacyDistinctRelation = slam::Relation<std::int32_t, double>;
static_assert(std::same_as<typename LegacyDistinctRelation::RelationVec::value_type, std::int32_t>);

// This external set has neither iterators nor an isValid() method.
struct ExternalSet
{
  using PositionType = Position;
  using ElementType = int;

  PositionType size() const { return 2; }
  bool empty() const { return false; }
  int at(PositionType pos) const { return pos == 0 ? 10 : 40; }
};

struct NonAdvancingSet : ExternalSet
{
  struct Iterator
  {
    int operator*() const;
    bool operator!=(Iterator) const;
  };
  Iterator begin() const;
  Iterator end() const;
};

struct ExtraSetAliases : ExternalSet
{
  using FirstSetType = void;
  using SecondSetType = void;
};

struct NonComparingSet : ExternalSet
{
  struct Iterator
  {
    int operator*() const;
    Iterator& operator++();
  };
  Iterator begin() const;
  Iterator end() const;
};

static_assert(slam::SetLike<ExternalSet>);
static_assert(slam::SetLike<ExtraSetAliases>);
static_assert(!slam::BivariateSetLike<ExtraSetAliases>);
static_assert(!slam::IterableSetLike<ExternalSet>);
static_assert(!slam::Validatable<ExternalSet>);
static_assert(!slam::IterableSetLike<NonAdvancingSet>);
static_assert(!slam::IterableSetLike<NonComparingSet>);
static_assert(!slam::IterableSetLike<MinimalBivariateSet>);

// An algorithm over the public core needs only row traversal and coordinate access.
template <slam::BivariateSetLike S>
void checkRowCoordinates(const S& set)
{
  typename S::PositionType flat = 0;
  for(typename S::FirstSetType::PositionType first = 0; first < set.getFirstSet()->size(); ++first)
  {
    const auto row = set.getElements(first);
    decltype(row.size()) count = 0;
    for(const auto second : row)
    {
      ASSERT_LT(flat, set.size());
      const auto coordinate = set.at(flat++);
      EXPECT_EQ(first, coordinate.first);
      EXPECT_EQ(second, coordinate.second);
      EXPECT_GE(second, 0);
      EXPECT_LT(second, set.getSecondSet()->size());
      ++count;
    }
    EXPECT_EQ(count, row.size());
  }
  EXPECT_EQ(flat, set.size());
}

template <slam::RelationLike R>
auto relationPositions(const R& relation)
{
  using Pair = std::pair<typename R::FromSetType::PositionType, typename R::ToSetType::PositionType>;
  std::vector<Pair> pairs;
  for(typename R::FromSetType::PositionType from = 0; from < relation.fromSet()->size(); ++from)
  {
    for(const auto to : relation[from])
    {
      pairs.emplace_back(from, to);
    }
  }
  return pairs;
}

// These tests use maps with no invalid entries within their extent.
// Writable access is constrained separately from MapLike.
template <slam::MapLike M>
  requires requires(M& map) { map.flatValue(0, 0) = 0.; }
void fillComponents(M& map)
{
  using P = typename M::PositionType;
  for(P pos = 0; pos < map.size(); ++pos)
  {
    for(P component = 0; component < map.numComp(); ++component)
    {
      map.flatValue(pos, component) = 100. * pos + component;
    }
  }
}

TEST(slam_concepts, external_models_use_only_core_operations)
{
  MinimalBivariateSet set;
  checkRowCoordinates(set);
  EXPECT_EQ(0, set.at(0).first);
  EXPECT_EQ(1, set.at(2).first);
  EXPECT_EQ(2, set.at(2).second);

  MinimalRelation relation;
  using Pair = std::pair<std::int32_t, std::int64_t>;
  const std::vector<Pair> expected {{0, 2}, {1, 0}, {1, 1}};
  EXPECT_EQ(expected, relationPositions(relation));

  MinimalBivariateMap map;
  fillComponents(map);
  const auto& const_map = map;
  EXPECT_DOUBLE_EQ(102., const_map.flatValue(1, 2));
  EXPECT_EQ(1, const_map.index(1).first);
  EXPECT_EQ(2, const_map.index(1).second);
  EXPECT_EQ(&map.values[5], &const_map.flatValue(1, 2));
}

TEST(slam_concepts, sets_of_coordinates_can_be_used_as_sets)
{
  ExternalSet external;
  slam::Map<double, ExternalSet> map(&external);
  static_assert(slam::MapOver<decltype(map), ExternalSet>);
  fillComponents(map);
  EXPECT_TRUE(map.isValid());
  EXPECT_EQ(40, map.index(1));
  EXPECT_DOUBLE_EQ(100., map.flatValue(1, 0));

  using ExternalProduct = typename slam::ProductSet<ExternalSet, ExternalSet>::ConcreteSet;
  ExternalProduct product(&external, &external);
  ASSERT_TRUE(product.isValid());
  checkRowCoordinates(product);
  slam::ProductSet<ExternalSet, ExternalSet> virtual_product(&external, &external);
  EXPECT_TRUE(virtual_product.isValid());
  checkRowCoordinates(virtual_product);
  slam::Map<double, ExternalProduct> coordinate_map(&product);
  static_assert(slam::MapOver<decltype(coordinate_map), ExternalProduct>);
  fillComponents(coordinate_map);
  EXPECT_TRUE(coordinate_map.isValid());
  EXPECT_EQ(std::make_pair(Position {1}, Position {0}), coordinate_map.index(2));
  EXPECT_DOUBLE_EQ(200., coordinate_map.flatValue(2, 0));

  using NestedProduct = typename slam::ProductSet<ExternalProduct, ExternalSet>::ConcreteSet;
  NestedProduct nested(&product, &external);
  EXPECT_TRUE(nested.isValid());
  checkRowCoordinates(nested);
  EXPECT_EQ(8, nested.size());

  std::vector<Position> begins {0, 1, 1, 3, 4};
  std::vector<Position> indices {1, 0, 1, 0};
  auto relation = slam::make_variable_relation(&product, &external, begins, indices);
  static_assert(slam::RelationLike<decltype(relation)>);
  EXPECT_TRUE(relation.isValid());
  const std::vector<std::pair<Position, Position>> expected {{0, 1}, {2, 0}, {2, 1}, {3, 0}};
  EXPECT_EQ(expected, relationPositions(relation));
  EXPECT_EQ(std::make_pair(Position {1}, Position {0}), relation.fromSet()->at(2));
  EXPECT_EQ(40, relation.toSet()->at(relation[2][1]));

  std::vector<Position> reverse_indices {3, 1};
  auto reverse = slam::make_constant_relation(&external, &product, 1, reverse_indices);
  static_assert(slam::RelationLike<decltype(reverse)>);
  EXPECT_TRUE(reverse.isValid());
  EXPECT_EQ(std::make_pair(Position {1}, Position {1}), reverse.toSet()->at(reverse[0][0]));
}

TEST(slam_concepts, uniform_component_access_preserves_shape_and_selection)
{
  using ShapePolicy = policies::MultiDimStride<Position, 2>;
  using TensorMap = slam::Map<double, ConcreteRange, ArrayIndirection, ShapePolicy>;
  using TensorBivariateMap = slam::BivariateMap<double, Product, ArrayIndirection, ShapePolicy>;
  ShapePolicy::ShapeType shape {{2, 3}};
  ConcreteRange set(10, 12);
  TensorMap map(&set, 0., shape);
  static_assert(slam::MapLike<TensorMap>);
  fillComponents(map);
  const auto& const_map = map;
  EXPECT_EQ(6, map.numComp());
  EXPECT_EQ(11, map.index(1));
  for(Position i = 0; i < 2; ++i)
  {
    for(Position j = 0; j < 3; ++j)
    {
      EXPECT_DOUBLE_EQ(100. + 3 * i + j, map.value(1, i, j));
      EXPECT_EQ(&const_map.value(1, i, j), &const_map.flatValue(1, 3 * i + j));
    }
  }

  ConcreteRange selected(1, 2);
  slam::SubMap<TensorMap, ConcreteRange> submap(&map, selected);
  slam::SubMap<decltype(submap), ConcreteRange> nested(&submap, ConcreteRange(1));
  static_assert(slam::MapLike<decltype(nested)>);
  static_assert(!slam::MapOver<decltype(nested), ConcreteRange>);
  EXPECT_EQ(1, submap.set()->at(0));
  EXPECT_EQ(11, nested.index(0));
  const auto& const_submap = nested;
  const_submap.flatValue(0, 5) = 42.;
  EXPECT_DOUBLE_EQ(42., map.value(1, 1, 2));

  Product product(&set, &set);
  TensorBivariateMap bimap(&product, 0., shape);
  static_assert(slam::MapOver<TensorBivariateMap, Product>);
  fillComponents(bimap);
  EXPECT_EQ(std::make_pair(Position {1}, Position {0}), bimap.index(2));
  auto it = bimap.begin();
  for(Position pos = 0; pos < bimap.size(); ++pos)
  {
    for(Position component = 0; component < bimap.numComp(); ++component, ++it)
    {
      EXPECT_DOUBLE_EQ(100. * pos + component, *it);
      EXPECT_EQ(&bimap.flatValue(pos, component), &*it);
    }
  }
  EXPECT_EQ(it, bimap.end());
  auto row = bimap(1);
  EXPECT_EQ(std::make_pair(Position {1}, Position {1}), row.index(1));
  EXPECT_EQ(&bimap.flatValue(3, 5), &row.flatValue(1, 5));
}

TEST(slam_concepts, scalar_dynamic_and_view_maps)
{
  DynamicSet set(2);
  DynamicMap dynamic(&set);
  fillComponents(dynamic);
  const auto& const_dynamic = dynamic;
  EXPECT_EQ(1, dynamic.numComp());
  EXPECT_EQ(set.at(1), dynamic.index(1));
  EXPECT_DOUBLE_EQ(100., const_dynamic.flatValue(1, 0));
  EXPECT_EQ(&dynamic[1], &const_dynamic.flatValue(1, 0));

  ConcreteRange range(2);
  double data[2] {};
  UnaryMap view(&range, axom::ArrayView<double>(data, 2));
  const auto& const_view = view;
  const_view.flatValue(1, 0) = 9.;
  EXPECT_DOUBLE_EQ(9., data[1]);
  using ReadOnlyMap = slam::Map<const double, ConcreteRange, ConstViewIndirection>;
  ReadOnlyMap readonly(&range, axom::ArrayView<const double>(data, 2));
  static_assert(slam::MapLike<ReadOnlyMap>);
  static_assert(std::same_as<decltype(readonly.flatValue(0, 0)), const double&>);
  EXPECT_EQ(&data[1], &readonly.flatValue(1, 0));
}
}  // namespace slam_concept_test
