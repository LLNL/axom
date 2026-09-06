// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file slam_policy_contracts.cpp
 * \brief External policies exercised through SLAM owners and adapters.
 */

#include "gtest/gtest.h"

#include "axom/slam/BivariateMap.hpp"
#include "axom/slam/DynamicSet.hpp"
#include "axom/slam/ProductSet.hpp"
#include "axom/slam/RelationSet.hpp"

#include <type_traits>
#include <vector>

namespace
{
namespace slam = axom::slam;
namespace policies = slam::policies;
using Position = slam::DefaultPositionType;
using Range = typename slam::RangeSet<Position, Position>::ConcreteSet;

// These policies provide named operations, without ValuePolicy, tags, or IntType.
struct Size
{
  static constexpr Position DEFAULT_VALUE = 0;
  explicit Size(Position n = DEFAULT_VALUE) : count(n) { }
  Position size() const { return count; }
  Position& size() { return count; }
  bool empty() const { return count == 0; }
  bool isValid(bool) const { return count >= 0; }
  Position count;
};

struct Offset
{
  static constexpr Position DEFAULT_VALUE = 0;
  explicit Offset(Position n = DEFAULT_VALUE) : start(n) { }
  Position offset() const { return start; }
  bool isValid(bool) const { return start >= 0; }
  Position start;
};

// OrderedSet needs a scalar stride, without map shape metadata.
struct Stride
{
  static constexpr Position DEFAULT_VALUE = 1;
  explicit Stride(Position n = DEFAULT_VALUE) : step(n) { }
  Position stride() const { return step; }
  bool isValid(bool) const { return step != 0; }
  Position step;
};

struct Indirection
{
  using IndirectionPtrType = double*;
  using IndirectionResult = double&;
  using ConstIndirectionResult = const double&;
  Indirection() = default;
  explicit Indirection(double* data) : values(data) { }
  double& indirection(Position pos) { return values[pos]; }
  const double& indirection(Position pos) const { return values[pos]; }
  bool isValid(Position size, Position offset, Position stride, bool) const
  {
    return size == 0 || (values != nullptr && offset >= 0 && stride > 0);
  }
  double* values {nullptr};
};

struct Subset
{
  using ParentSetType = Range;
  Subset() = default;
  explicit Subset(ParentSetType* parent) : parent(parent) { }
  const ParentSetType* parentSet() const { return parent; }
  bool isSubset() const { return parent != nullptr; }

  template <typename Iterator>
  bool isValid(Iterator begin, Iterator end, bool) const
  {
    for(; begin != end; ++begin)
    {
      if(parent == nullptr || *begin < 0 || *begin >= parent->size())
      {
        return false;
      }
    }
    return true;
  }
  ParentSetType* parent {nullptr};
};

struct ImmutableSize : Size
{
  using Size::Size;
  Position size() const { return count; }
};
struct FinalSize final : Size
{
  using Size::Size;
};
struct MissingEmpty : Size
{
  void empty() const;
};
struct MissingValidation : Indirection
{
  void isValid(Position, Position, Position, bool) const;
};
struct MissingSubsetParent : Subset
{
  void parentSet() const;
};

template <typename S>
concept CanFormOrderedSet = requires { sizeof(slam::OrderedSet<Position, Position, S>); };
template <typename S>
concept CanFormDynamicSet = requires { sizeof(slam::DynamicSet<Position, Position, S>); };
template <typename I>
concept CanFormIndirectSet =
  requires { sizeof(slam::OrderedSet<Position, double, Size, Offset, Stride, I>); };

static_assert(slam::SizePolicy<const Size&>);
static_assert(slam::StridePolicy<const Stride&>);
static_assert(slam::OffsetPolicy<const Offset&>);
static_assert(slam::SubsetPolicy<const Subset&>);
static_assert(CanFormOrderedSet<ImmutableSize>);
static_assert(!CanFormDynamicSet<ImmutableSize>);
static_assert(!CanFormDynamicSet<const Size>);
static_assert(!CanFormOrderedSet<const Size>);
static_assert(!CanFormOrderedSet<Size&>);
static_assert(!CanFormOrderedSet<FinalSize>);
static_assert(!CanFormOrderedSet<MissingEmpty>);
static_assert(!CanFormIndirectSet<MissingValidation>);
static_assert(!slam::detail::OrderedSetSubsetPolicy<MissingSubsetParent>);
static_assert(!slam::detail::MapStridePolicyFor<Stride, Position>);

TEST(slam_policy_contracts, ordered_and_dynamic_sets_use_named_operations)
{
  using Set = slam::OrderedSet<Position, double, Size, Offset, Stride, Indirection>;
  double values[] {1., 2., 3., 4.};
  Set set(Set::SetBuilder().size(2).offset(1).stride(2).data(values));
  ASSERT_TRUE(set.isValid());
  EXPECT_EQ(2, set.size());
  EXPECT_FALSE(set.empty());
  EXPECT_DOUBLE_EQ(2., set.at(0));
  EXPECT_DOUBLE_EQ(4., set.at(1));
  auto it = set.begin();
  EXPECT_DOUBLE_EQ(2., *it++);
  EXPECT_DOUBLE_EQ(4., *it++);
  EXPECT_EQ(set.end(), it);

  slam::DynamicSet<Position, Position, Size> dynamic(2);
  EXPECT_EQ(2, dynamic.insert(8));
  EXPECT_EQ(3, dynamic.size());
  EXPECT_EQ(8, dynamic.at(2));
  dynamic.reset(1);
  EXPECT_TRUE(dynamic.isValid());
  EXPECT_EQ(1, dynamic.size());

  using Selected =
    slam::OrderedSet<Position, Position, Size, Offset, Stride, policies::NoIndirection<Position, Position>, Subset>;
  Range parent(4);
  Selected selected(Selected::SetBuilder().size(2).offset(1).parent(&parent));
  EXPECT_TRUE(selected.isSubset());
  EXPECT_EQ(&parent, selected.parentSet());
  EXPECT_TRUE(selected.isValid());
  Selected invalid(Selected::SetBuilder().size(5).parent(&parent));
  EXPECT_FALSE(invalid.isValid());
}

// Map constructs a stride from a shape, without requiring default construction.
struct MapStride
{
  using IndexType = Position;
  using ShapeType = Position;
  static constexpr int NumDims = 1;
  static ShapeType DefaultSize() { return 1; }
  explicit MapStride(ShapeType shape) : components(shape) { }
  Position stride() const { return components; }
  ShapeType shape() const { return components; }
  Position components;
};

struct ScalarStrideWithUnrelatedStrides : MapStride
{
  using MapStride::MapStride;
  void strides() const { }
};

// This shape has subscripting, but no iterators or aggregate initialization.
struct Shape
{
  Shape() = default;
  Shape(Position rows, Position columns) : entries {rows, columns} { }
  Position operator[](int dim) const { return entries[dim]; }
  Position entries[2] {1, 1};
};

struct TensorStride
{
  using IndexType = Position;
  using ShapeType = Shape;
  static constexpr int NumDims = 2;
  static ShapeType DefaultSize() { return {}; }
  explicit TensorStride(Shape shape) : dimensions(shape) { }
  Position stride() const { return dimensions[0] * dimensions[1]; }
  Shape shape() const { return dimensions; }
  Shape strides() const { return {dimensions[1], 1}; }
  Shape dimensions;
};

struct UnindexableShape
{ };
struct BadTensorStride : TensorStride
{
  using ShapeType = UnindexableShape;
  explicit BadTensorStride(ShapeType);
  static ShapeType DefaultSize();
  ShapeType shape() const;
  ShapeType strides() const;
};
struct FinalStride final : MapStride
{
  using MapStride::MapStride;
};

// Storage is accessed through the descriptor, never through Buffer::operator[].
struct Buffer
{
  Buffer() = delete;
  explicit Buffer(Position count, double value = 0.) : values(count, value) { }
  std::size_t size() const { return values.size(); }
  bool empty() const { return values.empty(); }
  void resize(Position count) { values.resize(count); }
  std::vector<double> values;
};

struct Storage final
{
  using IndirectionBufferType = Buffer;
  using IndirectionResult = double&;
  using ConstIndirectionResult = const double&;
  using ResultPtr = double*;
  using ConstResultPtr = const double*;
  static constexpr bool IsMutableBuffer = true;
  static ResultPtr getIndirection(Buffer& buffer, Position pos = 0)
  {
    return buffer.values.data() + pos;
  }
  static ConstResultPtr getConstIndirection(const Buffer& buffer, Position pos = 0)
  {
    return buffer.values.data() + pos;
  }
  static Buffer create(Position size, const double& value, int) { return Buffer(size, value); }
};

template <typename S>
concept CanFormMapWithStride = requires { sizeof(slam::Map<double, Range, Storage, S>); };
template <typename I>
concept CanFormMapWithStorage = requires { sizeof(slam::Map<double, Range, I, MapStride>); };
template <typename M>
concept CanFill = requires(M& map) { map.fill(1.); };
template <typename M>
concept CanCopy = requires(M& map, const M& source) { map.copy(source); };
template <typename M>
concept HasPerDimensionStrides = requires(const M& map) { map.strides(); };

static_assert(CanFormMapWithStorage<Storage>);
static_assert(!CanFormMapWithStorage<const Storage>);
static_assert(!CanFormMapWithStorage<Storage&>);
static_assert(CanFormMapWithStride<MapStride>);
static_assert(CanFormMapWithStride<TensorStride>);
static_assert(!CanFormMapWithStride<const MapStride>);
static_assert(!CanFormMapWithStride<MapStride&>);
static_assert(!CanFormMapWithStride<FinalStride>);
static_assert(!CanFormMapWithStride<BadTensorStride>);

TEST(slam_policy_contracts, map_uses_static_storage_and_shaped_stride)
{
  using Map = slam::Map<double, Range, Storage, MapStride>;
  Range range(2);
  Map map(&range, 4., 3);
  EXPECT_TRUE(map.isValid());
  EXPECT_EQ(3, map.numComp());
  EXPECT_DOUBLE_EQ(4., map.flatValue(1, 2));
  map.fill(8.);
  EXPECT_DOUBLE_EQ(8., *map.begin());
  EXPECT_DOUBLE_EQ(8., map.set_begin().value(2));

  Map supplied(&range, Buffer(1), 3);
  EXPECT_TRUE(supplied.isValid());
  EXPECT_EQ(6u, supplied.data().size());
  supplied.copy(map);
  EXPECT_DOUBLE_EQ(8., supplied.flatValue(1, 2));
  supplied.clear();
  EXPECT_DOUBLE_EQ(0., supplied.flatValue(1, 2));

  double values[] {1., 2., 3., 4., 5., 6.};
  Map built(Map::MapBuilder().set(&range).stride(3).data(values));
  EXPECT_TRUE(built.isValid());
  EXPECT_DOUBLE_EQ(6., built.flatValue(1, 2));

  using Product = typename slam::ProductSet<Range, Range>::ConcreteSet;
  Product product(&range, &range);
  slam::BivariateMap<double, Product, Storage, ScalarStrideWithUnrelatedStrides> bivariate(&product,
                                                                                           6.,
                                                                                           3);
  static_assert(!HasPerDimensionStrides<decltype(bivariate)>);
  EXPECT_TRUE(bivariate.isValid());
  EXPECT_DOUBLE_EQ(6., bivariate(1).flatValue(1, 2));
  EXPECT_DOUBLE_EQ(6., bivariate.set_begin().value(2));
  EXPECT_EQ(&bivariate.flatValue(3, 2), bivariate.findValue(1, 1, 2));

  using TensorMap = slam::Map<double, Range, Storage, TensorStride>;
  TensorMap tensor(&range, 0., Shape(2, 3));
  tensor.value(1, 1, 2) = 9.;
  EXPECT_DOUBLE_EQ(9., tensor.flatValue(1, 5));
  EXPECT_DOUBLE_EQ(9., (tensor.set_begin() + 1).value(1, 2));
  slam::SubMap<TensorMap, Range> submap(&tensor, Range(1, 2));
  EXPECT_DOUBLE_EQ(9., submap.set_begin().value(5));
  EXPECT_DOUBLE_EQ(9., submap.set_begin().value(1, 2));

  using ReadOnly = slam::Map<double, Range, policies::ArrayViewIndirection<Position, const double>>;
  static_assert(slam::MapLike<ReadOnly>);
  static_assert(!CanFill<ReadOnly>);
  static_assert(!CanCopy<ReadOnly>);
  ReadOnly readOnly(&range, axom::ArrayView<const double>(values, 2));
  EXPECT_DOUBLE_EQ(2., readOnly.flatValue(1, 0));
}

struct ScalarParent
{
  using PositionType = Position;
  Position size() const { return 2; }
  Position numComp() const { return 2; }
  Position shape() const { return 2; }
  Position index(Position pos) const { return 10 + pos; }
  double& operator[](Position pos) { return values[pos]; }
  double values[4] {1., 2., 3., 4.};
};

struct SelectedPositions
{
  using PositionType = Position;
  using ElementType = Position;
  Position size() const { return 1; }
  bool empty() const { return false; }
  Position at(Position) const { return 1; }
};

struct TemporaryShapedValue : ScalarParent
{
  const double& operator[](Position pos) { return values[pos]; }
  double flatValue(Position, Position, Position) const { return 0.; }
};

struct BadRangeParent : ScalarParent
{
  struct Iterator
  {
    Iterator() = default;
    Iterator(BadRangeParent*, Position);
    double operator*() const;
    Position flatIndex() const;
  };
  Iterator set_begin();
};

// No set_end(), iterator value(), numComp(), or operator->() is needed.
struct RangeParent : ScalarParent
{
  struct Iterator
  {
    Iterator() = default;
    Iterator(RangeParent* parent, Position pos)
      : position(pos)
      , range(parent->values + pos * 2, pos < parent->size() ? 2 : 0)
    { }
    const axom::ArrayView<double>& operator*() const { return range; }
    Position flatIndex() const { return position; }
    Position position {};
    axom::ArrayView<double> range;
  };
  Iterator set_begin() { return {this, 0}; }
};

template <typename M>
concept HasRangeTraversal = requires(M& map) {
  map.set_begin();
  map.set_end();
};
template <typename M>
concept HasShapedAccess = requires(M& map) { map.flatValue(0, 1, 1); };

static_assert(slam::detail::SubMapSource<ScalarParent>);
static_assert(!slam::detail::SubMapRangeSource<ScalarParent>);
static_assert(!slam::detail::SubMapRangeSource<BadRangeParent>);
static_assert(slam::detail::SubMapRangeSource<RangeParent>);

TEST(slam_policy_contracts, scalar_submap_does_not_require_parent_range_iteration)
{
  ScalarParent parent;
  slam::SubMap<ScalarParent, SelectedPositions> submap(&parent, SelectedPositions {});
  static_assert(slam::MapLike<decltype(submap)>);
  static_assert(!HasRangeTraversal<decltype(submap)>);
  static_assert(!HasShapedAccess<decltype(submap)>);
  using TemporaryShapedSubMap = slam::SubMap<TemporaryShapedValue, SelectedPositions>;
  static_assert(!HasShapedAccess<TemporaryShapedSubMap>);
  EXPECT_TRUE(submap.isValid());
  EXPECT_EQ(11, submap.index(0));
  EXPECT_DOUBLE_EQ(4., submap.flatValue(0, 1));
  EXPECT_DOUBLE_EQ(3., *submap.begin());
  EXPECT_EQ(2, submap.end() - submap.begin());

  RangeParent rangeParent;
  slam::SubMap<RangeParent, Range> ranges(&rangeParent, Range(1, 2));
  auto it = ranges.set_begin();
  EXPECT_DOUBLE_EQ(4., it.value(1));
  EXPECT_DOUBLE_EQ(3., (*it)[0]);
  EXPECT_EQ(2, it->size());
  EXPECT_EQ(1, it.flatIndex());
  EXPECT_EQ(2, it.numComp());
  EXPECT_EQ(ranges.set_end(), ++it);
}
}  // namespace
