// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 * \file slam_map_SubMap.cpp
 *
 * \brief Unit tests for Slam's SubMap
 */

#include <iterator>
#include "gtest/gtest.h"

#include "axom/slic.hpp"

#include "axom/slam/Utilities.hpp"
#include "axom/slam/RangeSet.hpp"
#include "axom/slam/SubMap.hpp"
#include "axom/slam/ProductSet.hpp"
#include "axom/slam/BivariateMap.hpp"
#include "axom/slam/RelationBuilders.hpp"
#include "axom/slam/RelationSet.hpp"

#include <type_traits>
#include <utility>

namespace
{
namespace slam = axom::slam;

using SetBase = slam::Set<>;
using PositionType = SetBase::PositionType;
using ElementType = SetBase::ElementType;

using RangeSetType = slam::RangeSet<PositionType, ElementType>;

template <typename T>
using SuperMap = slam::Map<T, SetBase>;

using OrderedSetType =
  axom::slam::OrderedSet<PositionType,
                         ElementType,
                         slam::policies::RuntimeSize<PositionType>,
                         slam::policies::ZeroOffset<PositionType>,
                         slam::policies::StrideOne<PositionType>,
                         slam::policies::STLVectorIndirection<PositionType, ElementType>>;

constexpr double multFac = 1.0001;

PositionType const MAX_SET_SIZE = 10;

template <typename T>
T getValue(int idx)
{
  return static_cast<T>(idx * multFac);
}

//A struct to construct the Map for testing
template <typename T>
struct MapForTest
{
  MapForTest(int size)
    : set_data(size)
    , s(OrderedSetType::SetBuilder().size(size).data(&set_data))
    , m(&s)
  {
    SLIC_INFO("Initializing set of size " << s.size() << " and map on the set ");

    for(auto i : s.positions())
    {
      s[i] = 100 + i;
      m[i] = getValue<T>(i);
    }

    EXPECT_TRUE(s.isValid());
    EXPECT_TRUE(m.isValid());
  }

  std::vector<PositionType> set_data;
  OrderedSetType s;
  SuperMap<T> m;
};

}  // namespace

TEST(slam_map, construct_empty_subsetmap)
{
  slam::SubMap<SuperMap<int>, RangeSetType> m;
  EXPECT_TRUE(m.isValid(true));
  EXPECT_EQ(m.size(), 0);
  EXPECT_EQ(m.begin(), m.end());
  EXPECT_EQ(m.set_begin(), m.set_end());

  slam::SubMap<decltype(m), RangeSetType> nested(&m, RangeSetType {});
  EXPECT_TRUE(nested.isValid());
  EXPECT_EQ(nested.begin(), nested.end());
  EXPECT_EQ(nested.set_begin(), nested.set_end());
}

template <typename T>
bool constructAndTestSubMap()
{
  using MapType = SuperMap<T>;

  MapForTest<T> mft(MAX_SET_SIZE);
  MapType& m = mft.m;
  OrderedSetType& s = mft.s;

  int submapOffset = 3;
  int submapSize = 5;
  {
    using SubMapType = slam::SubMap<MapType, RangeSetType>;

    SLIC_INFO("Creating the Subset.");
    RangeSetType ss(submapOffset, submapOffset + submapSize);
    SubMapType ssm(&m, ss);
    EXPECT_TRUE(ssm.isValid(true));

    auto flatIter = ssm.begin();
    EXPECT_EQ(flatIter.operator->(), &ssm[0]);

    SLIC_INFO("Checking the elements.");
    for(auto idx = 0; idx < ssm.size(); ++idx)
    {
      auto expVal = getValue<T>(submapOffset + idx);
      EXPECT_EQ(expVal, ssm[idx]);
    }

    SLIC_INFO("Checking the elements using SubMap iterator.");
    int cnt = 0;
    for(auto it = ssm.begin(); it != ssm.end(); ++it, ++cnt)
    {
      // Check iterator's .index() function
      {
        auto subMapElt = it.index();

        auto setIdx = ss[cnt];
        auto expSetElt = s[setIdx];
        EXPECT_EQ(expSetElt, subMapElt);
      }

      // Check iterator's value access functions
      {
        auto expectedValue = getValue<T>(submapOffset + cnt);
        EXPECT_EQ(expectedValue, *it);
        EXPECT_EQ(expectedValue, it[0]);

        auto expVal = getValue<T>(submapOffset + cnt);
        EXPECT_EQ(expVal, expectedValue);
      }
    }
  }

  {
    using SubMapType = slam::SubMap<MapType, OrderedSetType>;

    SLIC_INFO("Creating Subset 2");
    std::vector<PositionType> subset_indices_data(submapSize);
    for(int i = 0; i < submapSize; i++)
    {
      subset_indices_data[i] = i * 2;
    }
    OrderedSetType subset_indices = OrderedSetType::SetBuilder().size(5).data(&subset_indices_data);

    SubMapType ssm(&m, subset_indices);

    SLIC_INFO("Checking the elements.");
    for(PositionType idx = 0; idx < submapSize; ++idx)
    {
      EXPECT_EQ(ssm[idx], getValue<T>(subset_indices[idx]));
      EXPECT_EQ(ssm.value(idx), getValue<T>(subset_indices[idx]));
      EXPECT_EQ(ssm.index(idx), s[subset_indices[idx]]);
    }

    SLIC_INFO("Checking the elements using SubMap range iterator.");
    {
      int cnt = 0;
      for(auto it = ssm.set_begin(); it != ssm.set_end(); ++it, ++cnt)
      {
        auto expectedValue = getValue<T>(subset_indices[cnt]);
        EXPECT_EQ(expectedValue, (*it)[0]);
        EXPECT_EQ(expectedValue, it(0));
        EXPECT_EQ(expectedValue, it.value(0));
        EXPECT_EQ(it.index(), s[subset_indices[cnt]]);
      }
      EXPECT_EQ(cnt, subset_indices.size());
    }
  }
  return true;
}

TEST(slam_map, construct_int_submap) { EXPECT_TRUE(constructAndTestSubMap<int>()); }

TEST(slam_map, construct_double_submap) { EXPECT_TRUE(constructAndTestSubMap<double>()); }

template <typename T>
bool constructBySubMap()
{
  //This tests modifying the values in the original map via a Submap
  //Create a Map, then create a SubMap on the Map, negate all values covered by
  //the Submap, then check the values are negative in the SubMap range.
  using MapType = SuperMap<T>;
  using SubMapType = slam::SubMap<MapType, RangeSetType>;

  MapForTest<T> mft(MAX_SET_SIZE);
  MapType& m = mft.m;

  int submapOffset = 3;
  int submapSize = 5;

  SLIC_INFO("Creating the Subset.");
  RangeSetType ss(submapOffset, submapOffset + submapSize);
  SubMapType ssm(&m, ss);
  EXPECT_TRUE(m.isValid(true));

  SLIC_INFO("Negating elements");
  for(PositionType idx = 0; idx < submapSize; ++idx)
  {
    ssm[idx] = -getValue<T>(submapOffset + idx);
  }

  SLIC_INFO("Checking the elements.");
  for(PositionType idx = 0; idx < m.size(); ++idx)
  {
    T val = getValue<T>(idx);
    if(idx >= submapOffset && idx < submapOffset + submapSize)
    {
      val = -getValue<T>(idx);
    }
    EXPECT_EQ(m[idx], val);
  }

  SLIC_INFO("Checking elements");

  return true;
}

TEST(slam_map, submap_of_submap)
{
  // A SubMap is itself a map, so we can take a SubMap from it
  using ProductSetType = typename slam::ProductSet<RangeSetType, RangeSetType>::ConcreteSet;
  using BMapType = slam::BivariateMap<double, ProductSetType>;
  using RowSubMap = typename BMapType::SubMapType;
  using NestedSubMap = slam::SubMap<RowSubMap, typename RowSubMap::IndexSetType>;

  static_assert(slam::detail::SubMapSource<BMapType>);
  static_assert(slam::detail::SubMapSource<RowSubMap>, "a SubMap can serve as a super-map");
  static_assert(slam::detail::SubMapSource<NestedSubMap>, "and composition does not bottom out");

  constexpr PositionType NROWS = 4, NCOLS = 5;
  RangeSetType rows(NROWS), cols(NCOLS);
  ProductSetType prod(&rows, &cols);
  BMapType bmap(prod, 0.0);

  // each value encodes its own coordinate, so a misindex cannot alias
  for(PositionType i = 0; i < NROWS; ++i)
  {
    for(PositionType j = 0; j < NCOLS; ++j)
    {
      bmap(i, j) = 100.0 * i + j;
    }
  }

  RowSubMap row2 = bmap(2);
  ASSERT_EQ(NCOLS, row2.size());

  // columns 1..3 of row 2
  constexpr PositionType FIRST = 1, COUNT = 3;
  auto inner = typename RowSubMap::IndexSetType::SetBuilder().size(COUNT).offset(FIRST);
  NestedSubMap mid(&row2, inner);

  ASSERT_EQ(COUNT, mid.size());
  for(PositionType k = 0; k < mid.size(); ++k)
  {
    EXPECT_EQ(200.0 + (k + FIRST), mid[k]);
  }

  // element iteration
  double sum = 0.0;
  PositionType visited = 0;
  for(auto it = mid.begin(); it != mid.end(); ++it, ++visited)
  {
    sum += *it;
  }
  EXPECT_EQ(COUNT, visited);
  EXPECT_EQ(201.0 + 202.0 + 203.0, sum);

  // Range iteration. This is a regression test for an index-space bug in
  // SubMap::RangeIterator::advance(): stepping the parent iterator by a
  // difference computed from the parent's flatIndex() only works when the
  // parent is a Map or BivariateMap.
  double rangeSum = 0.0;
  PositionType rangeVisited = 0;
  for(auto it = mid.set_begin(); it != mid.set_end(); ++it, ++rangeVisited)
  {
    rangeSum += (*it)[0];
  }
  EXPECT_EQ(COUNT, rangeVisited);
  EXPECT_EQ(201.0 + 202.0 + 203.0, rangeSum);

  // The index set contains positions in the immediate parent. index() follows
  // both selections to the bivariate coordinate associated with each value.
  for(PositionType k = 0; k < mid.size(); ++k)
  {
    EXPECT_EQ(k + FIRST, mid.set()->at(k));
    const auto coordinate = mid.index(k);
    EXPECT_EQ(coordinate, row2.index(k + FIRST));
    EXPECT_EQ(2, coordinate.first);
    EXPECT_EQ(k + FIRST, coordinate.second);
  }
}

TEST(slam_map, construct_with_int_submap) { EXPECT_TRUE(constructBySubMap<int>()); }

TEST(slam_map, construct_with_double_submap) { EXPECT_TRUE(constructBySubMap<double>()); }

namespace
{
template <typename Parent>
void checkNestedSubMapReferences(Parent& parent)
{
  using First = slam::SubMap<Parent, RangeSetType>;
  using Second = slam::SubMap<const First, RangeSetType>;
  using Third = slam::SubMap<Second, RangeSetType>;
  using Reference = decltype(parent[PositionType {}]);

  static_assert(std::same_as<typename First::reference, Reference>);
  static_assert(std::same_as<typename Second::reference, Reference>);
  static_assert(std::same_as<typename Third::const_reference, Reference>);
  static_assert(std::same_as<decltype(std::declval<const Third&>()[0]), Reference>);
  static_assert(std::same_as<decltype(std::declval<Third&>()(0)), Reference>);
  static_assert(std::same_as<decltype(std::declval<Third&>().set_begin().value(0)), Reference>);
  static_assert(std::same_as<decltype((*std::declval<Third&>().set_begin())[0]), Reference>);
  static_assert(std::random_access_iterator<typename Third::iterator>);
  static_assert(std::bidirectional_iterator<typename Third::range_iterator>);

  First first(&parent, RangeSetType(1, 8));
  const First& shallowConst = first;
  Second second(&shallowConst, RangeSetType(2, 6));
  Third third(&second, RangeSetType(1, 3));
  ASSERT_TRUE(third.isValid());

  for(PositionType i = 0; i < third.size(); ++i)
  {
    const PositionType originalPosition = i + 4;
    EXPECT_EQ(&third[i], &parent[originalPosition]);
    EXPECT_EQ(&third(i), &parent[originalPosition]);
    EXPECT_EQ(third.index(i), parent.index(originalPosition));
    auto it = third.set_begin() + i;
    EXPECT_EQ(&it.value(0), &parent[originalPosition]);
    EXPECT_EQ(it.index(), third.index(i));
    EXPECT_EQ(it.flatIndex(), originalPosition);
    if constexpr(!std::is_const_v<std::remove_reference_t<Reference>>)
    {
      third.value(i) = 900.0 + i;
      EXPECT_EQ(parent[originalPosition], 900.0 + i);
    }
  }
}

// A parent exposing only the operations consumed by SubMap: no storage-policy,
// value-type, or iterator aliases are needed on the parent.
class PolicyFreeParent
{
public:
  using PositionType = SetBase::PositionType;

  explicit PolicyFreeParent(double* values) : m_values(values) { }
  PositionType index(PositionType pos) const { return m_set.at(pos); }
  PositionType size() const { return m_set.size(); }
  PositionType numComp() const { return 2; }
  PositionType shape() const { return numComp(); }
  double& operator[](PositionType pos) const { return m_values[pos]; }

  class Rows
  {
  public:
    Rows() = default;
    Rows(const PolicyFreeParent* parent, PositionType pos)
      : m_position(pos)
      , m_values(parent->m_values + pos * parent->numComp(), parent->numComp())
    { }
    const axom::ArrayView<double>& operator*() const { return m_values; }
    const axom::ArrayView<double>* operator->() const { return &m_values; }
    double& operator()(PositionType component) const { return m_values[component]; }
    double& value(PositionType component) const { return m_values[component]; }
    PositionType flatIndex() const { return m_position; }
    PositionType numComp() const { return m_values.size(); }

  private:
    PositionType m_position {0};
    axom::ArrayView<double> m_values;
  };
  Rows set_begin() const { return Rows(this, 0); }
  Rows set_end() const { return Rows(this, size()); }

private:
  RangeSetType m_set {100, 104};
  double* m_values;
};

template <typename T>
concept HasStoragePolicies = requires {
  typename T::StridePolicyType;
  typename T::IndirectionPolicy;
};

template <typename T>
concept HasWritableStride = requires(T& map) { map.stride() = 2; };
}  // namespace

TEST(slam_map, nested_submaps_preserve_mutable_and_deep_const_references)
{
  MapForTest<double> data(MAX_SET_SIZE);
  checkNestedSubMapReferences(data.m);
  checkNestedSubMapReferences(std::as_const(data.m));

  using Product = typename slam::ProductSet<RangeSetType, RangeSetType>::ConcreteSet;
  RangeSetType rows(2), columns(5);
  Product product(&rows, &columns);
  slam::BivariateMap<double, Product> bmap(&product);
  checkNestedSubMapReferences(bmap);
  checkNestedSubMapReferences(std::as_const(bmap));

  using ViewMap =
    slam::Map<double, RangeSetType, slam::policies::ArrayViewIndirection<PositionType, double>>;
  RangeSetType domain(MAX_SET_SIZE);
  double values[MAX_SET_SIZE] {};
  ViewMap view(&domain, axom::ArrayView<double>(values, MAX_SET_SIZE));
  checkNestedSubMapReferences(view);
  checkNestedSubMapReferences(std::as_const(view));
}

TEST(slam_map, submap_accepts_parent_operations_without_storage_policies)
{
  double values[] {0, 1, 10, 11, 20, 21, 30, 31};
  PolicyFreeParent parent(values);
  using Subset = slam::SubMap<const PolicyFreeParent, RangeSetType>;
  static_assert(slam::detail::SubMapSource<PolicyFreeParent>);
  static_assert(!HasStoragePolicies<PolicyFreeParent>);
  static_assert(!HasStoragePolicies<Subset>);
  static_assert(std::same_as<Subset::reference, double&>);
  static_assert(!HasWritableStride<Subset>);

  const Subset subset(&parent, RangeSetType(1, 3));
  ASSERT_TRUE(subset.isValid());
  EXPECT_EQ(subset.index(0), 101);
  EXPECT_EQ(subset.index(1), 102);
  EXPECT_EQ(subset.numComp(), 2);
  subset(0, 1) = 42;
  EXPECT_EQ(values[3], 42);
  auto it = subset.set_end();
  --it;
  EXPECT_EQ(it.value(1), 21);
  --it;
  EXPECT_EQ(it(1), 42);
  EXPECT_EQ(it, subset.set_begin());
}

TEST(slam_map, submap_forwards_parent_component_state)
{
  using MapType = slam::Map<double,
                            RangeSetType,
                            slam::policies::ArrayIndirection<PositionType, double>,
                            slam::policies::RuntimeStride<PositionType>>;
  RangeSetType domain(4);
  MapType parent(&domain, 1.0, 2);
  slam::SubMap<MapType, RangeSetType> subset(&parent, RangeSetType(1, 3));
  slam::SubMap<decltype(subset), RangeSetType> nested(&subset, RangeSetType(1));
  EXPECT_EQ(nested.numComp(), 2);

  parent = MapType(&domain, 7.0, 3);
  ASSERT_TRUE(parent.isValid());
  EXPECT_EQ(subset.numComp(), 3);
  EXPECT_EQ(nested.numComp(), 3);
  EXPECT_EQ(nested.shape(), 3);
  EXPECT_EQ(nested.stride(), 3);
  EXPECT_EQ(nested.end() - nested.begin(), 3);
  nested(0, 2) = 42;
  EXPECT_EQ(parent(1, 2), 42);
  EXPECT_EQ(nested.set_begin()->size(), 3);
}

TEST(slam_map, nested_submaps_forward_multidimensional_access)
{
  using Shape = axom::StackArray<PositionType, 2>;
  using MapType = slam::Map<double,
                            RangeSetType,
                            slam::policies::ArrayIndirection<PositionType, double>,
                            slam::policies::MultiDimStride<PositionType, 2>>;
  RangeSetType domain(4);
  MapType parent(&domain, 0.0, Shape {{2, 3}});
  slam::SubMap<MapType, RangeSetType> subset(&parent, RangeSetType(1, 4));
  slam::SubMap<decltype(subset), RangeSetType> nested(&subset, RangeSetType(1, 3));
  nested(0, 1, 2) = 42;
  EXPECT_EQ(parent(2, 1, 2), 42);
  EXPECT_EQ(nested.set_begin().value(1, 2), 42);

  parent = MapType(&domain, 7.0, Shape {{3, 2}});
  EXPECT_EQ(nested.shape(), (Shape {{3, 2}}));
  nested(0, 2, 1) = 81;
  EXPECT_EQ(parent(2, 2, 1), 81);
  EXPECT_EQ(nested.set_begin().value(2, 1), 81);

  slam::SubMap<const MapType, RangeSetType> constSubset(&std::as_const(parent), RangeSetType(1, 4));
  slam::SubMap<decltype(constSubset), RangeSetType> constNested(&constSubset, RangeSetType(1, 3));
  static_assert(std::same_as<decltype(constNested(0, 2, 1)), const double&>);
  EXPECT_EQ(constNested(0, 2, 1), 81);
  EXPECT_EQ(constNested.set_begin().value(2, 1), 81);
}

TEST(slam_map, empty_sparse_rows_and_nested_subsets_have_empty_ranges)
{
  RangeSetType from(3), to(4);
  std::vector<PositionType> begins {0, 0, 2, 2};
  std::vector<PositionType> indices {1, 3};
  auto relation = slam::make_variable_relation(from, to, begins, indices);
  slam::RelationSet<decltype(relation)> domain(&relation);
  slam::BivariateMap<double, decltype(domain)> parent(&domain);
  ASSERT_TRUE(parent.isValid());

  for(PositionType row : {PositionType {0}, PositionType {2}})
  {
    auto subset = parent(row);
    EXPECT_EQ(subset.size(), 0);
    EXPECT_EQ(subset.begin(), subset.end());
    EXPECT_EQ(subset.set_begin(), subset.set_end());
    EXPECT_EQ(parent.set_begin(row), parent.set_end(row));
    auto constSubset = std::as_const(parent)(row);
    EXPECT_EQ(constSubset.set_begin(), constSubset.set_end());
    slam::SubMap<decltype(subset), RangeSetType> nested(&subset, RangeSetType {});
    EXPECT_EQ(nested.set_begin(), nested.set_end());
  }

  auto row = parent(1);
  slam::SubMap<decltype(row), RangeSetType> emptyTail(&row, RangeSetType(2, 2));
  EXPECT_EQ(emptyTail.set_begin(), emptyTail.set_end());

  std::vector<PositionType> emptyBegins {0, 0, 0, 0};
  std::vector<PositionType> emptyIndices;
  auto emptyRelation = slam::make_variable_relation(from, to, emptyBegins, emptyIndices);
  slam::RelationSet<decltype(emptyRelation)> emptyDomain(&emptyRelation);
  slam::BivariateMap<double, decltype(emptyDomain)> emptyParent(&emptyDomain);
  ASSERT_TRUE(emptyParent.isValid());
  EXPECT_EQ(emptyParent.set_begin(), emptyParent.set_end());
  for(PositionType i = 0; i < from.size(); ++i)
  {
    auto emptyRow = emptyParent(i);
    EXPECT_EQ(emptyRow.set_begin(), emptyRow.set_end());
  }
}

TEST(slam_map, reordered_nested_subsets_support_backward_range_iteration)
{
  MapForTest<double> data(MAX_SET_SIZE);
  std::vector<PositionType> firstIndices {7, 2, 5};
  OrderedSetType firstSet = OrderedSetType::SetBuilder().size(3).data(&firstIndices);
  slam::SubMap<SuperMap<double>, OrderedSetType> first(&data.m, firstSet);
  std::vector<PositionType> secondIndices {2, 0, 1};
  OrderedSetType secondSet = OrderedSetType::SetBuilder().size(3).data(&secondIndices);
  slam::SubMap<decltype(first), OrderedSetType> second(&first, secondSet);
  const PositionType originalPositions[] {5, 7, 2};

  auto it = second.set_end();
  for(PositionType i = second.size(); i > 0;)
  {
    --i;
    --it;
    EXPECT_EQ(&it.value(0), &data.m[originalPositions[i]]);
    EXPECT_EQ(it.index(), data.m.set()->at(originalPositions[i]));
    EXPECT_EQ(it.flatIndex(), originalPositions[i]);
  }
  EXPECT_EQ(it, second.set_begin());
  it += 2;
  EXPECT_EQ(&it.value(0), &data.m[2]);
  it -= 1;
  EXPECT_EQ(&it.value(0), &data.m[7]);
}

TEST(slam_map, submap_copies_index_metadata_and_preserves_element_type)
{
  using Set = typename slam::RangeSet<PositionType, double>::ConcreteSet;
  Set domain(10, 20);
  slam::Map<int, Set> parent(&domain);
  auto makeSubset = [&parent]() {
    auto indices = RangeSetType(2, 5);
    return slam::SubMap<decltype(parent), RangeSetType>(&parent, indices);
  };
  auto subset = makeSubset();
  auto copy = subset;
  auto moved = std::move(copy);
  static_assert(std::same_as<typename decltype(subset)::ProjectedElement, double>);
  EXPECT_EQ(moved.index(1), 13.0);
  EXPECT_NE(subset.set(), moved.set());

  auto iterators = [&parent]() {
    slam::SubMap<decltype(parent), RangeSetType> local(&parent, RangeSetType(2, 5));
    return std::make_pair(local.set_begin(), local.set_end());
  }();
  EXPECT_EQ(iterators.first.index(), 12.0);
  --iterators.second;
  EXPECT_EQ(iterators.second.index(), 14.0);
}

TEST(slam_map, submap_values_and_indices_follow_parent_rebinding)
{
  using Set = typename slam::RangeSet<PositionType, double>::ConcreteSet;
  Set original(10, 20), replacement(30, 40);
  slam::Map<int, Set> parent(&original, 1);
  slam::SubMap<decltype(parent), RangeSetType> subset(&parent, RangeSetType(2, 5));
  slam::SubMap<decltype(subset), RangeSetType> nested(&subset, RangeSetType(1, 3));
  EXPECT_EQ(nested.index(0), 13.0);

  parent = decltype(parent)(&replacement, 9);
  EXPECT_EQ(nested(0), 9);
  EXPECT_EQ(nested.index(0), 33.0);
  EXPECT_EQ(nested.set_begin().index(), 33.0);

  Set smaller(30, 33);
  parent = decltype(parent)(&smaller, 4);
  EXPECT_FALSE(subset.isValid());
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
#ifdef AXOM_DEBUG
  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
#endif

  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  int result = RUN_ALL_TESTS();

  return result;
}
