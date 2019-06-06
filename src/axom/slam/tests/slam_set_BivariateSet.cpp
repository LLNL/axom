// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 * \file slam_set_BivariateSet.cpp
 *
 * This file tests BivariateSet, ProductSet and RelationSet within Slam.
 */

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"

#include <type_traits>
#include <sstream>
#include <iostream>

namespace slam = axom::slam;
namespace policies = axom::slam::policies;

namespace
{

static const int FIRST_SET_SIZE = 5;
static const int SECOND_SET_SIZE = 16;

static const int FIRST_SET_OFFSET = 3;
static const int SECOND_SET_OFFSET = 2;

// Template aliases to simplify specifying some sets and relations
template<typename P, typename E, typename FromSet, typename ToSet>
using RelType =
        slam::StaticRelation<P,E,
                             policies::VariableCardinality<E>,
                             policies::STLVectorIndirection<P,E>,
                             FromSet, ToSet>;

template<typename SetType>
using PositionSetType = slam::PositionSet<typename SetType::PositionType,
                                          typename SetType::ElementType>;

template<typename SetType>
using RangeSetType = slam::RangeSet<typename SetType::PositionType,
                                    typename SetType::ElementType>;


// Factory function for SetType == slam::Set
// Caller is responsible for deallocating the returned object
template<typename SetType>
typename std::enable_if<std::is_abstract<SetType>::value,
                        SetType*>::type
constructSet(int sz)
{
  using P = typename SetType::PositionType;
  using E = typename SetType::ElementType;
  return new slam::PositionSet<P,E>(sz);
}

// Factory function for SetType == slam::PositionSet
// Caller is responsible for deallocating the returned object
template<typename SetType>
typename std::enable_if<std::is_same<SetType, PositionSetType<SetType> >::value,
                        SetType*>::type
constructSet(int sz)
{
  using P = typename SetType::PositionType;
  using E = typename SetType::ElementType;

  return new slam::PositionSet<P,E>(sz);
}


// Factory function for SetType == slam::RangeSet
// Caller is responsible for deallocating the returned object
template<typename SetType>
typename std::enable_if<std::is_same<SetType, RangeSetType<SetType> >::value,
                        SetType*>::type
constructSet(int sz)
{
  using P = typename SetType::PositionType;
  using E = typename SetType::ElementType;

  const int offset =
    (sz==FIRST_SET_SIZE) ? FIRST_SET_OFFSET : SECOND_SET_OFFSET;
  return new slam::RangeSet<P,E>(offset, sz+offset);
}

// Predicate to check if b%a is zero
// Note: Also returns true when a is zero to avoid division by zero
bool modCheck(int a, int b)
{
  return (a == 0)  || (b%a == 0);
}

} // end anonymous namespace


// Typed test fixture for testing slam's BivariateSet class
template<typename BSet>
class BivariateSetTester : public ::testing::Test
{
public:
  using BivariateSetType = BSet;
  using FirstSetType = typename BSet::FirstSetType;
  using SecondSetType = typename BSet::SecondSetType;
  using PositionType = typename BSet::PositionType;
  using ElementType = typename BSet::ElementType;

  using Vec = std::vector<PositionType>;
  using RelationType = ::RelType<PositionType,ElementType,
                                 FirstSetType,SecondSetType>;

  BivariateSetTester()
    : m_set1(nullptr)
    , m_set2(nullptr)
  {}

  virtual void SetUp()
  {
    m_set1 = constructSet<FirstSetType>(FIRST_SET_SIZE);
    EXPECT_TRUE(m_set1->isValid(true));

    m_set2 = constructSet<SecondSetType>(SECOND_SET_SIZE);
    EXPECT_TRUE(m_set2->isValid(true));

    initRelationData();
    constructRelation();
  }

  virtual void TearDown()
  {
    delete m_set1;
    m_set1 = nullptr;

    delete m_set2;
    m_set2 = nullptr;
  }

  void initRelationData()
  {
    relationBegins.clear();
    relationIndices.clear();

    int cnt = 0;
    for(int i = 0 ; i < m_set1->size() ; ++i)
    {
      auto outer = m_set1->at(i);
      relationBegins.push_back(cnt);

      for(int j=0 ; j < m_set2->size() ; ++j)
      {
        auto inner = m_set2->at(j);
        if( modCheck(outer, inner) )
        {
          relationIndices.push_back( j );
          ++cnt;
        }
      }
    }
    relationBegins.push_back(cnt);
  }

  void constructRelation()
  {
    using RelationBuilder = typename RelationType::RelationBuilder;
    modRelation = RelationBuilder()
                  .fromSet( m_set1 )
                  .toSet( m_set2 )
                  .begins( typename RelationBuilder::BeginsSetBuilder()
                           .size( relationBegins.size() )
                           .data( &relationBegins ) )
                  .indices( typename RelationBuilder::IndicesSetBuilder()
                            .size( relationIndices.size() )
                            .data( &relationIndices ) );

    EXPECT_TRUE(modRelation.isValid(true));
  }

  FirstSetType* m_set1;
  SecondSetType* m_set2;

  Vec relationBegins;
  Vec relationIndices;
  RelationType modRelation;
};


// Tests several types of BivariateSet that differ in their underlying set types
using MyTypes = ::testing::Types <
        slam::BivariateSet< >
        , slam::BivariateSet< slam::PositionSet<>, slam::PositionSet<> >
        , slam::BivariateSet< slam::PositionSet<>, slam::RangeSet<> >
        , slam::BivariateSet< slam::RangeSet<>, slam::PositionSet<> >
        , slam::BivariateSet< slam::RangeSet<>, slam::RangeSet<> >
        , slam::BivariateSet< slam::Set<>, slam::RangeSet<> >
        , slam::BivariateSet< slam::PositionSet<>, slam::Set<> >
        >;
TYPED_TEST_CASE(BivariateSetTester, MyTypes);

//-----------------------------------------------------------------------------

// Tests BivariateSeet::firstSetSize(), secondSetSize()
template<typename S1, typename S2>
void bSetSizesTest(slam::BivariateSet<S1,S2>* bset)
{
  const auto firstSetSize = bset->firstSetSize();
  const auto secondSetSize = bset->secondSetSize();

  EXPECT_EQ(FIRST_SET_SIZE, firstSetSize);
  EXPECT_EQ(SECOND_SET_SIZE, secondSetSize);

  EXPECT_EQ(firstSetSize, bset->getFirstSet()->size() );
  EXPECT_EQ(secondSetSize, bset->getSecondSet()->size() );
}

// Tests BivariateSeet::getFirstSet(), getSecondSet(), getElements(idx)
// and the ability to operate and iterate on the resulting sets
template<typename S1, typename S2>
void bSetTraverseTest(slam::BivariateSet<S1,S2>* bset, bool shouldCheckMod)
{
  const auto firstSetSize = bset->firstSetSize();

  const auto* firstSet = bset->getFirstSet();
  const auto* secondSet = bset->getSecondSet();

  for( int idx =0 ; idx < firstSetSize ; ++idx )
  {
    auto outer = firstSet->at(idx);

    std::stringstream sstr;
    sstr << idx << ": " << outer << " -> { ";

    auto elems = bset->getElements(idx);
    for(auto innerIdx : elems )
    {
      auto inner = secondSet->at(innerIdx);

      if(shouldCheckMod)
      {
        EXPECT_TRUE( modCheck(outer, inner));
      }

      sstr << "(" << outer <<"," << inner << ") ";
    }
    SLIC_INFO(sstr.str() << "}");
  }
}

TYPED_TEST(BivariateSetTester, smoke)
{
  using BSet = typename TestFixture::BivariateSetType;
  using S1 = typename TestFixture::FirstSetType;
  using S2 = typename TestFixture::SecondSetType;

  // Test a default-constructed BivariateSet
  {
    BSet* bset = new slam::ProductSet<S1,S2>();

    // We expect this object to be valid when the types of both S1 and S2 are
    // abstract slam::Sets, but not if either is specialized.
    bool isSet1Abstract = std::is_abstract<S1>::value;
    bool isSet2Abstract = std::is_abstract<S2>::value;

    EXPECT_EQ(isSet1Abstract && isSet2Abstract, bset->isValid(true));
    delete bset;
  }

  // Test a non-default-constructed BivariateSet
  {
    BSet* bset = new slam::ProductSet<S1,S2>(this->m_set1, this->m_set2);
    EXPECT_TRUE(bset->isValid(true));
    delete bset;
  }
}

TYPED_TEST(BivariateSetTester, sizes)
{
  using S1 = typename TestFixture::FirstSetType;
  using S2 = typename TestFixture::SecondSetType;

  // Test the size functions for a ProductSet
  {
    using PSet = slam::ProductSet<S1,S2>;

    PSet pset = PSet(this->m_set1, this->m_set2);
    EXPECT_TRUE(pset.isValid(true));

    bSetSizesTest<S1,S2>(&pset);
  }

  // Test the size functions for a RelationSet
  {
    using RType = typename TestFixture::RelationType;
    using RSet = slam::RelationSet<RType,S1,S2>;

    RSet rset = RSet(&this->modRelation);
    EXPECT_TRUE(rset.isValid(true));

    bSetSizesTest<S1,S2>(&rset);
  }

}

TYPED_TEST(BivariateSetTester, traverse)
{
  using S1 = typename TestFixture::FirstSetType;
  using S2 = typename TestFixture::SecondSetType;

  // Test traversal functions for a ProductSet
  SLIC_INFO("Traversing product set"
            << "\n       ----------------------");
  {
    using PSet = slam::ProductSet<S1,S2>;

    PSet pset = PSet(this->m_set1, this->m_set2);
    EXPECT_TRUE(pset.isValid(true));

    bool checkMod = false;
    bSetTraverseTest<S1,S2>(&pset, checkMod);
  }

  // Test traversal functions for a RelationSet
  std::cout<<std::endl;
  SLIC_INFO("Traversing relation set"
            << "\n       -----------------------");
  {
    using RType = typename TestFixture::RelationType;
    using RSet = slam::RelationSet<RType,S1,S2>;

    RSet rset = RSet(&this->modRelation);
    EXPECT_TRUE(rset.isValid(true));

    bool checkMod = true;
    bSetTraverseTest<S1,S2>(&rset, checkMod);
  }

}

//-----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  // create & initialize test logger. finalized when exiting main scope
  axom::slic::UnitTestLogger logger;

  int result = RUN_ALL_TESTS();

  return result;
}
