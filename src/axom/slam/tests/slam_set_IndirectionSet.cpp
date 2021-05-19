// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 * \file slam_set_indirectionset.cpp
 *
 * Most unit tests in this file are parametrized on an Indirection set type
 * using the IndirectionSetTester class
 */

#include <iterator>
#include <sstream>    // for std::stringstream
#include <algorithm>  // for std::next_permutation
#include <iterator>   // for std::ostream_iterator
#include "gtest/gtest.h"

#include "axom/slic.hpp"

#include "axom/slam/Utilities.hpp"
#include "axom/slam/IndirectionSet.hpp"
#include "axom/slam/RangeSet.hpp"  // for PositionSet

namespace slam = axom::slam;

namespace
{
static const int MAX_SET_SIZE = 10;

}  // end anonymous namespace

/**
 * Simple test class for IndirectionSets
 * Initializes three vectors of data
 *
 * \note Uses specialized variants of getDataBuffer() to access
 * the data in an appropriated way for the indirectionSet \a TheSet
 */
template <typename TheSet>
class IndirectionSetTester : public ::testing::Test
{
public:
  using SetType = TheSet;
  using IndirectionType = typename SetType::IndirectionPolicyType;
  using BufferType = typename IndirectionType::IndirectionBufferType;

  using PosType = typename SetType::PositionType;
  using ElemType = typename SetType::ElementType;

  virtual void SetUp()
  {
    // Allocate and initialize three incrementing vectors
    // vecs 0 and 1 have the same size; vec 2 is shorter

    mDataVec[0].reserve(MAX_SET_SIZE);
    mDataVec[1].reserve(MAX_SET_SIZE);
    mDataVec[2].reserve(MAX_SET_SIZE / 2);

    for(int i = 0; i < MAX_SET_SIZE; ++i)
    {
      mDataVec[0].push_back(static_cast<ElemType>(i));
      mDataVec[1].push_back(static_cast<ElemType>(i));
    }
    for(int i = 0; i < MAX_SET_SIZE / 2; ++i)
    {
      mDataVec[2].push_back(static_cast<ElemType>(i));
    }
  }

  /**
   * Permute the data buffer with index \a index
   * \pre index must be 0, 1 or 2
   */
  void permuteData(int index = 0)
  {
    checkIndex(index);

    for(int i = 0; i < 1000; ++i)
    {
      std::next_permutation(mDataVec[index].begin(), mDataVec[index].end());
    }
  }

  /**
   * Get a pointer to the data buffer
   * Specialization for ArrayIndirection
   * \pre index must be 0, 1 or 2
   */
  void getDataBuffer(ElemType*& ptr, int index = 0)
  {
    checkIndex(index);

    ptr = mDataVec[index].data();
  }

  /**
   * Get a pointer to the data buffer
   * Specialization for VectorIndirection
   * \pre index must be 0, 1 or 2
   */
  void getDataBuffer(std::vector<ElemType>*& ptr, int index = 0)
  {
    checkIndex(index);

    ptr = &mDataVec[index];
  }

private:
  /** Check that the index is valid; fail test if not */
  void checkIndex(int index)
  {
    ASSERT_GE(index, 0);
    ASSERT_LT(index, 3);
  }

private:
  std::vector<ElemType> mDataVec[3];
};

// Tests several types of indirection sets
using MyTypes =
  ::testing::Types<slam::ArrayIndirectionSet<axom::int32, axom::int64>,
                   slam::VectorIndirectionSet<axom::int32, axom::int64>>;

TYPED_TEST_SUITE(IndirectionSetTester, MyTypes);

TYPED_TEST(IndirectionSetTester, constuct)
{
  SLIC_INFO("Testing simple indirection Set constructors");

  using SetType = typename TestFixture::SetType;

  // Empty set
  {
    SetType s_empty;
    EXPECT_TRUE(s_empty.empty());
    EXPECT_EQ(0, s_empty.size());
    EXPECT_TRUE(s_empty.isValid());
    EXPECT_FALSE(s_empty.hasIndirection());
  }

  // Incrementally build the set
  {
    SetType s(MAX_SET_SIZE);

    // set is originally empty
    EXPECT_FALSE(s.empty());
    EXPECT_EQ(MAX_SET_SIZE, s.size());
    EXPECT_FALSE(s.hasIndirection());
    EXPECT_FALSE(s.isValid())
      << "IndirectionSet not valid until we set the indirection buffer";

    // Add data to the set
    this->getDataBuffer(s.data());

    EXPECT_FALSE(s.empty());
    EXPECT_EQ(MAX_SET_SIZE, s.size());
    EXPECT_TRUE(s.hasIndirection());
    EXPECT_TRUE(s.isValid());
  }
}

TYPED_TEST(IndirectionSetTester, set_builder)
{
  SLIC_INFO("Testing construction of IndirectionSet using SetBuilders");

  using SetType = typename TestFixture::SetType;
  using BufferType = typename TestFixture::BufferType;
  using SetBuilder = typename SetType::SetBuilder;

  BufferType* bufPtr = nullptr;
  this->getDataBuffer(bufPtr);

  bool bVerbose = false;

  // Only set the size, but not the data
  {
    SetType s_empty(SetBuilder().size(MAX_SET_SIZE));

    EXPECT_FALSE(s_empty.empty());
    EXPECT_EQ(MAX_SET_SIZE, s_empty.size());
    EXPECT_FALSE(s_empty.hasIndirection());
    EXPECT_FALSE(s_empty.isValid(bVerbose));
  }

  // Set the data but not the size
  {
    SetType s_empty2(SetBuilder().data(bufPtr));

    EXPECT_TRUE(s_empty2.empty());
    EXPECT_EQ(0, s_empty2.size());
    EXPECT_TRUE(s_empty2.hasIndirection());
    EXPECT_TRUE(s_empty2.isValid(bVerbose));
  }

  // Set the data and the size
  {
    SetType s(SetBuilder().size(MAX_SET_SIZE).data(bufPtr));

    EXPECT_FALSE(s.empty());
    EXPECT_EQ(MAX_SET_SIZE, s.size());
    EXPECT_EQ(0, s.offset());
    EXPECT_EQ(1, s.stride());
    EXPECT_TRUE(s.hasIndirection());
    EXPECT_TRUE(s.isValid(bVerbose));
  }

  // Set the data and the size
  // also set the offset to 0 and the stride to 1
  {
    SetType s(SetBuilder().size(MAX_SET_SIZE).offset(0).stride(1).data(bufPtr));

    EXPECT_FALSE(s.empty());
    EXPECT_EQ(MAX_SET_SIZE, s.size());
    EXPECT_EQ(0, s.offset());
    EXPECT_EQ(1, s.stride());
    EXPECT_TRUE(s.hasIndirection());
    EXPECT_TRUE(s.isValid(bVerbose));
  }
}

TYPED_TEST(IndirectionSetTester, iterate)
{
  using SetType = typename TestFixture::SetType;
  using BufferType = typename TestFixture::BufferType;
  using SetPosition = typename TestFixture::PosType;
  using SetElement = typename TestFixture::ElemType;

  BufferType* bufPtr = nullptr;
  this->getDataBuffer(bufPtr);

  SetType s(typename SetType::SetBuilder().size(MAX_SET_SIZE).data(bufPtr));

  SLIC_INFO("Iterating through set of size " << s.size());
  EXPECT_EQ(MAX_SET_SIZE, s.size());
  EXPECT_TRUE(s.isValid());

  SLIC_INFO("Using random access -- operator[] and at()");
  {
    std::stringstream sstr;
    for(SetPosition pos = SetPosition(); pos < s.size(); ++pos)
    {
      SetElement exp = static_cast<SetElement>(pos);
      EXPECT_EQ(exp, s[pos]);
      EXPECT_EQ(exp, s.at(pos));

      sstr << s[pos] << "\t";
    }
    SLIC_INFO("Data using operator[]:\t" << sstr.str());

    // using range for over OrderedSet::positions()
    int count = 0;
    for(auto pos : s.positions())
    {
      SetElement exp = static_cast<SetElement>(pos);
      EXPECT_EQ(exp, s[pos]);
      EXPECT_EQ(exp, s.at(pos));
      ++count;
    }
    EXPECT_EQ(s.size(), count);
  }

  SLIC_INFO("Using iterator interface");
  {
    // also tests default .ctor and operator=()

    std::stringstream sstr;
    int count = 0;
    typename SetType::iterator it, itEnd;
    EXPECT_EQ(it, itEnd);

    for(it = s.begin(), itEnd = s.end(); it != itEnd; ++it)
    {
      EXPECT_EQ(std::distance(s.begin(), it), *it);
      ++count;

      sstr << *it << "\t";
    }
    EXPECT_EQ(s.size(), count);
    SLIC_INFO("Data from iterator:\t" << sstr.str());
  }

  const int NUM_PERMUTATIONS = 3;
  SLIC_INFO("Printing a few permutations of the data.");
  for(int i = 0; i < NUM_PERMUTATIONS; ++i)
  {
    this->permuteData();
    std::stringstream sstr;
    std::copy(s.begin(),
              s.end(),
              std::ostream_iterator<SetPosition>(sstr, "\t"));
    SLIC_INFO("Permutation " << i << ":\t" << sstr.str());
  }
}

TYPED_TEST(IndirectionSetTester, equality)
{
  SLIC_INFO("Testing equality/inequality for several sets");

  using SetType = typename TestFixture::SetType;
  using BufferType = typename TestFixture::BufferType;
  using SetPosition = typename TestFixture::PosType;
  using SetElement = typename TestFixture::ElemType;

  // Get pointers to the data
  BufferType* bufPtr[3] = {nullptr, nullptr, nullptr};
  this->getDataBuffer(bufPtr[0], 0);
  this->getDataBuffer(bufPtr[1], 1);
  this->getDataBuffer(bufPtr[2], 2);
  EXPECT_NE(bufPtr[0], bufPtr[1]);
  EXPECT_NE(bufPtr[0], bufPtr[2]);
  EXPECT_NE(bufPtr[1], bufPtr[2]);

  // Initialize the first set
  SetType s0(typename SetType::SetBuilder()  //
               .size(MAX_SET_SIZE)           //
               .data(bufPtr[0]));

  EXPECT_TRUE(s0.isValid());

  // Test self-equality
  EXPECT_EQ(s0, s0);

  // Test equality against a set with the same array
  SetType s0_b(typename SetType::SetBuilder()  //
                 .size(MAX_SET_SIZE)           //
                 .data(s0.data()));
  EXPECT_TRUE(s0_b.isValid());
  EXPECT_EQ(s0, s0_b);
  EXPECT_EQ(s0.data(), s0_b.data());

  // Test equality against a set with a different array (same values)

  SetType s1(typename SetType::SetBuilder()  //
               .size(MAX_SET_SIZE)           //
               .data(bufPtr[1]));
  EXPECT_TRUE(s1.isValid());
  EXPECT_EQ(s0, s1);
  EXPECT_EQ(s1, s0);

  // Check that sets are not equal if we permute the data
  this->permuteData(1);
  EXPECT_NE(s0, s1);

  // Test against a smaller set
  SetType s2(typename SetType::SetBuilder()  //
               .size(MAX_SET_SIZE / 2)       //
               .data(bufPtr[2]));
  EXPECT_TRUE(s2.isValid());
  EXPECT_NE(s0, s2);
  EXPECT_NE(s1, s2);

  // Check against a PositionSet with the same elements
  slam::PositionSet<SetPosition, SetElement> posSet(MAX_SET_SIZE);
  EXPECT_TRUE(posSet.isValid());
  EXPECT_EQ(posSet, s0);
}

TYPED_TEST(IndirectionSetTester, out_of_bounds)
{
  SLIC_INFO("Testing out of bounds access on initialized set"
            << "-- code is expected to assert and die.");

  using SetType = typename TestFixture::SetType;
  using BufferType = typename TestFixture::BufferType;

  BufferType* bufPtr = nullptr;
  this->getDataBuffer(bufPtr);

  SetType s(typename SetType::SetBuilder()  //
              .size(MAX_SET_SIZE)           //
              .data(bufPtr));

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode,
  // so this test will only fail in debug mode
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED(s[MAX_SET_SIZE], "");
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif
}

TEST(slam_set_indirectionset, compare_array_and_vector)
{
  SLIC_INFO("Comparing Array and Vector indirection Sets");

  using P = int;
  using E = int;

  using VecSet = slam::VectorIndirectionSet<P, E>;
  using ArrSet = slam::ArrayIndirectionSet<P, E>;
  using PosSet = slam::PositionSet<P, E>;

  const int SZ = MAX_SET_SIZE;
  std::vector<E> vVals(SZ);
  std::vector<E> aVals(SZ);

  PosSet ps(SZ);

  // Create vector and array sets without initializing their elements
  VecSet vs(VecSet::SetBuilder()  //
              .size(SZ)           //
              .data(&vVals));

  ArrSet as(ArrSet::SetBuilder()  //
              .size(SZ)           //
              .data(aVals.data()));

  EXPECT_NE(vs, ps);
  EXPECT_NE(as, ps);

  // Use position set to initialize vector and array set elements
  for(auto i : ps)
  {
    vs[i] = i;
    as[i] = i;
  }

  EXPECT_EQ(vs, as);
  EXPECT_EQ(as, ps);
  EXPECT_EQ(ps, vs);
}

TEST(slam_set_indirectionset, negative_stride)
{
  SLIC_INFO("Testing negative strides on an indirection set");

  namespace policies = axom::slam::policies;

  using SetPosition = slam::DefaultPositionType;
  using SetElement = slam::DefaultElementType;

  using SizePol = policies::RuntimeSize<SetPosition>;
  ;
  using OffPol = policies::RuntimeOffset<SetPosition>;
  using StridePol = policies::RuntimeStride<SetPosition>;
  using VecIndPol = policies::STLVectorIndirection<SetPosition, SetElement>;
  using ArrIndPol = policies::ArrayIndirection<SetPosition, SetElement>;

  using VecSet =
    slam::OrderedSet<SetPosition, SetElement, SizePol, OffPol, StridePol, VecIndPol>;

  using ArrSet =
    slam::OrderedSet<SetPosition, SetElement, SizePol, OffPol, StridePol, ArrIndPol>;

  // Set up data -- an array of incrementing integers
  std::vector<SetElement> intVec(MAX_SET_SIZE);
  for(auto i : slam::PositionSet<>(MAX_SET_SIZE))
  {
    intVec[i] = i;
  }

  const int setSize = MAX_SET_SIZE / 2;
  const int setOffset = MAX_SET_SIZE - 1;
  const int setStride = -2;
  const bool bVerbose = false;

  // Setup the VectorIndirectionSet and test basic functionality
  VecSet vSet(VecSet::SetBuilder()  //
                .size(setSize)      //
                .offset(setOffset)  //
                .stride(setStride)  //
                .data(&intVec));
  {
    EXPECT_TRUE(vSet.isValid());
    EXPECT_FALSE(vSet.empty());
    EXPECT_EQ(setSize, vSet.size());
    EXPECT_TRUE(vSet.hasIndirection());

    SLIC_INFO("Ordered vector set has:"
              << "{ size: " << vSet.size() << ", stride: " << vSet.stride()
              << ", offset: " << vSet.offset() << ", first elt: " << vSet[0]
              << ", last elt: " << vSet[vSet.size() - 1] << "}");

    // Test the elements
    EXPECT_EQ(intVec[setOffset], vSet[0]);
    for(int i = 0; i < vSet.size(); ++i)
    {
      EXPECT_EQ(setOffset + setStride * i, vSet[i]);
    }

    /// Several checks that sets with bad offsets and strides are invalid
    SLIC_DEBUG_IF(bVerbose,
                  "--- Checking isValid() on several sets with "
                    << "bad sizes, offsets and strides.");
    VecSet noDataVSet(VecSet::SetBuilder()  // Note: Missing a data pointer
                        .size(setSize)
                        .offset(setOffset)
                        .stride(setStride));
    EXPECT_FALSE(noDataVSet.isValid(bVerbose));

    VecSet outOfBoundsVSet(VecSet::SetBuilder()
                             .size(setSize + 1)  // Note: This will cause the last
                                                 // index to be out of bounds
                             .offset(setOffset)
                             .stride(setStride)
                             .data(&intVec));
    EXPECT_FALSE(outOfBoundsVSet.isValid(bVerbose));

    VecSet outOfBoundsVSet2(VecSet::SetBuilder()
                              .size(setSize)
                              .offset(-1)  // Note: This will cause the first index
                                           // to be out of bounds
                              .stride(1)
                              .data(&intVec));
    EXPECT_FALSE(outOfBoundsVSet2.isValid(bVerbose));

    VecSet zeroStrideVSet(VecSet::SetBuilder()
                            .size(setSize)
                            .offset(setOffset)
                            .stride(0)  // Note: A stride of zero is not valid
                            .data(&intVec));
    EXPECT_FALSE(zeroStrideVSet.isValid(bVerbose));

    SLIC_DEBUG_IF(bVerbose, "--- Done.");
  }

  // Setup the VectorIndirectionSet and test basic functionality
  ArrSet aSet(ArrSet::SetBuilder()  //
                .size(setSize)      //
                .offset(setOffset)  //
                .stride(setStride)  //
                .data(intVec.data()));
  {
    EXPECT_TRUE(aSet.isValid());
    EXPECT_FALSE(aSet.empty());
    EXPECT_EQ(setSize, aSet.size());
    EXPECT_TRUE(aSet.hasIndirection());

    EXPECT_EQ(intVec[setOffset], aSet[0]);

    SLIC_INFO("Ordered array set has:"
              << "{ size: " << aSet.size() << ", stride: " << aSet.stride()
              << ", offset: " << aSet.offset() << ", first elt: " << aSet[0]
              << ", last elt: " << aSet[aSet.size() - 1] << "}");

    for(int i = 0; i < aSet.size(); ++i)
    {
      EXPECT_EQ(setOffset + setStride * i, aSet[i]);
    }

    /// Several checks that sets with bad offsets and strides are invalid
    SLIC_DEBUG_IF(bVerbose,
                  "--- Checking isValid() on several sets with "
                    << "bad sizes, offsets and strides.");

    ArrSet noDataASet(ArrSet::SetBuilder()  // Note: Missing a data pointer
                        .size(setSize)
                        .offset(setOffset)
                        .stride(setStride));
    EXPECT_FALSE(noDataASet.isValid(bVerbose));

    ArrSet outOfBoundsASet1(ArrSet::SetBuilder()
                              .size(setSize + 1)  // Note: This will cause the last
                                                  // index to be out of bounds
                              .offset(setOffset)
                              .stride(setStride)
                              .data(intVec.data()));
    EXPECT_FALSE(outOfBoundsASet1.isValid(bVerbose));

    ArrSet outOfBoundsASet2(ArrSet::SetBuilder()
                              .size(setSize)
                              .offset(-1)  // Note: This will cause the first index
                                           // to be out of bounds
                              .stride(1)
                              .data(intVec.data()));
    EXPECT_FALSE(outOfBoundsASet2.isValid(bVerbose));

    ArrSet zeroStrideASet(ArrSet::SetBuilder()
                            .size(setSize)
                            .offset(setOffset)
                            .stride(0)  // Note: A stride of zero is not valid
                            .data(intVec.data()));
    EXPECT_FALSE(zeroStrideASet.isValid(bVerbose));

    SLIC_DEBUG_IF(bVerbose, "--- Done.");
  }

  // check that vset and aset are equivalent
  EXPECT_EQ(vSet, aSet);
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  // create & initialize test logger. finalized when exiting main scope
  axom::slic::SimpleLogger logger;

  int result = RUN_ALL_TESTS();

  return result;
}
