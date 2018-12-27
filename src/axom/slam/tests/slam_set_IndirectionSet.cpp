// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


/**
 * \file slam_set_indirectionset.cpp
 */


#include <iterator>
#include <sstream>                  // for std::stringstream
#include <algorithm>                // for std::next_permutation
#include <iterator>                 // for std::ostream_iterator
#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/slam/Utilities.hpp"
#include "axom/slam/IndirectionSet.hpp"
#include "axom/slam/RangeSet.hpp"        // for PositionSet

namespace slam = axom::slam;

namespace
{
typedef slam::ArrayIndirectionSet SetType;
typedef SetType::PositionType SetPosition;
typedef SetType::ElementType SetElement;

static const SetPosition MAX_SET_SIZE = 10;

SetPosition* allocIncrementingArray(int size)
{
  SLIC_ASSERT(size > 0);

  SetPosition* arr = new SetPosition[size];

  for(int i = 0 ; i< size ; ++i)
    arr[i] = i;

  return arr;
}

void permuteArray(SetPosition* arr, int size)
{
  for(int i = 0 ; i<1000 ; ++i)
  {
    std::next_permutation(arr, arr + size);
  }
}
}

TEST(slam_set_indirectionset,construct)
{
  SLIC_INFO("Testing constructors for IndirectionSet");

  SetType s_empty;
  EXPECT_TRUE(s_empty.empty());
  EXPECT_EQ(0, s_empty.size() );
  EXPECT_TRUE(s_empty.isValid());


  SetType s(MAX_SET_SIZE);
  EXPECT_FALSE(s.empty());
  EXPECT_EQ(MAX_SET_SIZE, s.size() );
  EXPECT_FALSE(s.isValid())
    << "IndirectionSet not valid until we set the indirection buffer";

  s.data() = allocIncrementingArray(MAX_SET_SIZE);

  EXPECT_FALSE(s.empty());
  EXPECT_EQ(MAX_SET_SIZE, s.size() );
  EXPECT_TRUE(s.isValid());

  delete [] s.data();
}

TEST(slam_set_indirectionset,set_builder)
{
  SLIC_INFO("Testing construction of IndirectionSet using SetBuilders");

  typedef SetType::SetBuilder SetBuilder;

  bool bVerbose = true;

  {
    SetType s_empty(SetBuilder()
                    .size(MAX_SET_SIZE) );

    EXPECT_FALSE(s_empty.empty());
    EXPECT_EQ(MAX_SET_SIZE, s_empty.size() );
    EXPECT_FALSE(s_empty.isValid(bVerbose));
  }

  {
    SetType s_empty2(SetBuilder()
                     .data(allocIncrementingArray(MAX_SET_SIZE) ) );

    EXPECT_TRUE(s_empty2.empty());
    EXPECT_EQ(0, s_empty2.size() );
    EXPECT_TRUE(s_empty2.isValid(bVerbose));

    delete [] s_empty2.data();
  }

  {
    SetType s(SetBuilder()
              .size(MAX_SET_SIZE)
              .data(allocIncrementingArray(MAX_SET_SIZE) ) );

    EXPECT_FALSE(s.empty());
    EXPECT_EQ(MAX_SET_SIZE, s.size() );
    EXPECT_TRUE(s.isValid(bVerbose));

    delete [] s.data();
  }
}

TEST(slam_set_indirectionset,iterate)
{
  SetType s(SetType::SetBuilder()
            .size(MAX_SET_SIZE)
            .data(allocIncrementingArray(MAX_SET_SIZE) ) );

  SLIC_INFO("Iterating through set of size " << s.size());
  EXPECT_EQ(MAX_SET_SIZE, s.size());
  EXPECT_TRUE(s.isValid());

  SLIC_INFO("Using random access -- operator[] and at()");
  {
    std::stringstream sstr;
    for(SetPosition pos = SetPosition() ; pos < s.size() ; ++pos)
    {
      SetElement exp = static_cast<SetElement>(pos);
      EXPECT_EQ(exp,s[pos]);
      EXPECT_EQ(exp,s.at(pos));

      sstr << s[pos] << "\t";
    }
    SLIC_INFO("Data using operator[]:\t" << sstr.str());
  }

#ifdef AXOM_USE_CXX11
  SLIC_INFO("Using iterator interface");
  {
    std::stringstream sstr;
    typedef SetType::iterator SetIterator;
    for(SetIterator it = s.begin(), itEnd = s.end() ; it != itEnd ; ++it)
    {
      EXPECT_EQ( std::distance(s.begin(), it), *it );

      sstr << *it << "\t";
    }
    SLIC_INFO("Data from iterator:\t" << sstr.str());
  }

  const int NUM_PERMUTATIONS = 3;
  SLIC_INFO("Printing a few permutations of the data.");
  for(int i = 0 ; i< NUM_PERMUTATIONS ; ++i)
  {
    permuteArray(s.data(), s.size());
    std::stringstream sstr;
    std::copy(s.begin(), s.end(),
              std::ostream_iterator<SetPosition>(sstr, "\t") );
    SLIC_INFO("Permutation " << i << ":\t" << sstr.str());
  }
#endif

  delete [] s.data();
}

TEST(slam_set_indirectionset,equality)
{
  SLIC_INFO("Testing equality/inequality for several sets");

  SetType s1(SetType::SetBuilder()
             .size(MAX_SET_SIZE)
             .data(allocIncrementingArray(MAX_SET_SIZE) ) );
  EXPECT_TRUE(s1.isValid());

  // Test self-equality
  EXPECT_EQ(s1, s1);

  // Test equality against a set with the same array
  SetType s1b(SetType::SetBuilder()
              .size(MAX_SET_SIZE)
              .data(s1.data() ) );
  EXPECT_TRUE(s1b.isValid());
  EXPECT_EQ(s1, s1b);

  // Test equality against a set with a different array (same values)
  SetType s2(SetType::SetBuilder()
             .size(MAX_SET_SIZE)
             .data(allocIncrementingArray(MAX_SET_SIZE) ) );
  EXPECT_TRUE(s2.isValid());
  EXPECT_EQ(s1, s2);
  EXPECT_EQ(s2, s1);

  // Check that sets are not equal if we permute the data
  permuteArray(s2.data(), s2.size());
  EXPECT_NE(s1, s2);

  // Test against a smaller set
  const int S3_SIZE = MAX_SET_SIZE / 2;
  SetType s3(SetType::SetBuilder()
             .size(S3_SIZE)
             .data(allocIncrementingArray(S3_SIZE) ) );
  EXPECT_TRUE(s3.isValid());
  EXPECT_NE(s1, s3);
  EXPECT_NE(s2, s3);


  // Check against a PositionSet with the same elements
  slam::PositionSet posSet(MAX_SET_SIZE);
  EXPECT_TRUE(posSet.isValid() );
  EXPECT_EQ(posSet, s1);


  // Reclaim memory
  delete [] s1.data();
  delete [] s2.data();
  delete [] s3.data();
}


TEST(slam_set_indirectionset,out_of_bounds)
{
  SLIC_INFO("Testing out of bounds access on initialized set"
            <<"-- code is expected to assert and die.");

  SetType s(SetType::SetBuilder()
            .size(MAX_SET_SIZE)
            .data(allocIncrementingArray(MAX_SET_SIZE) ) );
  EXPECT_TRUE(s.isValid());

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode, so this test will only fail
  // in debug mode
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED( s[MAX_SET_SIZE], "");
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif

  delete [] s.data();
}


TEST(slam_set_indirectionset,vector_indirection)
{
  SLIC_INFO("Testing basic set operations on a VectorIndirectionSet");

  typedef slam::VectorIndirectionSet VectorIndirectionSet;
  typedef VectorIndirectionSet::SetBuilder Builder;

  // Set up data -- an array of incrementing integers
  std::vector<int> intVec;
  intVec.reserve(MAX_SET_SIZE);
  for(int i = 0 ; i<MAX_SET_SIZE ; ++i)
  {
    intVec.push_back(i);
  }

  // Setup the VectorIndirectionSet and test basic functionality
  VectorIndirectionSet vecSet(Builder()
                              .size(MAX_SET_SIZE)
                              .data(&intVec));
  EXPECT_TRUE(vecSet.isValid());
  EXPECT_FALSE(vecSet.empty());
  EXPECT_EQ(MAX_SET_SIZE, vecSet.size());
  EXPECT_TRUE(vecSet.hasIndirection());

  // Compare the values
  for(int i = 0 ; i<MAX_SET_SIZE ; ++i)
  {
    EXPECT_EQ(  i,  vecSet[i]);
    EXPECT_EQ(  i,  vecSet.at(i));
  }

  // Test equality with an array indirection set
  {
    SetType arrSet(SetType::SetBuilder()
                   .size(MAX_SET_SIZE)
                   .data(allocIncrementingArray(MAX_SET_SIZE) ) );
    EXPECT_TRUE(arrSet.isValid());

    EXPECT_EQ(vecSet, arrSet);
    EXPECT_EQ(arrSet, vecSet);

    delete [] arrSet.data();
  }

  // Test equality with a position set
  {
    slam::PositionSet posSet(MAX_SET_SIZE);
    EXPECT_TRUE(posSet.isValid() );
    EXPECT_EQ(vecSet, posSet);
    EXPECT_EQ(posSet, vecSet);
  }
}

TEST(slam_set_indirectionset,negative_stride)
{
  SLIC_INFO("Testing negative strides on an indirection set");

  namespace policies = axom::slam::policies;

  typedef slam::Set::PositionType SetPosition;

  typedef policies::RuntimeSize<SetPosition>                          SizePol;
  typedef policies::RuntimeOffset<SetPosition>                        OffPol;
  typedef policies::RuntimeStride<SetPosition>                        StridePol;
  typedef policies::STLVectorIndirection<SetPosition, SetPosition>    VecIndPol;
  typedef policies::ArrayIndirection<SetPosition, SetPosition>        ArrIndPol;

  typedef slam::OrderedSet<SetPosition, SizePol,OffPol,StridePol,
                           VecIndPol> VecSet;
  typedef VecSet::SetBuilder VecSetBuilder;

  typedef slam::OrderedSet<SetPosition, SizePol,OffPol,StridePol,
                           ArrIndPol> ArrSet;
  typedef ArrSet::SetBuilder ArrSetBuilder;

  // Set up data -- an array of incrementing integers
  std::vector<int> intVec;
  intVec.reserve(MAX_SET_SIZE);
  for(int i = 0 ; i<MAX_SET_SIZE ; ++i)
  {
    intVec.push_back(i);
  }

  const int setSize = MAX_SET_SIZE / 2;
  const int setOffset = MAX_SET_SIZE - 1;
  const int setStride = -2;

  // Setup the VectorIndirectionSet and test basic functionality
  {
    VecSet vSet(VecSetBuilder()
                .size(setSize)
                .offset(setOffset)
                .stride(setStride)
                .data(&intVec));

    EXPECT_TRUE(vSet.isValid());
    EXPECT_FALSE(vSet.empty());
    EXPECT_EQ(setSize, vSet.size());
    EXPECT_TRUE(vSet.hasIndirection());

    SLIC_INFO("Ordered vector set has:"
              << "\n\t --size "   << vSet.size()
              << "\n\t --stride " << vSet.stride()
              << "\n\t --offset " << vSet.offset()
              << "\n\t --first elt " << vSet[0]
              << "\n\t --last elt " << vSet[ vSet.size() - 1 ]
              );

    // Test the elements
    EXPECT_EQ(intVec[setOffset], vSet[0]);
    for(int i = 0 ; i< vSet.size() ; ++i)
    {
      EXPECT_EQ( setOffset + setStride * i, vSet[i] );
    }

    /// Several tests to check that sets with bad offsets and strides are not
    // valid

    VecSet noDataVSet(VecSetBuilder()  // Note: Missing a data pointer
                      .size(setSize)
                      .offset(setOffset)
                      .stride(setStride));
    EXPECT_FALSE(noDataVSet.isValid(true));

    VecSet outOfBoundsVSet(VecSetBuilder()
                           .size(setSize + 1) // Note: This will cause the last
                                              // index to be out of bounds
                           .offset(setOffset)
                           .stride(setStride)
                           .data(&intVec));
    EXPECT_FALSE(outOfBoundsVSet.isValid(true));

    VecSet outOfBoundsVSet2(VecSetBuilder()
                            .size(setSize)
                            .offset(-1) // Note: This will cause the first index
                                        // to be out of bounds
                            .stride(1)
                            .data(&intVec));
    EXPECT_FALSE(outOfBoundsVSet2.isValid(true));

    VecSet zeroStrideVSet(VecSetBuilder()
                          .size(setSize)
                          .offset(setOffset)
                          .stride(0) // Note: A stride of zero is not valid
                          .data(&intVec));
    EXPECT_FALSE(zeroStrideVSet.isValid(true));
  }

  // Setup the VectorIndirectionSet and test basic functionality
  {
    ArrSet aSet(ArrSetBuilder()
                .size(setSize)
                .offset(setOffset)
                .stride(setStride)
                .data(&intVec[0]));

    EXPECT_TRUE(aSet.isValid());
    EXPECT_FALSE(aSet.empty());
    EXPECT_EQ(setSize, aSet.size());
    EXPECT_TRUE(aSet.hasIndirection());

    EXPECT_EQ(intVec[setOffset], aSet[0]);

    SLIC_INFO("Ordered array set has:"
              << "\n\t --size "   << aSet.size()
              << "\n\t --stride " << aSet.stride()
              << "\n\t --offset " << aSet.offset()
              << "\n\t --first elt " << aSet[0]
              << "\n\t --last elt " << aSet[ aSet.size() - 1 ]
              );

    for(int i = 0 ; i< aSet.size() ; ++i)
    {
      EXPECT_EQ( setOffset + setStride * i, aSet[i] );
    }

    /// Several tests to check that sets with bad offsets and strides are not
    // valid

    ArrSet noDataASet(ArrSetBuilder()  // Note: Missing a data pointer
                      .size(setSize)
                      .offset(setOffset)
                      .stride(setStride));
    EXPECT_FALSE(noDataASet.isValid(true));

    ArrSet outOfBoundsASet1(ArrSetBuilder()
                            .size(setSize + 1) // Note: This will cause the last
                                               // index to be out of bounds
                            .offset(setOffset)
                            .stride(setStride)
                            .data(&intVec[0]));
    EXPECT_FALSE(outOfBoundsASet1.isValid(true));

    ArrSet outOfBoundsASet2(ArrSetBuilder()
                            .size(setSize)
                            .offset(-1) // Note: This will cause the first index
                                        // to be out of bounds
                            .stride(1)
                            .data(&intVec[0]));
    EXPECT_FALSE(outOfBoundsASet2.isValid(true));

    ArrSet zeroStrideASet(ArrSetBuilder()
                          .size(setSize)
                          .offset(setOffset)
                          .stride(0) // Note: A stride of zero is not valid
                          .data(&intVec[0]));
    EXPECT_FALSE(zeroStrideASet.isValid(true));
  }

}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  // create & initialize test logger. finalized when exiting main scope
  axom::slic::UnitTestLogger logger;

  int result = RUN_ALL_TESTS();

  return result;
}
