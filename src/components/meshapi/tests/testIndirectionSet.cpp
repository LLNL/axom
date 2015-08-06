/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/**
 * \file
 */

#include <iterator>
#include "gtest/gtest.h"

#include "meshapi/Utilities.hpp"
#include "meshapi/IndirectionSet.hpp"

typedef asctoolkit::meshapi::ArrayIndirectionSet SetType;
typedef SetType::iterator                   SetIterator;
typedef SetType::PositionType               SetPosition;
typedef SetType::ElementType                SetElement;

static const SetPosition MAX_SET_SIZE = 10;


TEST(gtest_meshapi_indirection_set,construct_indirection_set)
{

  SetType s(MAX_SET_SIZE);
  SetPosition arr[MAX_SET_SIZE];
  s.data() = arr;
  for(int i=0;i<MAX_SET_SIZE; ++i)
      arr[i] = i;


  EXPECT_TRUE(s.isValid());

    if(MAX_SET_SIZE > SetPosition())
        EXPECT_FALSE(s.empty());
/*
    std::cout<<"Iterating through set of size " << s.size() << std::endl;
    EXPECT_EQ(s.size(), MAX_SET_SIZE);


    std::cout<<"\n --Using begin/end" << std::endl;
//    for(SetIterator it=s.begin(), itEnd=s.end(); it != itEnd; ++it)
//    {
//        EXPECT_EQ( std::distance(s.begin(), it), *it )
//                << "Iterator dereference should be equal to its position in the set";
//        std::cout << "\t" << *it <<"\n";
//    }

    std::cout<<"\n --Using random access -- operator[]" << std::endl;
    for(SetPosition pos = SetPosition(); pos < s.size(); ++pos)
    {
        SetElement elt = static_cast<SetElement>(pos);
        EXPECT_EQ(elt,s[pos])
                <<"Random access iterator dereference to equal its position in the set";
        std::cout << "\t" << s[pos] <<"\n";
    }

    std::cout<<"\n --Using checked random access -- at()" << std::endl;
    for(SetPosition pos = SetPosition(); pos < s.size(); ++pos)
    {
        SetElement elt = static_cast<SetElement>(pos);
        EXPECT_EQ(elt,s.at(pos))
                <<"Expected checked random access iterator dereference to equal its position in the set";

        std::cout << "\t" << s.at(pos) <<"\n";
    }

   #ifdef ATK_DEBUG
    // NOTE: ATK_ASSSERT is disabled in release mode, so this test will only fail in debug mode

    std::cout<<"\n --Using checked random access -- at() with invalid address" << std::endl;

    // add this line to avoid a warning in the output about thread safety
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    ASSERT_DEATH(s.at(MAX_SET_SIZE),"") << "tried to access out of range element";
   #else
    std::cout <<"Did not check for assertion failure since assertions are compiled out in release mode." << std::endl;
   #endif


  std::cout << "--\ndone." << std::endl;
*/
}


   TEST(gtest_meshapi_range_set,test_range_set_out_of_bounds)
   {
/*
    std::cout<<"\n****** Testing out of bounds access on initialized set-- code is expected to assert and die." << std::endl;


    SetType s(MAX_SET_SIZE);

   #ifdef ATK_DEBUG
    // NOTE: ATK_ASSSERT is disabled in release mode, so this test will only fail in debug mode
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    ASSERT_DEATH( s[MAX_SET_SIZE], "");
   #else
    std::cout <<"Did not check for assertion failure since assertions are compiled out in release mode." << std::endl;
   #endif
 */
   }

