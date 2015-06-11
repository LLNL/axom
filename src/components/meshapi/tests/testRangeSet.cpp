/*
 * testSet.cxx
 *
 *  Created on: Apr 23, 2015
 *      Author: weiss27
 */

#include <iterator>
#include "gtest/gtest.h"

#include "meshapi/Utilities.hpp"
#include "meshapi/RangeSet.hpp"


typedef asctoolkit::meshapi::RangeSet SetType;
typedef SetType::iterator SetIterator;
typedef SetType::SetPosition SetPosition;
typedef SetType::SetElement SetElement;

static const SetPosition MAX_SET_SIZE = 10;


TEST(gtest_meshapi_range_set,construct_range_set)
{

    SetType s(MAX_SET_SIZE);

    EXPECT_TRUE(s.isValid());

    std::cout<<"Iterating through set of size " << s.size() << std::endl;
    EXPECT_EQ(s.size(), MAX_SET_SIZE);


    std::cout<<"\n --Using begin/end" << std::endl;
    for(SetIterator it=s.begin(), itEnd=s.end(); it != itEnd; ++it)
    {
        EXPECT_EQ( std::distance(s.begin(), it), *it )
                << "Iterator dereference should be equal to its position in the set";
        std::cout << "\t" << *it <<"\n";
    }

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


    std::cout<<"\n --Using checked random access -- at() with invalid address" << std::endl;

    // add this line to avoid a warning in the output about thread safety
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    ASSERT_DEATH(s.at(MAX_SET_SIZE),"") << "tried to access out of range element";


    std::cout << "--\ndone." << std::endl;

}

TEST(gtest_meshapi_range_set,test_range_set_out_of_bounds)
{
    std::cout<<"\n****** Testing out of bounds access on initialized set-- code is expected to assert and die." << std::endl;

    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    SetType s(MAX_SET_SIZE);

    ASSERT_DEATH( s[MAX_SET_SIZE], "");
}


