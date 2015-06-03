/*
 * testSet.cxx
 *
 *  Created on: Apr 23, 2015
 *      Author: weiss27
 */

#include <iterator>
#include "gtest/gtest.h"

#include "meshapi/Utilities.hpp"
#include "meshapi/OrderedSet.hpp"


typedef asctoolkit::meshapi::OrderedSet SetType;
typedef SetType::iterator SetIterator;
static SetType::SizeType const MAX_SET_SIZE = 10;


TEST(gtest_meshapi_ordered_set,construct_ordered_set)
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
    for(SetType::SetPosition pos = SetType::SetPosition(); pos < static_cast<SetType::SetPosition>(s.size()); ++pos)
    {
        SetType::SetIndex idx = static_cast<SetType::SetIndex>(pos);
        EXPECT_EQ(idx,s[pos])
                <<"Random access iterator dereference to equal its position in the set";
        std::cout << "\t" << s[pos] <<"\n";
    }

    std::cout<<"\n --Using checked random access -- at()" << std::endl;
    for(SetType::SetPosition pos = SetType::SetPosition(); pos < static_cast<SetType::SetPosition>(s.size()); ++pos)
    {
        SetType::SetIndex idx = static_cast<SetType::SetIndex>(pos);
        EXPECT_EQ(idx,s.at(pos))
                <<"Expected checked random access iterator dereference to equal its position in the set";

        std::cout << "\t" << s.at(pos) <<"\n";
    }


    std::cout<<"\n --Using checked random access -- at() with invalid address" << std::endl;

    // add this line to avoid a warning in the output about thread safety
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    ASSERT_DEATH(s.at(MAX_SET_SIZE),"") << "tried to access out of range element";


    std::cout << "--\ndone." << std::endl;

}

TEST(gtest_meshapi_ordered_set,test_ordered_set_out_of_bounds)
{
    std::cout<<"\n****** Testing out of bounds access on initialized set-- code is expected to assert and die." << std::endl;

    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    SetType s(MAX_SET_SIZE);

    ASSERT_DEATH( s[MAX_SET_SIZE], "");
}


