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
static SetType::size_type const MAX_SET_SIZE = 10;


TEST(gtest_meshapi_ordered_set,construct_ordered_set)
{

    SetType s(MAX_SET_SIZE);

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
    for(SetType::size_type i = SetType::size_type(); i < s.size(); ++i)
    {
        EXPECT_EQ(i,s[i])
                <<"Random access iterator dereference to equal its position in the set";
        std::cout << "\t" << s[i] <<"\n";
    }

    std::cout<<"\n --Using checked random access -- at()" << std::endl;
    for(SetType::size_type i = SetType::size_type(); i < s.size(); ++i)
    {
        EXPECT_EQ(i,s.at(i))
                <<"Expected checked random access iterator dereference to equal its position in the set";

        std::cout << "\t" << s.at(i) <<"\n";
    }


    std::cout<<"\n --Using checked random access -- at() with invalid address" << std::endl;
    try{

        s.at(MAX_SET_SIZE);


        EXPECT_TRUE(false)  << "OrderedSet should have thrown an out_of_order exception when accessing element "
                            << MAX_SET_SIZE <<".";
    }
    catch(std::out_of_range oor)
    {
        EXPECT_TRUE(true);
    }
    catch(...)
    {
        EXPECT_TRUE(false) <<" Unexpected exception with message when accessing out of range element using checked access at()";
    }

    std::cout << "--\ndone." << std::endl;

}
