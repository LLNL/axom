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
#include "meshapi/Map.hpp"


typedef asctoolkit::meshapi::OrderedSet SetType;
typedef asctoolkit::meshapi::Map<int> IntMap;
typedef asctoolkit::meshapi::Map<double> RealMap;
typedef SetType::SetIndex  SetIndex;

typedef SetType::iterator SetIterator;
static SetType::SizeType const MAX_SET_SIZE = 10;

TEST(gtest_meshapi_map,construct_empty_map)
{
    IntMap m;

    EXPECT_TRUE(m.isValid(true));
}

template<typename T>
bool constructAndTestMap()
{
    SetType s(MAX_SET_SIZE);

    std::cout<<"\nCreating set of size " << s.size() << std::endl;
    EXPECT_EQ(s.size(), MAX_SET_SIZE);
    EXPECT_TRUE(s.isValid());

    std::cout<<"\nCreating " << asctoolkit::meshapi::util::TypeToString<T>::to_string() << " map on the set " << std::endl;
    asctoolkit::meshapi::Map<T> m(&s);
    EXPECT_TRUE(m.isValid());

    std::cout<<"\nSetting the elements.";
    double multFac = 1.0001;
    for(SetIndex idx=0; idx < static_cast<SetIndex>(m.size()); ++idx)
    {
        m[idx] = static_cast<T>(idx * multFac);
    }

    std::cout<<"\nChecking the elements.";
    for(SetIndex idx=0; idx < static_cast<SetIndex>(m.size()); ++idx)
    {
        EXPECT_EQ(m[idx], static_cast<T>(idx * multFac) );
    }

    EXPECT_TRUE(m.isValid(true));

    return true;
}

TEST(gtest_meshapi_map,construct_int_map)
{
    EXPECT_TRUE( constructAndTestMap<int>() );
}

TEST(gtest_meshapi_map,construct_double_map)
{
    EXPECT_TRUE( constructAndTestMap<double>());
}
