/*
 * testRelation.cxx
 *
 *  Created on: Apr 29, 2015
 *      Author: weiss27
 */


#include <iostream>
#include <iterator>

#include "gtest/gtest.h"


#include "meshapi/OrderedSet.hpp"
#include "meshapi/Relation.hpp"

using asctoolkit::meshapi::OrderedSet;
using asctoolkit::meshapi::Relation;

typedef asctoolkit::meshapi::MeshIndexType IndexType;
const IndexType FROMSET_SIZE = 5;
const IndexType TOSET_SIZE = 8;


TEST(gtest_meshapi_relation,empty_relation)
{
    Relation emptyRel(NULL, NULL);

    EXPECT_TRUE(emptyRel.isValid()) << "Empty relation was not valid";



}

TEST(gtest_meshapi_relation,test_uninitialized_relation)
{
    OrderedSet fromSet(FROMSET_SIZE);
    OrderedSet toSet(TOSET_SIZE);

    Relation emptyRel(&fromSet, &toSet);

    EXPECT_FALSE(emptyRel.isValid(true)) << "Empty relation was not initialized";



}

template<typename StrType, typename VecType>
void printVector(StrType const& msg, VecType const& vec)
{
    std::cout<< "\n** " << msg << "\n\t";
    std::cout<<"Array of size " << vec.size() <<": ";
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<IndexType>(std::cout, " "));
}

template<typename VecType>
void generateRelations(VecType* begins, VecType* offsets)
{
    VecType& beginsVec = *begins;
    VecType& offsetsVec = *offsets;

    IndexType curIdx = IndexType();
    for(IndexType i=0; i < FROMSET_SIZE; ++i)
    {
        beginsVec[i] = curIdx;
        for(IndexType j=0; j <= i; ++j)
        {
            offsetsVec.push_back( j % TOSET_SIZE );
            ++curIdx;
        }
    }
    beginsVec[FROMSET_SIZE] = curIdx;
}

TEST(gtest_meshapi_relation,simple_relation)
{
    OrderedSet fromSet(FROMSET_SIZE);
    OrderedSet toSet(TOSET_SIZE);

    Relation incrementingRel(&fromSet, &toSet);

    typedef Relation::RelationVec IndexVec;
    IndexVec begins(FROMSET_SIZE +1);
    IndexVec offsets;

    printVector("begins vector", begins);
    printVector("offsets vector", offsets);

    generateRelations(&begins, &offsets);

    printVector("begins vector", begins);
    printVector("offsets vector", offsets);


    incrementingRel.setRelation(begins, offsets);


    EXPECT_TRUE(incrementingRel.isValid(true)) << "Incrementing relation was not valid";

}


