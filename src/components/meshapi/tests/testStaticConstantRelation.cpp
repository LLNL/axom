/**
 * \file testStaticConstantRelation.cxx
 *
 *  Created on: Apr 29, 2015
 *      Author: weiss27
 */


#include <iostream>
#include <iterator>

#include "gtest/gtest.h"


#include "meshapi/OrderedSet.hpp"
#include "meshapi/Relation.hpp"
#include "meshapi/StaticConstantRelation.hpp"

using asctoolkit::meshapi::OrderedSet;
using asctoolkit::meshapi::StaticConstantRelation;

typedef asctoolkit::meshapi::MeshIndexType IndexType;
typedef asctoolkit::meshapi::MeshSizeType SizeType;

const IndexType FROMSET_SIZE = 5;
const IndexType TOSET_SIZE = 6;


TEST(gtest_meshapi_static_constant_relation,empty_relation)
{
    std::cout<<"\n****** Testing empty relation.  isValid() should be true." << std::endl;

    StaticConstantRelation emptyRel;

    EXPECT_TRUE(emptyRel.isValid(true)) << "Empty relation was not valid";

    std::cout<<"\n****** done."<<std::endl;
}

TEST(gtest_meshapi_static_constant_relation,empty_relation_out_of_bounds)
{
    std::cout<<"\n****** Testing access on empty relation -- code is expected to assert and die." << std::endl;
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    StaticConstantRelation emptyRel;

    ASSERT_DEATH( emptyRel[FROMSET_SIZE], "");
}

TEST(gtest_meshapi_static_constant_relation,test_uninitialized_relation)
{
    std::cout<<"\n****** Testing uninitialized relation.  isValid() should be TRUE since stride is 0 by default." << std::endl;

    OrderedSet fromSet(FROMSET_SIZE);
    OrderedSet toSet(TOSET_SIZE);

    StaticConstantRelation emptyRel(&fromSet, &toSet);

    EXPECT_TRUE(emptyRel.isValid(true)) << "Constant relation with stride 0 should be valid";

    std::cout<<"\n****** done."<<std::endl;
}

template<typename StrType, typename VecType>
void printVector(StrType const& msg, VecType const& vec)
{
    std::cout<< "\n** " << msg << "\n\t";
    std::cout<<"Array of size " << vec.size() <<": ";
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<IndexType>(std::cout, " "));
}

template<typename VecType>
void generateIncrementingRelations(IndexType stride, VecType* offsets)
{
    VecType& offsetsVec = *offsets;

    IndexType curIdx = IndexType();
    for(IndexType i=0; i < FROMSET_SIZE; ++i)
    {
        for(IndexType j=0; j < stride; ++j)
        {
            offsetsVec.push_back( (i+j) % TOSET_SIZE );
            ++curIdx;
        }
    }
}

TEST(gtest_meshapi_static_constant_relation,simple_relation)
{
    std::cout<<"\n****** Testing simple incrementing relation.  isValid() should be true." << std::endl;

    OrderedSet fromSet(FROMSET_SIZE);
    OrderedSet toSet(TOSET_SIZE);

    StaticConstantRelation incrementingRel(&fromSet, &toSet);

    typedef StaticConstantRelation::RelationVec IndexVec;
    IndexVec offsets;

    printVector("offsets vector", offsets);

    IndexType const ELEM_STRIDE = 5;

    generateIncrementingRelations(ELEM_STRIDE, &offsets);

    printVector("offsets vector", offsets);


    incrementingRel.setRelation(offsets, ELEM_STRIDE);


    EXPECT_TRUE(incrementingRel.isValid(true)) << "Incrementing relation was not valid";

    typedef OrderedSet::iterator SetIter;
    typedef StaticConstantRelation::RelationVecConstIterator RelSetConstIter;

    std::cout<<"\n\tLooking at relation's stored values...";
    for(SetIter sIt = fromSet.begin(), sItEnd = fromSet.end(); sIt != sItEnd; ++sIt)
    {
        std::cout<<"\n\tInspecting element " << *sIt << " of first set.";

        SizeType actualSize = incrementingRel.size( *sIt);
        SizeType expectedSize = ELEM_STRIDE;

        std::cout <<"\n\t\tExpected: " << expectedSize;
        std::cout <<"\n\t\tActual: " <<  actualSize <<"\n";

        EXPECT_EQ( expectedSize, actualSize ) << "relation for this element was incorrect size.";

        IndexType fromSetEltNum = std::distance(fromSet.begin(), sIt);

        RelSetConstIter toSetBegin = incrementingRel.begin(*sIt);
        RelSetConstIter toSetEnd = incrementingRel.end(*sIt);
        for(RelSetConstIter innerIt = toSetBegin; innerIt != toSetEnd; ++innerIt)
        {
            IndexType toSetEltNum = std::distance(toSetBegin, innerIt);

            std::cout <<"\n\t\t " << toSetEltNum <<": " << *innerIt ;

            IndexType expectedVal =  (fromSetEltNum + toSetEltNum) % TOSET_SIZE;
            IndexType actualVal = *innerIt;
            ASSERT_EQ( expectedVal, actualVal) << "incrementing relation's value was incorrect";
        }
    }
    std::cout<<"\n****** done."<<std::endl;
}

TEST(gtest_meshapi_static_constant_relation,initialized_rel_out_of_bounds)
{
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    std::cout<<"\n****** Testing out of bounds access on initialized relation.  Code is expected to assert and die." << std::endl;

    OrderedSet fromSet(FROMSET_SIZE);
    OrderedSet toSet(TOSET_SIZE);
    StaticConstantRelation incrementingRel(&fromSet, &toSet);

    typedef StaticConstantRelation::RelationVec IndexVec;
    IndexVec offsets;

    IndexType const ELEM_STRIDE = 5;
    generateIncrementingRelations(ELEM_STRIDE, &offsets);
    incrementingRel.setRelation(offsets, ELEM_STRIDE);

    ASSERT_DEATH( incrementingRel[FROMSET_SIZE], "");
}


TEST(gtest_meshapi_static_constant_relation,test_iterator_range)
{
    std::cout<<"\n****** Testing range function on incrementing relation." << std::endl;

    OrderedSet fromSet(FROMSET_SIZE);
    OrderedSet toSet(TOSET_SIZE);

    StaticConstantRelation incrementingRel(&fromSet, &toSet);

    typedef StaticConstantRelation::RelationVec IndexVec;
    IndexVec offsets;

    IndexType const ELEM_STRIDE = 5;
    generateIncrementingRelations(ELEM_STRIDE, &offsets);
    incrementingRel.setRelation(offsets, ELEM_STRIDE);

    EXPECT_TRUE(incrementingRel.isValid(true)) << "Incrementing relation was not valid";

    typedef OrderedSet::iterator SetIter;
    typedef StaticConstantRelation::RelationVecConstIterator RelSetConstIter;
    typedef StaticConstantRelation::RelationVecConstIteratorPair RelSetConstIterPair;

    std::cout<<"\n\tLooking at relation's stored values...";
    for(SetIter sIt = fromSet.begin(), sItEnd = fromSet.end(); sIt != sItEnd; ++sIt)
    {
        std::cout<<"\n\tInspecting element " << *sIt << " of first set.";
        IndexType fromSetEltNum = std::distance(fromSet.begin(), sIt);

        RelSetConstIterPair toSetItPair = incrementingRel.range(*sIt);
        for(RelSetConstIter it = toSetItPair.first; it < toSetItPair.second; ++it)
        {
            IndexType toSetEltNum = std::distance(toSetItPair.first, it);

            std::cout <<"\n\t\t " << toSetEltNum <<": " << *it ;

            IndexType expectedVal =  (fromSetEltNum + toSetEltNum) % TOSET_SIZE;
            ASSERT_EQ( expectedVal, *it) << "incrementing relation's value was incorrect";
        }
    }

    std::cout<<"\n****** done."<<std::endl;
}


TEST(gtest_meshapi_static_constant_relation,double_subscript_test)
{
    std::cout<<"\n****** Testing access via double subscript." << std::endl;

    OrderedSet fromSet(FROMSET_SIZE);
    OrderedSet toSet(TOSET_SIZE);

    StaticConstantRelation incrementingRel(&fromSet, &toSet);

    typedef StaticConstantRelation::RelationVec IndexVec;
    IndexVec offsets;

    IndexType const ELEM_STRIDE = 5;
    generateIncrementingRelations(ELEM_STRIDE, &offsets);
    incrementingRel.setRelation(offsets, ELEM_STRIDE);

    EXPECT_TRUE(incrementingRel.isValid(true)) << "Incrementing relation was not valid";

    typedef OrderedSet::iterator SetIter;
    typedef StaticConstantRelation::RelationVecConstIterator RelSetConstIter;

    std::cout<<"\n\tLooking at relation's stored values...";
    for(SetIter sIt = fromSet.begin(), sItEnd = fromSet.end(); sIt != sItEnd; ++sIt)
    {
        std::cout<<"\n\tInspecting element " << *sIt << " of first set.";

        for(IndexType idx=0; idx< incrementingRel.size(*sIt); ++idx)
        {
            IndexType expectedVal =  (*sIt + idx) % TOSET_SIZE;
            IndexType actualVal = incrementingRel[*sIt][idx];
            EXPECT_EQ( expectedVal, actualVal) << "incrementing relation's value was incorrect";
        }
    }

    std::cout<<"\n****** done."<<std::endl;
}



