/**
 * \file testStaticVariableRelation.cxx
 *
 *  Created on: Apr 29, 2015
 *      Author: weiss27
 */


#include <iostream>
#include <iterator>

#include "gtest/gtest.h"


#include "meshapi/RangeSet.hpp"
#include "meshapi/Relation.hpp"
#include "meshapi/StaticVariableRelation.hpp"

using asctoolkit::meshapi::RangeSet;
using asctoolkit::meshapi::StaticVariableRelation;

typedef asctoolkit::meshapi::MeshIndexType IndexType;
typedef asctoolkit::meshapi::MeshSizeType SizeType;

const IndexType FROMSET_SIZE = 5;
const IndexType TOSET_SIZE = 8;


TEST(gtest_meshapi_static_variable_relation,empty_relation)
{
    std::cout<<"\n****** Testing empty relation.  isValid() should be true." << std::endl;

    StaticVariableRelation emptyRel;

    EXPECT_TRUE(emptyRel.isValid(true)) << "Empty relation was not valid";



    std::cout<<"\n****** done."<<std::endl;
}

TEST(gtest_meshapi_static_variable_relation,empty_relation_out_of_bounds)
{
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    StaticVariableRelation emptyRel;

    ASSERT_DEATH( emptyRel[FROMSET_SIZE], "");
}



TEST(gtest_meshapi_static_variable_relation,test_uninitialized_relation)
{
    std::cout<<"\n****** Testing uninitialized relation.  isValid() should be false." << std::endl;

    RangeSet fromSet(FROMSET_SIZE);
    RangeSet toSet(TOSET_SIZE);

    StaticVariableRelation emptyRel(&fromSet, &toSet);

    EXPECT_FALSE(emptyRel.isValid(true)) << "Empty relation was not initialized";

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
void generateIncrementingRelations(VecType* begins, VecType* offsets)
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

TEST(gtest_meshapi_static_variable_relation,simple_relation)
{
    std::cout<<"\n****** Testing simple incrementing relation.  isValid() should be true." << std::endl;

    RangeSet fromSet(FROMSET_SIZE);
    RangeSet toSet(TOSET_SIZE);

    StaticVariableRelation incrementingRel(&fromSet, &toSet);

    typedef StaticVariableRelation::RelationVec IndexVec;
    IndexVec begins(FROMSET_SIZE +1);
    IndexVec offsets;

    printVector("begins vector", begins);
    printVector("offsets vector", offsets);

    generateIncrementingRelations(&begins, &offsets);

    printVector("begins vector", begins);
    printVector("offsets vector", offsets);


    incrementingRel.setRelation(begins, offsets);


    EXPECT_TRUE(incrementingRel.isValid(true)) << "Incrementing relation was not valid";

    typedef RangeSet::iterator SetIter;
    typedef StaticVariableRelation::RelationVecConstIterator RelSetConstIter;

    std::cout<<"\n\tLooking at relation's stored values...";
    for(SetIter sIt = fromSet.begin(), sItEnd = fromSet.end(); sIt != sItEnd; ++sIt)
    {
        std::cout<<"\n\tInspecting element " << *sIt << " of first set.";

        SizeType actualSize = incrementingRel.size( *sIt);
        SizeType expectedSize = std::distance(fromSet.begin(), sIt) +1;

        std::cout <<"\n\t\tExpected: " << expectedSize;
        std::cout <<"\n\t\tActual: " <<  actualSize <<"\n";

        EXPECT_EQ( expectedSize, actualSize ) << "relation for this element was incorrect size.";

        RelSetConstIter toSetBegin = incrementingRel.begin(*sIt);
        RelSetConstIter toSetEnd = incrementingRel.end(*sIt);
        for(RelSetConstIter innerIt = toSetBegin; innerIt != toSetEnd; ++innerIt)
        {
            IndexType eltNum = std::distance(toSetBegin, innerIt);

            std::cout <<"\n\t\t " << eltNum <<": " << *innerIt ;

            IndexType expectedVal =  (eltNum ) % TOSET_SIZE;
            IndexType actualVal = *innerIt;
            ASSERT_EQ( expectedVal, actualVal) << "incrementing relation's value was incorrect";
        }
    }

    std::cout<<"\n****** done."<<std::endl;
}

TEST(gtest_meshapi_static_variable_relation,initialized_rel_out_of_bounds)
{
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    std::cout<<"\n****** Testing simple incrementing relation.  isValid() should be true." << std::endl;

    RangeSet fromSet(FROMSET_SIZE);
    RangeSet toSet(TOSET_SIZE);
    StaticVariableRelation incrementingRel(&fromSet, &toSet);

    typedef StaticVariableRelation::RelationVec IndexVec;
    IndexVec begins(FROMSET_SIZE +1);
    IndexVec offsets;

    generateIncrementingRelations(&begins, &offsets);
    incrementingRel.setRelation(begins, offsets);

    ASSERT_DEATH( incrementingRel[FROMSET_SIZE], "");
}

TEST(gtest_meshapi_static_variable_relation,test_iterator_range)
{
    std::cout<<"\n****** Testing range function on incrementing relation." << std::endl;

    RangeSet fromSet(FROMSET_SIZE);
    RangeSet toSet(TOSET_SIZE);

    StaticVariableRelation incrementingRel(&fromSet, &toSet);

    typedef StaticVariableRelation::RelationVec IndexVec;
    IndexVec begins(FROMSET_SIZE +1);
    IndexVec offsets;
    generateIncrementingRelations(&begins, &offsets);
    incrementingRel.setRelation(begins, offsets);

    EXPECT_TRUE(incrementingRel.isValid(true)) << "Incrementing relation was not valid";

    typedef RangeSet::iterator SetIter;
    typedef StaticVariableRelation::RelationVecConstIterator RelSetConstIter;
    typedef StaticVariableRelation::RelationVecConstIteratorPair RelSetConstIterPair;

    std::cout<<"\n\tLooking at relation's stored values...";
    for(SetIter sIt = fromSet.begin(), sItEnd = fromSet.end(); sIt != sItEnd; ++sIt)
    {
        std::cout<<"\n\tInspecting element " << *sIt << " of first set.";

        RelSetConstIterPair toSetItPair = incrementingRel.range(*sIt);
        for(RelSetConstIter it = toSetItPair.first; it < toSetItPair.second; ++it)
        {
            IndexType eltNum = std::distance(toSetItPair.first, it);

            std::cout <<"\n\t\t " << eltNum <<": " << *it ;

            IndexType expectedVal =  (eltNum ) % TOSET_SIZE;
            IndexType actualVal = *it;
            ASSERT_EQ( expectedVal, actualVal) << "incrementing relation's value was incorrect";
        }
    }

    std::cout<<"\n****** done."<<std::endl;
}


TEST(gtest_meshapi_static_variable_relation,double_subscript_test)
{
    std::cout<<"\n****** Testing access via double subscript." << std::endl;

    RangeSet fromSet(FROMSET_SIZE);
    RangeSet toSet(TOSET_SIZE);

    StaticVariableRelation incrementingRel(&fromSet, &toSet);

    typedef StaticVariableRelation::RelationVec IndexVec;
    IndexVec begins(FROMSET_SIZE +1);
    IndexVec offsets;
    generateIncrementingRelations(&begins, &offsets);
    incrementingRel.setRelation(begins, offsets);

    EXPECT_TRUE(incrementingRel.isValid(true)) << "Incrementing relation was not valid";

    typedef RangeSet::iterator SetIter;
    typedef StaticVariableRelation::RelationVecConstIterator RelSetConstIter;

    std::cout<<"\n\tLooking at relation's stored values...";
    for(SetIter sIt = fromSet.begin(), sItEnd = fromSet.end(); sIt != sItEnd; ++sIt)
    {
        std::cout<<"\n\tInspecting element " << *sIt << " of first set.";

        for(IndexType idx=0; idx< static_cast<IndexType>(incrementingRel.size(*sIt ) ); ++idx)
        {
            IndexType actualVal = incrementingRel[*sIt][idx];
            IndexType expectedVal = idx;
            EXPECT_EQ( expectedVal, actualVal) << "incrementing relation's value was incorrect";
        }
    }

    std::cout<<"\n****** done."<<std::endl;
}
