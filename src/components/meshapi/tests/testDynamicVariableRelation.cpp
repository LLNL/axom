/**
 * \file testDynamicVariableRelation.cxx
 *
 *  Created on: Apr 29, 2015
 *      Author: weiss27
 */


#include <iostream>
#include <iterator>

#include "gtest/gtest.h"


#include "meshapi/OrderedSet.hpp"
#include "meshapi/Relation.hpp"
#include "meshapi/DynamicVariableRelation.hpp"

using asctoolkit::meshapi::OrderedSet;
using asctoolkit::meshapi::DynamicVariableRelation;

typedef asctoolkit::meshapi::MeshIndexType IndexType;
typedef asctoolkit::meshapi::MeshSizeType SizeType;
// typedef asctoolkit::meshapi::Set::SetPosition PositionType;

const IndexType FROMSET_SIZE = 5;
const IndexType TOSET_SIZE = 8;


TEST(gtest_meshapi_dynamic_variable_relation,empty_relation)
{
    std::cout<<"\n****** Testing empty relation.  isValid() should be true." << std::endl;

    DynamicVariableRelation emptyRel;

    EXPECT_TRUE(emptyRel.isValid(true)) << "Empty relation was not valid";

    std::cout<<"\n****** done."<<std::endl;
}

TEST(gtest_meshapi_dynamic_variable_relation,test_uninitialized_relation)
{
    std::cout<<"\n****** Testing uninitialized relation.  isValid() should be false." << std::endl;

    OrderedSet fromSet(FROMSET_SIZE);
    OrderedSet toSet(TOSET_SIZE);

    DynamicVariableRelation emptyRel(&fromSet, &toSet);

    EXPECT_TRUE(emptyRel.isValid(true)) << "Empty relation was not initialized";

    std::cout<<"\n****** done."<<std::endl;
}

template<typename StrType, typename VecType>
void printVector(StrType const& msg, VecType const& vec)
{
    std::cout<< "\n** " << msg << "\n\t";
    std::cout<<"Array of size " << vec.size() <<": ";
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<IndexType>(std::cout, " "));
}

void generateIncrementingRelations(DynamicVariableRelation* rel)
{
    IndexType curIdx = IndexType();
    for(IndexType i=0; i < FROMSET_SIZE; ++i)
    {
        for(IndexType j=0; j <= i; ++j)
        {
            rel->insert(i, j % TOSET_SIZE );
            ++curIdx;
        }
    }
}

TEST(gtest_meshapi_dynamic_variable_relation,simple_relation)
{
    std::cout<<"\n****** Testing simple incrementing relation.  isValid() should be true." << std::endl;

    OrderedSet fromSet(FROMSET_SIZE);
    OrderedSet toSet(TOSET_SIZE);

    DynamicVariableRelation incrementingRel(&fromSet, &toSet);
    generateIncrementingRelations(&incrementingRel);

    for(IndexType idx=0; idx< static_cast<IndexType>(fromSet.size()); ++idx)
    {
        std::stringstream sstr;
        sstr << "Related to index " << idx;
        printVector(sstr.str(), incrementingRel[idx]);
    }

    EXPECT_TRUE(incrementingRel.isValid(true)) << "Incrementing relation was not valid";

    typedef OrderedSet::iterator SetIter;
    typedef DynamicVariableRelation::RelationVecConstIterator RelSetConstIter;

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

TEST(gtest_meshapi_dynamic_variable_relation,test_iterator_range)
{
    std::cout<<"\n****** Testing range function on incrementing relation." << std::endl;

    OrderedSet fromSet(FROMSET_SIZE);
    OrderedSet toSet(TOSET_SIZE);

    DynamicVariableRelation incrementingRel(&fromSet, &toSet);
    generateIncrementingRelations(&incrementingRel);

    for(IndexType idx=0; idx< static_cast<IndexType>(fromSet.size()); ++idx)
    {
        std::stringstream sstr;
        sstr << "Related to index " << idx;
        printVector(sstr.str(), incrementingRel[idx]);
    }

    EXPECT_TRUE(incrementingRel.isValid(true)) << "Incrementing relation was not valid";

    typedef OrderedSet::iterator SetIter;
    typedef DynamicVariableRelation::RelationVecConstIterator RelSetConstIter;
    typedef DynamicVariableRelation::RelationVecConstIteratorPair RelSetConstIterPair;

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


TEST(gtest_meshapi_dynamic_variable_relation,double_subscript_test)
{
    std::cout<<"\n****** Testing access via double subscript." << std::endl;

    OrderedSet fromSet(FROMSET_SIZE);
    OrderedSet toSet(TOSET_SIZE);

    DynamicVariableRelation incrementingRel(&fromSet, &toSet);
    generateIncrementingRelations(&incrementingRel);

    for(IndexType idx=0; idx< static_cast<IndexType>(fromSet.size()); ++idx)
    {
        std::stringstream sstr;
        sstr << "Related to index " << idx;
        printVector(sstr.str(), incrementingRel[idx]);
    }

    EXPECT_TRUE(incrementingRel.isValid(true)) << "Incrementing relation was not valid";

    typedef OrderedSet::iterator SetIter;
    typedef DynamicVariableRelation::RelationVecConstIterator RelSetConstIter;

    std::cout<<"\n\tLooking at relation's stored values...";
    for(SetIter sIt = fromSet.begin(), sItEnd = fromSet.end(); sIt != sItEnd; ++sIt)
    {
        std::cout<<"\n\tInspecting element " << *sIt << " of first set.";

        for(IndexType idx=0; idx< static_cast<IndexType>(incrementingRel.size(*sIt)); ++idx)
        {
            EXPECT_EQ( idx, incrementingRel[*sIt][idx]) << "incrementing relation's value was incorrect";
        }
    }

    std::cout<<"\n****** done."<<std::endl;
}




