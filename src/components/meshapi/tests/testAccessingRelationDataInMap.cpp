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
#include "meshapi/Map.hpp"


using asctoolkit::meshapi::RangeSet;
using asctoolkit::meshapi::StaticVariableRelation;

typedef RangeSet::ElementType ElementType;
typedef RangeSet::PositionType SetPosition;

const SetPosition FROMSET_SIZE = 10;
const SetPosition TOSET_SIZE = 8;

template<typename StrType, typename VecType>
void printVector(StrType const& msg, VecType const& vec)
{
    std::cout<< "\n** " << msg << "\n\t";
    std::cout<<"Array of size " << vec.size() <<": ";
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<SetPosition>(std::cout, " "));
}

template<typename VecType>
void generateIncrementingRelations(VecType* begins, VecType* offsets)
{
    VecType& beginsVec = *begins;
    VecType& offsetsVec = *offsets;

    SetPosition curIdx = SetPosition();
    for(SetPosition i=0; i < FROMSET_SIZE; ++i)
    {
        beginsVec[i] = curIdx;
        for(SetPosition j=0; j <= i; ++j)
        {
            offsetsVec.push_back( j % TOSET_SIZE );
            ++curIdx;
        }
    }
    beginsVec[FROMSET_SIZE] = curIdx;
}

TEST(gtest_meshapi_set_relation_map,access_pattern)
{
    std::cout<<"\n****** Testing accessing relation data." << std::endl;

    RangeSet fromSet(FROMSET_SIZE);
    RangeSet toSet(TOSET_SIZE);

    StaticVariableRelation incrementingRel(&fromSet, &toSet);

    typedef StaticVariableRelation::RelationVec IndexVec;
    IndexVec begins(FROMSET_SIZE +1);
    IndexVec offsets;
    generateIncrementingRelations(&begins, &offsets);
    incrementingRel.bindRelationData(begins, offsets);


    // Note: Nothing requires the relations elements to be unique -- the relation can still be valid with duplicates
    EXPECT_TRUE(incrementingRel.isValid(true)) << "Incrementing relation was not valid";

    typedef RangeSet::iterator SetIter;
    typedef StaticVariableRelation::RelationVecConstIterator RelSetConstIter;

    std::cout<<"\n\tLooking at relation's stored values...";
    for(SetPosition fromPos = SetPosition(); fromPos < fromSet.size(); ++fromPos)
    {
        std::cout<<"\n\tInspecting element " << fromSet[fromPos]
                <<" in position "<< fromPos << " of first set.";

        for(SetPosition idx=0; idx< incrementingRel.size( fromPos ); ++idx)
        {
            SetPosition posInToSet_actual = incrementingRel[fromPos][idx];
            SetPosition posInToSet_expected = idx % toSet.size();
            EXPECT_EQ( posInToSet_expected, posInToSet_actual) << "incrementing relation's value was incorrect";

            std::cout<<"\n\t\t pos: " << idx
                    <<" ToSet position: " << incrementingRel[fromPos][idx]
                    << " ToSet element " << toSet[ incrementingRel[fromPos][idx] ]
                    ;
        }
    }
    std::cout<<"\n****** done."<<std::endl;
}
