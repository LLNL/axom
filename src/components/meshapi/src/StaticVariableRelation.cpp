/*
 * StaticVariableRelation.cpp
 *
 *  Created on: Apr 29, 2015
 *      Author: weiss27
 */

#include "StaticVariableRelation.hpp"

#include <sstream>
#include <iostream>
#include <iterator>

namespace asctoolkit {
namespace meshapi {

    StaticVariableRelation::StaticVariableRelation (Set* fromSet, Set* toSet)
    : m_fromSet(fromSet), m_toSet(toSet)
{
}

void StaticVariableRelation::setRelation(RelationVec const& beginsVec, RelationVec const& toOffsets)
{
    m_fromSetBeginsVec.clear();
    m_fromSetBeginsVec.reserve(beginsVec.size());
    std::copy(beginsVec.begin(), beginsVec.end(), std::back_inserter(m_fromSetBeginsVec));

    m_toSetIndicesVec.clear();
    m_toSetIndicesVec.reserve(toOffsets.size());
    std::copy(toOffsets.begin(), toOffsets.end(), std::back_inserter(m_toSetIndicesVec));
}

bool StaticVariableRelation::isValid(bool verboseOutput) const
{
    bool bValid = true;

    std::stringstream sstr;


    if( m_fromSet == NULL || m_toSet == NULL)
    {
        if(!m_fromSetBeginsVec.empty())
        {
            if(verboseOutput)
            {
                sstr << "\n\t* fromSetBeginsVec was not empty "
                    <<" -- fromSet was " << (m_fromSet == NULL ? "" : " not ") << "null"
                    <<" , toSet was " << (m_toSet == NULL ? "" : " not ") << "null";
            }

            bValid = false;
        }
        if(!m_toSetIndicesVec.empty())
        {
            if(verboseOutput)
            {
                sstr << "\n\t* toSetIndicesVec was not empty "
                    <<" -- fromSet was " << (m_fromSet == NULL ? "" : " not ") << "null"
                    <<" , toSet was " << (m_toSet == NULL ? "" : " not ") << "null";
            }

            bValid = false;
        }
    }
    else
    {
        if(verboseOutput)
            sstr << "\n\t* Neither set was null";

        // Check that the fromSet vector has the correct size
        if( m_fromSetBeginsVec.size() != (m_fromSet->size() +1) )
        {
            if(verboseOutput)
            {
                sstr << "\n\t* fromSetBeginsVec was the wrong size."
                 << " -- expected " << (m_fromSet->size() +1) << ", actual " << m_fromSetBeginsVec.size() << ")";
            }

            bValid = false;
        }

        // Check that no element of the fromSetBegins vector points outside of the toSetIndices vector
        for(RelationVecConstIterator it = m_fromSetBeginsVec.begin(), itEnd = m_fromSetBeginsVec.end(); it != itEnd; ++it)
        {
            if( *it > m_toSetIndicesVec.size() )
            {
                if(verboseOutput)
                {
                    sstr << "\n\t* fromSetBeginsVec had an out-of-range element."
                         << " -- value of element " << std::distance(m_fromSetBeginsVec.begin(), it) << " was " << *it
                         << ". Max possible value should be " << m_toSetIndicesVec.size() <<"." ;
                }

                bValid = false;
            }
        }

        // Check that all elements of the toSetIndices vector point to valid set elements
        for(RelationVecConstIterator it = m_toSetIndicesVec.begin(), itEnd = m_toSetIndicesVec.end(); it != itEnd; ++it)
        {
            if( *it >= m_toSet->size() )
            {
                if(verboseOutput)
                {
                    sstr << "\n\t* toSetIndices had an out-of-range element."
                         << " -- value of element " << std::distance(m_toSetIndicesVec.begin(), it) << " was " << *it
                         << ". Max possible value should be " << m_toSet->size() <<"." ;
                }

                bValid = false;
            }
        }
    }


    if(verboseOutput)
    {
        std::cout<<"\n*** Detailed results of isValid on the relation.\n";
        if(bValid)
        {
            std::cout<<"Relation was valid."<< std::endl;
        }
        else
        {
            std::cout<<"Relation was NOT valid.\n"
                     << sstr.str()
                     << std::endl;
        }

        std::cout<< "\n** fromSetBeginsVec vec w/ size " << m_fromSetBeginsVec.size() <<": ";
        std::copy(m_fromSetBeginsVec.begin(), m_fromSetBeginsVec.end(), std::ostream_iterator<SetIndex>(std::cout, " "));

        std::cout<< "\n** toSetIndices vec w/ size " << m_toSetIndicesVec.size() <<": ";
        std::copy(m_toSetIndicesVec.begin(), m_toSetIndicesVec.end(), std::ostream_iterator<SetIndex>(std::cout, " "));

    }

    return bValid;
}

} // namespace meshapi
} // namespace asctoolkit
