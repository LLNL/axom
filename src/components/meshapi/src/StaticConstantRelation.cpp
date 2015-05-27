/**
 * \file StaticConstantRelation.cpp
 *
 *  Created on: Apr 29, 2015
 *      Author: weiss27
 */

#include "StaticConstantRelation.hpp"

#include <sstream>
#include <iostream>
#include <iterator>

namespace asctoolkit {
namespace meshapi {

    StaticConstantRelation::StaticConstantRelation (Set* fromSet, Set* toSet)
    : m_stride(SetIndex()), m_fromSet(fromSet), m_toSet(toSet)
{
}

void StaticConstantRelation::setRelation(RelationVec const& toOffsets, SetIndex stride)
{
    m_stride = stride;

    m_toSetIndicesVec.clear();
    m_toSetIndicesVec.reserve(toOffsets.size());
    std::copy(toOffsets.begin(), toOffsets.end(), std::back_inserter(m_toSetIndicesVec));
}

bool StaticConstantRelation::isValid(bool verboseOutput) const
{
    bool bValid = true;

    std::stringstream sstr;

    if( m_fromSet == NULL || m_toSet == NULL)
    {
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

        // Check that the toSetIndices vector has the right size
        if( m_toSetIndicesVec.size() != (m_stride * m_fromSet->size()) )
        {
            if(verboseOutput)
            {
                sstr << "\n\t* toSetIndices has the wrong size."
                     << "\n\t-- from set size is: " << m_fromSet->size()
                     << "\n\t-- constant stride is: " << m_stride
                     << "\n\t-- expected relation size: " << (m_stride * m_fromSet->size())
                     << "\n\t-- actual size: " << m_toSetIndicesVec.size()
                     ;
            }
            bValid = false;
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
        if(bValid)
        {
            std::cout<<"(static,constant) Relation with stride " << m_stride << " was valid."<< std::endl;
        }
        else
        {
            std::cout<<"Relation was NOT valid.\n"
                     << sstr.str()
                     << std::endl;
        }

        std::cout<<"\n*** Detailed results of isValid on the relation.\n";
        if(m_fromSet)
            std::cout<< "\n** fromSet has size " << m_fromSet->size() <<": ";
        if(m_toSet)
            std::cout<< "\n** toSet has size " << m_toSet->size() <<": ";

        std::cout<< "\n** toSetIndices vec w/ size " << m_toSetIndicesVec.size() <<": ";
        std::copy(m_toSetIndicesVec.begin(), m_toSetIndicesVec.end(), std::ostream_iterator<SetIndex>(std::cout, " "));

    }

    return bValid;
}

} // namespace meshapi
} // namespace asctoolkit
