/*
 * OrderedSet.cpp
 *
 *  Created on: Apr 23, 2015
 *      Author: weiss27
 */

#include <stdexcept>
#include <sstream>
#include <iostream>
#include "RangeSet.hpp"

namespace asctoolkit {
namespace meshapi {


const NullSet RangeSet::s_nullSet;

bool RangeSet::isValid(bool verboseOutput) const
{
    bool bValid = true;

    std::stringstream errStr;

    // Ensure 0 <= m_lowerIdx <= m_upperIdx
#if 0
    // Compiled out when SizeType is unsigned..
    if( m_lowerIdx < SizeType() )
    {
        bValid = false;

        if(verboseOutput)
            errStr <<"Lower index ("<< m_lowerIdx << ") was less than " << SizeType() <<"\n";
    }
#endif

    if( m_lowerIdx > m_upperIdx )
    {
        bValid = false;

        if(verboseOutput)
            errStr <<"Lower index ("<< m_lowerIdx << ") must be less than or equal to upper index (" << m_upperIdx <<")";
    }

    if(verboseOutput)
    {
        std::stringstream sstr;

        sstr<<"\n*** Detailed results of isValid on the RangeSet.\n";
        if(bValid)
        {
            sstr<<"Set was valid."<< std::endl;
        }
        else
        {
            sstr<<"Set was NOT valid.\n"
                     << errStr.str()
                     << std::endl;
        }

        sstr<<"\n** RangeSet [" << m_lowerIdx <<"," << m_upperIdx << ") has " << size() <<" elements.";
        std::cout<< sstr.str() << std::endl;
    }

    return bValid;
}

} /* namespace meshapi */
} /* namespace asctoolkit */

