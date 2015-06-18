/**
 * \file Set.cpp
 *
 *  \brief Implementation of the IndirectionSet class
 */

#include <stdexcept>
#include <sstream>

#include "IndirectionSet.hpp"

namespace asctoolkit {
namespace meshapi {

const NullSet IndirectionSet::s_nullSet;

bool IndirectionSet::isValid(bool verboseOutput) const
{
    bool bValid = true;

    std::stringstream errStr;

    if(verboseOutput)
    {
        errStr << "hello";
    }



    if(verboseOutput)
    {
        if( !bValid)
            std::cout<<" There was a problem: " << errStr.str() << std::endl;
    }

    return bValid;
}

} /* namespace meshapi */
} /* namespace asctoolkit */

