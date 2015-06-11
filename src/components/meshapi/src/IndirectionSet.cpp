/*
 * Set.cpp
 *
 *  Created on: Apr 23, 2015
 *      Author: weiss27
 */

#include <stdexcept>
#include <sstream>

#include "meshapi/Set.hpp"

namespace asctoolkit {
namespace meshapi {

bool IndirectionSet::isValid(bool verboseOutput)
{
    bool valid;

    std::stringstream errStr;

    if(verboseOutput)
    {
        errStr << "hello";
    }



    if(verboseOutput)
    {
        if( !valid)
            std::cout<<" There was a problem: " << errStr.str() << std::endl;
    }

    return valid;
}

} /* namespace meshapi */
} /* namespace asctoolkit */

