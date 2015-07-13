#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <string>

#ifndef ENABLE_CXX11
#include <sstream>
#endif

namespace asctoolkit {
namespace lumberjack {

std::string intToString(int intValue)
{
    std::string stringValue = "";
#ifdef ENABLE_CXX11
    stringValue += std::to_string(intValue);
#else
    std::ostringstream ss;
    ss << intValue;
    stringValue += ss.str();
#endif
    return stringValue;
}

int stringToInt(const std::string& stringValue)
{
    int intValue = 0;
#ifdef ENABLE_CXX11
    intValue = stoi(stringValue);
#else
    std::istringstream(stringValue) >> intValue;
#endif
    return intValue;
}

}
}

#endif