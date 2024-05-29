#ifndef SINA_CONDUITTESTUTILS_HPP
#define SINA_CONDUITTESTUTILS_HPP

#include <string>

#include "gmock/gmock.h"

#include "conduit.hpp"

namespace axom
{
namespace sina
{
namespace testing
{

/**
 * Parse a JSON value.
 *
 * @param valueAsString the value as a string
 * @return the value as a Conduit node
 */
conduit::Node parseJsonValue(std::string const &valueAsString);

// A matcher which verifies that a given Conduit node produces the expected
// JSON string
MATCHER_P(MatchesJson, expectedJsonString, "") {
    conduit::Node expected = parseJsonValue(expectedJsonString);
    *result_listener << "Given node is " << arg.to_json_default();
    conduit::Node diff;
    bool differ = expected.diff(arg, diff);
    return !differ;
}

}  // end testing namespace
}  // end sina namespace
}  // end axom namespace

#endif //SINA_CONDUITTESTUTILS_HPP
