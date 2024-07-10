// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


#include "axom/sina/tests/include/ConduitTestUtils.hpp"

namespace axom
{
namespace sina
{
namespace testing
{

conduit::Node parseJsonValue(std::string const &valueAsString) {
    // If we just try to do node.parse(valueAsString, "json"), then passing
    // in strings does not work. We need to create a document with a key
    // so that valueAsString can be parsed as a value.
    conduit::Node node;
    std::string fullContents = "{\"TEST_KEY\": ";
    fullContents += valueAsString;
    fullContents += "}";
    node.parse(fullContents, "json");
    return node.child("TEST_KEY");
}

}  // end testing namespace
}  // end sina namespace
}  // end axom namespace
