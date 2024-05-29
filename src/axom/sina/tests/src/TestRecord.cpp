#include "axom/sina/tests/include/TestRecord.hpp"

namespace axom
{
namespace sina
{
namespace testing
{

template<>
TestRecord<std::string>::TestRecord(conduit::Node const &asNode) :
        Record{asNode},
        value{getRequiredString(TEST_RECORD_VALUE_KEY, asNode, "TestRecord")} {}

template<>
TestRecord<int>::TestRecord(conduit::Node const &asNode) :
        Record{asNode},
        value{getRequiredField(TEST_RECORD_VALUE_KEY, asNode,
              "TestRecord").as_int()} {}

}  // end testing namespace
}  // end sina namespace
}  // end axom namespace
