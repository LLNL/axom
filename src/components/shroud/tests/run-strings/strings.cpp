//
// tutorial.hpp - wrapped routines
//

#include "strings.hpp"

static std::string last_function_called;

// These variables exist to avoid warning errors
static std::string static_str = std::string("dog");
static std::string global_str;


const std::string& getName()
{
    return static_str;
}
