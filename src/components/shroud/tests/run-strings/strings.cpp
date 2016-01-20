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

void acceptStringConstReference(const std::string & arg1)
{
    global_str = arg1;
}

void acceptStringReference(std::string & arg1)
{
    global_str = arg1;
}

void acceptStringPointer(std::string * arg1)
{
    global_str = *arg1;
}

