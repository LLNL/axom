//
// tutorial.hpp - wrapped routines
//

#include "strings.hpp"

static std::string last_function_called;

// These variables exist to avoid warning errors
static const char * static_char = "bird";
static std::string static_str = std::string("dog");
static std::string global_str;


const char * getChar1()
{
    return static_char;
}

const char * getChar2()
{
    return static_char;
}

const char * getChar3()
{
    return static_char;
}

//----------------------------------------

const std::string& getString1()
{
    return static_str;
}

const std::string& getString2()
{
    return static_str;
}

const std::string& getString3()
{
    return static_str;
}

//----------------------------------------

void acceptStringConstReference(const std::string & arg1)
{
    global_str = arg1;
}

void acceptStringReference(std::string & arg1)
{
    arg1.append("dog");
}

void acceptStringPointer(std::string * arg1)
{
    global_str = *arg1;
}

