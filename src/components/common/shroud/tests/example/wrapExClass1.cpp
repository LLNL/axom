// blah blah
// yada yada
//
// wrapExClass1.cpp
#define EXAMPLE_WRAPPER_IMPL
#include "wrapExClass1.h"
#include "ExClass1.hpp"

extern "C" {
namespace example {
namespace nested {
// splicer push class.ExClass1.method

AA_exclass1 * AA_exclass1_new(const char * name)
{
ExClass1 *selfobj = new ExClass1(name);
// splicer begin AA_exclass1_new
return (AA_exclass1 *) selfobj;
// splicer end AA_exclass1_new
}

void AA_exclass1_delete(AA_exclass1 * self)
{
ExClass1 *selfobj = static_cast<ExClass1 *>(self);
// splicer begin AA_exclass1_delete
delete selfobj;
// splicer end AA_exclass1_delete
}

int AA_exclass1_increment_count(AA_exclass1 * self, int incr)
{
ExClass1 *selfobj = static_cast<ExClass1 *>(self);
// splicer begin AA_exclass1_increment_count
int rv = selfobj->incrementCount(incr);
return rv;
// splicer end AA_exclass1_increment_count
}

const char * AA_exclass1_get_name(const AA_exclass1 * self)
{
const ExClass1 *selfobj = static_cast<const ExClass1 *>(self);
// splicer begin AA_exclass1_get_name
const std::string & rv = selfobj->getName();
return rv.c_str();
// splicer end AA_exclass1_get_name
}

const int AA_exclass1_get_name_length(const AA_exclass1 * self)
{
const ExClass1 *selfobj = static_cast<const ExClass1 *>(self);
// splicer begin AA_exclass1_get_name_length
return selfobj->getName().length();
// splicer end AA_exclass1_get_name_length
}

AA_exclass2 * AA_exclass1_get_root(AA_exclass1 * self)
{
ExClass1 *selfobj = static_cast<ExClass1 *>(self);
// splicer begin AA_exclass1_get_root
ExClass2 * rv = selfobj->getRoot();
return rv;
// splicer end AA_exclass1_get_root
}

int AA_exclass1_get_value_from_int(AA_exclass1 * self, int value)
{
ExClass1 *selfobj = static_cast<ExClass1 *>(self);
// splicer begin AA_exclass1_get_value_from_int
int rv = selfobj->getValue(value);
return rv;
// splicer end AA_exclass1_get_value_from_int
}

long AA_exclass1_get_value_1(AA_exclass1 * self, long value)
{
ExClass1 *selfobj = static_cast<ExClass1 *>(self);
// splicer begin AA_exclass1_get_value_1
long rv = selfobj->getValue(value);
return rv;
// splicer end AA_exclass1_get_value_1
}

void * AA_exclass1_get_addr(AA_exclass1 * self)
{
ExClass1 *selfobj = static_cast<ExClass1 *>(self);
// splicer begin AA_exclass1_get_addr
void * rv = selfobj->getAddr();
return rv;
// splicer end AA_exclass1_get_addr
}

bool AA_exclass1_has_addr(AA_exclass1 * self, bool in)
{
ExClass1 *selfobj = static_cast<ExClass1 *>(self);
// splicer begin AA_exclass1_has_addr
bool rv = selfobj->hasAddr(in);
return rv;
// splicer end AA_exclass1_has_addr
}

// splicer pop.class.ExClass1 method

}  // namespace example
}  // namespace nested
}  // extern "C"
