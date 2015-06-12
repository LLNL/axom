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

AA_exclass1 * AA_exclass1_new(const char * name)
{
ExClass1 *selfobj = new ExClass1(name);
// splicer begin
return (AA_exclass1 *) selfobj;
// splicer end
}

void AA_exclass1_delete(AA_exclass1 * self)
{
ExClass1 *selfobj = static_cast<ExClass1 *>(self);
// splicer begin
delete selfobj;
// splicer end
}

int AA_exclass1_increment_count(AA_exclass1 * self, int incr)
{
ExClass1 *selfobj = static_cast<ExClass1 *>(self);
// splicer begin
int rv = selfobj->incrementCount(incr);
return rv;
// splicer end
}

const char * AA_exclass1_get_name(const AA_exclass1 * self)
{
const ExClass1 *selfobj = static_cast<const ExClass1 *>(self);
// splicer begin
const std::string & rv = selfobj->getName();
return rv.c_str();
// splicer end
}

const int AA_exclass1_get_name_length(const AA_exclass1 * self)
{
const ExClass1 *selfobj = static_cast<const ExClass1 *>(self);
// splicer begin
return selfobj->getName().length();
// splicer end
}

AA_exclass2 * AA_exclass1_get_root(AA_exclass1 * self)
{
ExClass1 *selfobj = static_cast<ExClass1 *>(self);
// splicer begin
ExClass2 * rv = selfobj->getRoot();
return rv;
// splicer end
}

int AA_exclass1_get_value_from_int(AA_exclass1 * self, int value)
{
ExClass1 *selfobj = static_cast<ExClass1 *>(self);
// splicer begin
int rv = selfobj->getValue(value);
return rv;
// splicer end
}

long AA_exclass1_get_value_1(AA_exclass1 * self, long value)
{
ExClass1 *selfobj = static_cast<ExClass1 *>(self);
// splicer begin
long rv = selfobj->getValue(value);
return rv;
// splicer end
}

void * AA_exclass1_get_addr(AA_exclass1 * self)
{
ExClass1 *selfobj = static_cast<ExClass1 *>(self);
// splicer begin
void * rv = selfobj->getAddr();
return rv;
// splicer end
}

}  // namespace example
}  // namespace nested
}  // extern "C"
