// blah blah
// yada yada
//
// wrapExClass2.cpp
#define EXAMPLE_WRAPPER_IMPL
#include "wrapExClass2.h"
#include "ExClass2.hpp"

extern "C" {
namespace example {
namespace nested {

AA_exclass2 * AA_exclass2_ex_class2(const char * name)
{
ExClass2 *selfobj = new ExClass2(name);
// splicer begin
return (AA_exclass2 *) selfobj;
// splicer end
}

void AA_exclass2_ex_class1(AA_exclass2 * self)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(self);
// splicer begin
delete selfobj;
// splicer end
}

const char * AA_exclass2_get_name(const AA_exclass2 * self)
{
const ExClass2 *selfobj = static_cast<const ExClass2 *>(self);
// splicer begin
const std::string & rv = selfobj->getName();
return rv.c_str();
// splicer end
}

const int AA_exclass2_get_name_length(const AA_exclass2 * self)
{
const ExClass2 *selfobj = static_cast<const ExClass2 *>(self);
// splicer begin
return selfobj->getName().length();
// splicer end
}

}  // namespace example
}  // namespace nested
}  // extern "C"
