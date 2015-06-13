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

AA_exclass1 * AA_exclass2_get_class1(AA_exclass2 * self, AA_exclass1 * in)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(self);
// splicer begin
ExClass1 * rv = selfobj->get_class1(static_cast<ExClass1 *>(in));
return rv;
// splicer end
}

void * AA_exclass2_declare(AA_exclass2 * self, int type, ATK_SidreLength len)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(self);
// splicer begin
void * rv = selfobj->declare(getTypeID(type), len);
return rv;
// splicer end
}

}  // namespace example
}  // namespace nested
}  // extern "C"
