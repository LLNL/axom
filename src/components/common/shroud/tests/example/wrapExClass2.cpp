// blah blah
// yada yada
//
// wrapExClass2.cpp
#define EXAMPLE_WRAPPER_IMPL
#include "wrapExClass2.h"
#include "ExClass2.hpp"
#include "sidre/SidreWrapperHelpers.hpp"

extern "C" {
namespace example {
namespace nested {
// splicer push class.ExClass2.method

AA_exclass2 * AA_exclass2_ex_class2(const char * name)
{
ExClass2 *selfobj = new ExClass2(name);
// splicer begin AA_exclass2_ex_class2
return (AA_exclass2 *) selfobj;
// splicer end AA_exclass2_ex_class2
}

void AA_exclass2_ex_class1(AA_exclass2 * self)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(self);
// splicer begin AA_exclass2_ex_class1
delete selfobj;
// splicer end AA_exclass2_ex_class1
}

const char * AA_exclass2_get_name(const AA_exclass2 * self)
{
const ExClass2 *selfobj = static_cast<const ExClass2 *>(self);
// splicer begin AA_exclass2_get_name
const std::string & rv = selfobj->getName();
return isNameValid(rv) ? rv.c_str() : NULL;
// splicer end AA_exclass2_get_name
}

const int AA_exclass2_get_name_length(const AA_exclass2 * self)
{
const ExClass2 *selfobj = static_cast<const ExClass2 *>(self);
// splicer begin AA_exclass2_get_name_length
return selfobj->getName().length();
// splicer end AA_exclass2_get_name_length
}

AA_exclass1 * AA_exclass2_get_class1(AA_exclass2 * self, AA_exclass1 * in)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(self);
// splicer begin AA_exclass2_get_class1
ExClass1 * rv = selfobj->get_class1(static_cast<ExClass1 *>(in));
return rv;
// splicer end AA_exclass2_get_class1
}

void AA_exclass2_declare(AA_exclass2 * self, int type, ATK_SidreLength len)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(self);
// splicer begin AA_exclass2_declare
selfobj->declare(getTypeID(type), len);
return;
// splicer end AA_exclass2_declare
}

void AA_exclass2_destroyall(AA_exclass2 * self)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(self);
// splicer begin AA_exclass2_destroyall
selfobj->destroyall();
return;
// splicer end AA_exclass2_destroyall
}

int AA_exclass2_get_type_id(AA_exclass2 * self)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(self);
// splicer begin AA_exclass2_get_type_id
TypeID rv = selfobj->getTypeID();
return rv;
// splicer end AA_exclass2_get_type_id
}

// splicer pop.class.ExClass2 method

}  // namespace example
}  // namespace nested
}  // extern "C"
