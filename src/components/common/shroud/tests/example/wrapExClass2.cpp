// wrapExClass2.cpp
// This is generated code, do not edit
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

AA_exclass2 * AA_exclass2_ex_class2(const char * name)
{
ExClass2 *selfobj = new ExClass2(name);
// splicer begin class.ExClass2.method.ExClass2
return (AA_exclass2 *) selfobj;
// splicer end class.ExClass2.method.ExClass2
}

void AA_exclass2_ex_class1(AA_exclass2 * self)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(self);
// splicer begin class.ExClass2.method.ExClass1
delete selfobj;
// splicer end class.ExClass2.method.ExClass1
}

const char * AA_exclass2_get_name(const AA_exclass2 * self)
{
const ExClass2 *selfobj = static_cast<const ExClass2 *>(self);
// splicer begin class.ExClass2.method.getName
const std::string & rv = selfobj->getName();
return isNameValid(rv) ? rv.c_str() : NULL;
// splicer end class.ExClass2.method.getName
}

const int AA_exclass2_get_name_length(const AA_exclass2 * self)
{
const ExClass2 *selfobj = static_cast<const ExClass2 *>(self);
// splicer begin class.ExClass2.method.get_name_length
return selfobj->getName().length();
// splicer end class.ExClass2.method.get_name_length
}

AA_exclass1 * AA_exclass2_get_class1(AA_exclass2 * self, AA_exclass1 * in)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(self);
// splicer begin class.ExClass2.method.get_class1
ExClass1 * rv = selfobj->get_class1(static_cast<ExClass1 *>(in));
return rv;
// splicer end class.ExClass2.method.get_class1
}

void AA_exclass2_declare(AA_exclass2 * self, int type, ATK_SidreLength len)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(self);
// splicer begin class.ExClass2.method.declare
selfobj->declare(getTypeID(type), len);
return;
// splicer end class.ExClass2.method.declare
}

void AA_exclass2_destroyall(AA_exclass2 * self)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(self);
// splicer begin class.ExClass2.method.destroyall
selfobj->destroyall();
return;
// splicer end class.ExClass2.method.destroyall
}

int AA_exclass2_get_type_id(AA_exclass2 * self)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(self);
// splicer begin class.ExClass2.method.getTypeID
TypeID rv = selfobj->getTypeID();
return rv;
// splicer end class.ExClass2.method.getTypeID
}

// splicer begin class.ExClass2.additional_functions
// splicer end class.ExClass2.additional_functions

}  // namespace example
}  // namespace nested
}  // extern "C"
