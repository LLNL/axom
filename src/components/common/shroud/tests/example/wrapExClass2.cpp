// wrapExClass2.cpp
// This is generated code, do not edit
// blah blah
// yada yada
//
// wrapExClass2.cpp
#include "wrapExClass2.h"
#include "ExClass2.hpp"
#include "sidre/SidreWrapperHelpers.hpp"

extern "C" {
namespace example {
namespace nested {

AA_exclass2 * AA_exclass2_ex_class2(const char * name)
{
ExClass2 *selfobj = new ExClass2(name);
// splicer begin class.ExClass2.method.ex_class2
return static_cast<AA_exclass2 *>(static_cast<void *>(selfobj));
// splicer end class.ExClass2.method.ex_class2
}

AA_exclass2 * AA_exclass2_ex_class2_bufferify(const char * name, int Lname)
{
ExClass2 *selfobj = new ExClass2(std::string(name, Lname));
// splicer begin class.ExClass2.method.ex_class2_bufferify
return static_cast<AA_exclass2 *>(static_cast<void *>(selfobj));
// splicer end class.ExClass2.method.ex_class2_bufferify
}

void AA_exclass2_delete(AA_exclass2 * self)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.delete
delete selfobj;
// splicer end class.ExClass2.method.delete
}

const char * AA_exclass2_get_name(const AA_exclass2 * self)
{
const ExClass2 *selfobj = static_cast<const ExClass2 *>(static_cast<const void *>(self));
// splicer begin class.ExClass2.method.get_name
const std::string & rv = selfobj->getName();
return rv.c_str();
// splicer end class.ExClass2.method.get_name
}

int AA_exclass2_get_name_length(AA_exclass2 * self)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.get_name_length
return selfobj->getName().length();
// splicer end class.ExClass2.method.get_name_length
}

AA_exclass1 * AA_exclass2_get_class1(AA_exclass2 * self, AA_exclass1 * in)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.get_class1
ExClass1 * rv = selfobj->get_class1(static_cast<ExClass1 *>(static_cast<void *>(in)));
return static_cast<AA_exclass1 *>(static_cast<void *>(rv));
// splicer end class.ExClass2.method.get_class1
}

void AA_exclass2_declare(AA_exclass2 * self, int type, ATK_SidreLength len)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.declare
selfobj->declare(getTypeID(type), len);
return;
// splicer end class.ExClass2.method.declare
}

void AA_exclass2_destroyall(AA_exclass2 * self)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.destroyall
selfobj->destroyall();
return;
// splicer end class.ExClass2.method.destroyall
}

int AA_exclass2_get_type_id(AA_exclass2 * self)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.get_type_id
TypeID rv = selfobj->getTypeID();
return rv;
// splicer end class.ExClass2.method.get_type_id
}

void AA_exclass2_set_value_int(AA_exclass2 * self, int value)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.set_value_int
selfobj->setValue<int>(value);
return;
// splicer end class.ExClass2.method.set_value_int
}

void AA_exclass2_set_value_long(AA_exclass2 * self, long value)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.set_value_long
selfobj->setValue<long>(value);
return;
// splicer end class.ExClass2.method.set_value_long
}

void AA_exclass2_set_value_float(AA_exclass2 * self, float value)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.set_value_float
selfobj->setValue<float>(value);
return;
// splicer end class.ExClass2.method.set_value_float
}

void AA_exclass2_set_value_double(AA_exclass2 * self, double value)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.set_value_double
selfobj->setValue<double>(value);
return;
// splicer end class.ExClass2.method.set_value_double
}

int AA_exclass2_get_value_int(AA_exclass2 * self)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.get_value_int
int rv = selfobj->getValue<int>();
return rv;
// splicer end class.ExClass2.method.get_value_int
}

double AA_exclass2_get_value_double(AA_exclass2 * self)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.get_value_double
double rv = selfobj->getValue<double>();
return rv;
// splicer end class.ExClass2.method.get_value_double
}

void AA_exclass2_testoptional(AA_exclass2 * self, int i, long j)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.testoptional
selfobj->testoptional(i, j);
return;
// splicer end class.ExClass2.method.testoptional
}

// splicer begin class.ExClass2.additional_functions
// splicer end class.ExClass2.additional_functions

}  // namespace example
}  // namespace nested
}  // extern "C"
