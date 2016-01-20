// wrapExClass1.cpp
// This is generated code, do not edit
// blah blah
// yada yada
//
// wrapExClass1.cpp
#include "wrapExClass1.h"
#include "ExClass1.hpp"
#include "shroud/shroudrt.hpp"

extern "C" {
namespace example {
namespace nested {

// ExClass1 * new(const string * name+intent(in))+constructor
// function_index=0
/**
 * \brief constructor
 *
 * longer description
 * usually multiple lines
 *
 * \return return new instance
 */
AA_exclass1 * AA_exclass1_new(const char * name)
{

// splicer begin class.ExClass1.method.new
ExClass1 * rv = new ExClass1(name);
return static_cast<AA_exclass1 *>(static_cast<void *>(rv));
// splicer end class.ExClass1.method.new
}

// ExClass1 * new(const string * name+intent(in)+len_trim(Lname))+constructor
// function_index=13
/**
 * \brief constructor
 *
 * longer description
 * usually multiple lines
 *
 * \return return new instance
 */
AA_exclass1 * AA_exclass1_new_bufferify(const char * name, int Lname)
{

// splicer begin class.ExClass1.method.new_bufferify
ExClass1 * rv = new ExClass1(std::string(name, Lname));
return static_cast<AA_exclass1 *>(static_cast<void *>(rv));
// splicer end class.ExClass1.method.new_bufferify
}

// void delete()+destructor
// function_index=1
/**
 * longer description joined with previous line
 */
void AA_exclass1_delete(AA_exclass1 * self)
{
ExClass1 *selfobj = static_cast<ExClass1 *>(static_cast<void *>(self));
// splicer begin class.ExClass1.method.delete
delete selfobj;
// splicer end class.ExClass1.method.delete
}

// int incrementCount(int incr+intent(in)+value)
// function_index=2
int AA_exclass1_increment_count(AA_exclass1 * self, int incr)
{
ExClass1 *selfobj = static_cast<ExClass1 *>(static_cast<void *>(self));
// splicer begin class.ExClass1.method.increment_count
int rv = selfobj->incrementCount(incr);
return rv;
// splicer end class.ExClass1.method.increment_count
}

// const string & getName() const
// function_index=3
const char * AA_exclass1_get_name(const AA_exclass1 * self)
{
const ExClass1 *selfobj = static_cast<const ExClass1 *>(static_cast<const void *>(self));
// splicer begin class.ExClass1.method.get_name
const std::string & rv = selfobj->getName();
if (! isNameValid(rv)) {
    return NULL;
}

return rv.c_str();
// splicer end class.ExClass1.method.get_name
}

// int GetNameLength() const
// function_index=4
/**
 * \brief helper function for Fortran
 *
 */
int AA_exclass1_get_name_length(AA_exclass1 * self)
{
ExClass1 *selfobj = static_cast<ExClass1 *>(static_cast<void *>(self));
// splicer begin class.ExClass1.method.get_name_length
return selfobj->getName().length();
// splicer end class.ExClass1.method.get_name_length
}

// const string & getNameErrorCheck() const
// function_index=5
const char * AA_exclass1_get_name_error_check(const AA_exclass1 * self)
{
const ExClass1 *selfobj = static_cast<const ExClass1 *>(static_cast<const void *>(self));
// splicer begin class.ExClass1.method.get_name_error_check
const std::string & rv = selfobj->getNameErrorCheck();
return rv.c_str();
// splicer end class.ExClass1.method.get_name_error_check
}

// const string & getNameArg() const
// function_index=6
const char * AA_exclass1_get_name_arg(const AA_exclass1 * self)
{
const ExClass1 *selfobj = static_cast<const ExClass1 *>(static_cast<const void *>(self));
// splicer begin class.ExClass1.method.get_name_arg
const std::string & rv = selfobj->getNameArg();
return rv.c_str();
// splicer end class.ExClass1.method.get_name_arg
}

// void getNameArg(string_result_as_arg & name+intent(out)+len(Lname)) const
// function_index=15
void AA_exclass1_get_name_arg_bufferify(AA_exclass1 * self, char * name, int Lname)
{
ExClass1 *selfobj = static_cast<ExClass1 *>(static_cast<void *>(self));
// splicer begin class.ExClass1.method.get_name_arg_bufferify
const std::string & rv = selfobj->getNameArg();
asctoolkit::shroud::FccCopy(name, Lname, rv.c_str());
return;
// splicer end class.ExClass1.method.get_name_arg_bufferify
}

// ExClass2 * getRoot()
// function_index=7
AA_exclass2 * AA_exclass1_get_root(AA_exclass1 * self)
{
ExClass1 *selfobj = static_cast<ExClass1 *>(static_cast<void *>(self));
// splicer begin class.ExClass1.method.get_root
ExClass2 * rv = selfobj->getRoot();
return static_cast<AA_exclass2 *>(static_cast<void *>(rv));
// splicer end class.ExClass1.method.get_root
}

// int getValue(int value+intent(in)+value)
// function_index=8
int AA_exclass1_get_value_from_int(AA_exclass1 * self, int value)
{
ExClass1 *selfobj = static_cast<ExClass1 *>(static_cast<void *>(self));
// splicer begin class.ExClass1.method.get_value_from_int
int rv = selfobj->getValue(value);
return rv;
// splicer end class.ExClass1.method.get_value_from_int
}

// long getValue(long value+intent(in)+value)
// function_index=9
long AA_exclass1_get_value_1(AA_exclass1 * self, long value)
{
ExClass1 *selfobj = static_cast<ExClass1 *>(static_cast<void *>(self));
// splicer begin class.ExClass1.method.get_value_1
long rv = selfobj->getValue(value);
return rv;
// splicer end class.ExClass1.method.get_value_1
}

// void * getAddr()
// function_index=10
void * AA_exclass1_get_addr(AA_exclass1 * self)
{
ExClass1 *selfobj = static_cast<ExClass1 *>(static_cast<void *>(self));
// splicer begin class.ExClass1.method.get_addr
void * rv = selfobj->getAddr();
return rv;
// splicer end class.ExClass1.method.get_addr
}

// bool hasAddr(bool in+intent(in)+value)
// function_index=11
bool AA_exclass1_has_addr(AA_exclass1 * self, bool in)
{
ExClass1 *selfobj = static_cast<ExClass1 *>(static_cast<void *>(self));
// splicer begin class.ExClass1.method.has_addr
bool rv = selfobj->hasAddr(in);
return rv;
// splicer end class.ExClass1.method.has_addr
}

// void SplicerSpecial()
// function_index=12
void AA_exclass1_splicer_special(AA_exclass1 * self)
{
ExClass1 *selfobj = static_cast<ExClass1 *>(static_cast<void *>(self));
// splicer begin class.ExClass1.method.splicer_special
//   splicer for SplicerSpecial
// splicer end class.ExClass1.method.splicer_special
}

// splicer begin class.ExClass1.additional_functions
// splicer end class.ExClass1.additional_functions

}  // namespace example
}  // namespace nested
}  // extern "C"
