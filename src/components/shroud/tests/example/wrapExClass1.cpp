// wrapExClass1.cpp
// This is generated code, do not edit
// blah blah
// yada yada
//
// wrapExClass1.cpp
#include "wrapExClass1.h"
#include <string>
#include "ExClass1.hpp"
#include "shroudrt.hpp"

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
    const std::string SH_name(name);
    ExClass1 * SH_rv = new ExClass1(SH_name);
    return static_cast<AA_exclass1 *>(static_cast<void *>(SH_rv));
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
    const std::string SH_name(name, Lname);
    ExClass1 * SH_rv = new ExClass1(SH_name);
    return static_cast<AA_exclass1 *>(static_cast<void *>(SH_rv));
// splicer end class.ExClass1.method.new_bufferify
}

// void delete()+destructor
// function_index=1
/**
 * longer description joined with previous line
 */
void AA_exclass1_delete(AA_exclass1 * self)
{
// splicer begin class.ExClass1.method.delete
    ExClass1 *SH_this = static_cast<ExClass1 *>(static_cast<void *>(self));
    delete SH_this;
    return;
// splicer end class.ExClass1.method.delete
}

// int incrementCount(int incr+intent(in)+value)
// function_index=2
int AA_exclass1_increment_count(AA_exclass1 * self, int incr)
{
// splicer begin class.ExClass1.method.increment_count
    ExClass1 *SH_this = static_cast<ExClass1 *>(static_cast<void *>(self));
    int SH_rv = SH_this->incrementCount(incr);
    return SH_rv;
// splicer end class.ExClass1.method.increment_count
}

// const string & getName() const
// function_index=3
const char * AA_exclass1_get_name(const AA_exclass1 * self)
{
// splicer begin class.ExClass1.method.get_name
    const ExClass1 *SH_this = static_cast<const ExClass1 *>(static_cast<const void *>(self));
    const std::string & SH_rv = SH_this->getName();
    // check for error
    if (! isNameValid(SH_rv)) {
        return NULL;
    }

    return SH_rv.c_str();
// splicer end class.ExClass1.method.get_name
}

// void getName(string & SH_F_rv+intent(out)+len(LSH_F_rv)) const
// function_index=14
void AA_exclass1_get_name_bufferify(const AA_exclass1 * self, char * SH_F_rv, int LSH_F_rv)
{
// splicer begin class.ExClass1.method.get_name_bufferify
    const ExClass1 *SH_this = static_cast<const ExClass1 *>(static_cast<const void *>(self));
    const std::string & SH_rv = SH_this->getName();
    shroud_FccCopy(SH_F_rv, LSH_F_rv, SH_rv.c_str());
    return;
// splicer end class.ExClass1.method.get_name_bufferify
}

// int GetNameLength() const
// function_index=4
/**
 * \brief helper function for Fortran to get length of name.
 *
 */
int AA_exclass1_get_name_length(const AA_exclass1 * self)
{
// splicer begin class.ExClass1.method.get_name_length
    const ExClass1 *SH_this = static_cast<const ExClass1 *>(static_cast<const void *>(self));
    return SH_this->getName().length();

// splicer end class.ExClass1.method.get_name_length
}

// const string & getNameErrorCheck() const
// function_index=5
const char * AA_exclass1_get_name_error_check(const AA_exclass1 * self)
{
// splicer begin class.ExClass1.method.get_name_error_check
    const ExClass1 *SH_this = static_cast<const ExClass1 *>(static_cast<const void *>(self));
    const std::string & SH_rv = SH_this->getNameErrorCheck();
    return SH_rv.c_str();
// splicer end class.ExClass1.method.get_name_error_check
}

// void getNameErrorCheck(string & SH_F_rv+intent(out)+len(LSH_F_rv)) const
// function_index=15
void AA_exclass1_get_name_error_check_bufferify(const AA_exclass1 * self, char * SH_F_rv, int LSH_F_rv)
{
// splicer begin class.ExClass1.method.get_name_error_check_bufferify
    const ExClass1 *SH_this = static_cast<const ExClass1 *>(static_cast<const void *>(self));
    const std::string & SH_rv = SH_this->getNameErrorCheck();
    shroud_FccCopy(SH_F_rv, LSH_F_rv, SH_rv.c_str());
    return;
// splicer end class.ExClass1.method.get_name_error_check_bufferify
}

// const string & getNameArg() const
// function_index=6
const char * AA_exclass1_get_name_arg(const AA_exclass1 * self)
{
// splicer begin class.ExClass1.method.get_name_arg
    const ExClass1 *SH_this = static_cast<const ExClass1 *>(static_cast<const void *>(self));
    const std::string & SH_rv = SH_this->getNameArg();
    return SH_rv.c_str();
// splicer end class.ExClass1.method.get_name_arg
}

// void getNameArg(string & name+intent(out)+len(Lname)) const
// function_index=16
void AA_exclass1_get_name_arg_bufferify(const AA_exclass1 * self, char * name, int Lname)
{
// splicer begin class.ExClass1.method.get_name_arg_bufferify
    const ExClass1 *SH_this = static_cast<const ExClass1 *>(static_cast<const void *>(self));
    const std::string & SH_rv = SH_this->getNameArg();
    shroud_FccCopy(name, Lname, SH_rv.c_str());
    return;
// splicer end class.ExClass1.method.get_name_arg_bufferify
}

// ExClass2 * getRoot()
// function_index=7
AA_exclass2 * AA_exclass1_get_root(AA_exclass1 * self)
{
// splicer begin class.ExClass1.method.get_root
    ExClass1 *SH_this = static_cast<ExClass1 *>(static_cast<void *>(self));
    ExClass2 * SH_rv = SH_this->getRoot();
    return static_cast<AA_exclass2 *>(static_cast<void *>(SH_rv));
// splicer end class.ExClass1.method.get_root
}

// int getValue(int value+intent(in)+value)
// function_index=8
int AA_exclass1_get_value_from_int(AA_exclass1 * self, int value)
{
// splicer begin class.ExClass1.method.get_value_from_int
    ExClass1 *SH_this = static_cast<ExClass1 *>(static_cast<void *>(self));
    int SH_rv = SH_this->getValue(value);
    return SH_rv;
// splicer end class.ExClass1.method.get_value_from_int
}

// long getValue(long value+intent(in)+value)
// function_index=9
long AA_exclass1_get_value_1(AA_exclass1 * self, long value)
{
// splicer begin class.ExClass1.method.get_value_1
    ExClass1 *SH_this = static_cast<ExClass1 *>(static_cast<void *>(self));
    long SH_rv = SH_this->getValue(value);
    return SH_rv;
// splicer end class.ExClass1.method.get_value_1
}

// void * getAddr()
// function_index=10
void * AA_exclass1_get_addr(AA_exclass1 * self)
{
// splicer begin class.ExClass1.method.get_addr
    ExClass1 *SH_this = static_cast<ExClass1 *>(static_cast<void *>(self));
    void * SH_rv = SH_this->getAddr();
    return SH_rv;
// splicer end class.ExClass1.method.get_addr
}

// bool hasAddr(bool in+intent(in)+value)
// function_index=11
bool AA_exclass1_has_addr(AA_exclass1 * self, bool in)
{
// splicer begin class.ExClass1.method.has_addr
    ExClass1 *SH_this = static_cast<ExClass1 *>(static_cast<void *>(self));
    bool SH_rv = SH_this->hasAddr(in);
    return SH_rv;
// splicer end class.ExClass1.method.has_addr
}

// void SplicerSpecial()
// function_index=12
void AA_exclass1_splicer_special(AA_exclass1 * self)
{
// splicer begin class.ExClass1.method.splicer_special
//   splicer for SplicerSpecial
// splicer end class.ExClass1.method.splicer_special
}

// splicer begin class.ExClass1.additional_functions
// splicer end class.ExClass1.additional_functions

}  // namespace example
}  // namespace nested
}  // extern "C"
