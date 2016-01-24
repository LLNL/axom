//
// strings.hpp - wrapped routines
//

#ifndef STRINGS_HPP
#define STRINGS_HPP

#include <string>

const char * getChar1();
//const char * getChar2();
//const char * getChar3();

const std::string& getString1();
const std::string& getString2();
const std::string& getString3();

void acceptName_instance(std::string arg1);

void acceptStringConstReference(const std::string & arg1);

void acceptStringReference(std::string & arg1);

void acceptStringPointer(std::string * arg1);



#endif // STRINGS_HPP
