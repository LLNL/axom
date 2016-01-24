//
// strings.hpp - wrapped routines
//

#ifndef STRINGS_HPP
#define STRINGS_HPP

#include <string>

const std::string& getName();
const std::string& getName1();
const std::string& getName2();

void acceptName_instance(std::string arg1);

void acceptStringConstReference(const std::string & arg1);

void acceptStringReference(std::string & arg1);

void acceptStringPointer(std::string * arg1);



#endif // STRINGS_HPP
