//
// strings.hpp - wrapped routines
//

#ifndef STRINGS_HPP
#define STRINGS_HPP

#include <string>

const std::string& getName();

void acceptName_instance(std::string arg1);

void acceptStringConstReference(const std::string & arg1);

void acceptStringReference(std::string & arg1);

void acceptStringPointer(std::string * arg1);



#endif // STRINGS_HPP
