#pragma once

#define __STRINGIZE(x) __STRINGIZE2(x)
#define __STRINGIZE2(x) #x
//!@brief String literal for code location
#define __WHERE \
  __FILE__ ":" __STRINGIZE(__LINE__) "(" + std::string(__func__) + ") "
