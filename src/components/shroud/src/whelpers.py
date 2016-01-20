"""
Helper functions for C and Fortran wrappers.
"""

FccHeaders = """
#ifndef SHROUDRT_HPP_
#define SHROUDRT_HPP_

// Standard C++ headers
#include <cstring>

namespace asctoolkit
{
namespace shroud
{

void FccCopy(char *a, int la, const char *s)
{
   int ls,nm;
   ls = std::strlen(s);
   nm = ls < la ? ls : la;
   memcpy(a,s,nm);
   if(la > nm) { memset(a+nm,' ',la-nm);}
}

} /* end namespace shroud */
} /* end namespace asctoolkit */

#endif /* SHROUDRT_HPP_ */
"""
