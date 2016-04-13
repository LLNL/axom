// shroudrt.hpp
// This is generated code, do not edit
//
// Copyright (c) 2015, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
//
// All rights reserved.
//
// This source code cannot be distributed without permission and
// further review from Lawrence Livermore National Laboratory.
//

#ifndef SHROUDRT_HPP_
#define SHROUDRT_HPP_

// Standard C++ headers
#include <cstring>

namespace asctoolkit
{
namespace shroud
{

static inline void FccCopy(char *a, int la, const char *s)
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

