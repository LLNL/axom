// wrapSidre.cpp
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
// wrapSidre.cpp
#define EXAMPLE_WRAPPER_IMPL
#include "wrapSidre.h"
#include "sidre/SidreTypes.hpp"

extern "C" {
namespace asctoolkit {
namespace sidre {

bool ATK_is_name_valid(const char * name)
{
// splicer begin function.is_name_valid
return name != NULL;
// splicer end function.is_name_valid
}

// splicer begin additional_functions
// splicer end additional_functions

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"
