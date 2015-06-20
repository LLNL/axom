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

bool ATK_is_name_valid(const char * name)
{
// splicer begin function.isNameValid
return name != NULL;
// splicer end function.isNameValid
}

}  // extern "C"
