// wrapSidre.cpp
// This is generated code, do not edit
//
// Copyright (c) 2017, Lawrence Livermore National Security, LLC.
//
// Produced at the Lawrence Livermore National Laboratory
//
// LLNL-CODE-741217
//
// All rights reserved.
//
// This file is part of Axom.
//
// For details about use and distribution, please read axom/LICENSE.
//
// wrapSidre.cpp
#include "wrapSidre.h"
#include <string>

namespace axom
{
namespace sidre
{

// splicer begin CXX_definitions
// splicer end CXX_definitions

extern "C" {

// splicer begin C_definitions
// splicer end C_definitions

bool SIDRE_name_is_valid(const char * name)
{
// splicer begin function.name_is_valid
  return name != NULL;
// splicer end function.name_is_valid
}

}  // extern "C"

}  // namespace sidre
}  // namespace axom
