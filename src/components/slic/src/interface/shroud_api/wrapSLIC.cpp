// wrapSLIC.cpp
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
// wrapSLIC.cpp
#include "wrapSLIC.h"
#include "slic/slic.hpp"

extern "C" {
namespace asctoolkit {
namespace slic {

void SLIC_initialize()
{
// splicer begin function.initialize
initialize();
return;
// splicer end function.initialize
}

bool SLIC_is_initialized()
{
// splicer begin function.is_initialized
bool rv = isInitialized();
return rv;
// splicer end function.is_initialized
}

void SLIC_finalize()
{
// splicer begin function.finalize
finalize();
return;
// splicer end function.finalize
}

// splicer begin additional_functions
// splicer end additional_functions

}  // namespace asctoolkit
}  // namespace slic
}  // extern "C"
