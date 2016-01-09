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

void SLIC_set_abort_on_assert(bool willAbort)
{
// splicer begin function.set_abort_on_assert
setAbortOnAssert(willAbort);
return;
// splicer end function.set_abort_on_assert
}

bool SLIC_get_abort_on_assert()
{
// splicer begin function.get_abort_on_assert
bool rv = getAbortOnAssert();
return rv;
// splicer end function.get_abort_on_assert
}

void SLIC_set_abort_on_error(bool willAbort)
{
// splicer begin function.set_abort_on_error
setAbortOnError(willAbort);
return;
// splicer end function.set_abort_on_error
}

bool SLIC_get_abort_on_error()
{
// splicer begin function.get_abort_on_error
bool rv = getAbortOnError();
return rv;
// splicer end function.get_abort_on_error
}

// splicer begin additional_functions
// splicer end additional_functions

}  // namespace asctoolkit
}  // namespace slic
}  // extern "C"
