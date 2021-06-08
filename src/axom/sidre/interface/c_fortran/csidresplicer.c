// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file csidresplicer.c
 *
 * \brief   C code that will be inserted into code via shroud splicer blocks
 *
 ******************************************************************************
 */

// splicer begin C_definitions
// equivalent to C_LOC
// called from Fortran
// https://gcc.gnu.org/bugzilla/show_bug.cgi?id=53945
// Work around a problem with gfortran 4.7 where C_LOC does not work
// with assumed shape array.  Passing the first element of the
// array to a function without an interface will force the compiler
// to use f77 semantics and pass the address of the data, essentially
// the same as C_LOC.
// XXX Pass the first element, not the entire array, to avoid getting
// XXX a copy of the array.
//
// The result must be an argument because some compilers (Intel)
// cannot return type(C_PTR)
void sidre_c_loc(void* addr, void** out) { *out = addr; }
void sidre_c_loc_(void* addr, void** out) { *out = addr; }

// splicer end C_definitions
