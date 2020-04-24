// typesSLIC.h
// This is generated code, do not edit
//
// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
// For C users and C++ implementation

#ifndef TYPESSLIC_H
#define TYPESSLIC_H


#ifdef __cplusplus
extern "C" {
#endif

struct s_SLIC_SHROUD_capsule_data
{
  void* addr;       /* address of C++ memory */
  int idtor;        /* index of destructor */
};
typedef struct s_SLIC_SHROUD_capsule_data SLIC_SHROUD_capsule_data;

void SLIC_SHROUD_memory_destructor(SLIC_SHROUD_capsule_data* cap);

#ifdef __cplusplus
}
#endif

#endif  // TYPESSLIC_H
