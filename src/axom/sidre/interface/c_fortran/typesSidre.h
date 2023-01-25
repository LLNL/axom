// typesSidre.h
// This file is generated by Shroud 0.12.2. Do not edit.
//
// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
// For C users and C++ implementation

#ifndef TYPESSIDRE_H
#define TYPESSIDRE_H

#ifdef __cplusplus
extern "C" {
#endif

// helper capsule_SIDRE_Buffer
struct s_SIDRE_Buffer
{
  void *addr; /* address of C++ memory */
  int idtor;  /* index of destructor */
};
typedef struct s_SIDRE_Buffer SIDRE_Buffer;

// helper capsule_SIDRE_DataStore
struct s_SIDRE_DataStore
{
  void *addr; /* address of C++ memory */
  int idtor;  /* index of destructor */
};
typedef struct s_SIDRE_DataStore SIDRE_DataStore;

// helper capsule_SIDRE_Group
struct s_SIDRE_Group
{
  void *addr; /* address of C++ memory */
  int idtor;  /* index of destructor */
};
typedef struct s_SIDRE_Group SIDRE_Group;

// helper capsule_SIDRE_View
struct s_SIDRE_View
{
  void *addr; /* address of C++ memory */
  int idtor;  /* index of destructor */
};
typedef struct s_SIDRE_View SIDRE_View;

// helper capsule_data_helper
struct s_SIDRE_SHROUD_capsule_data
{
  void *addr; /* address of C++ memory */
  int idtor;  /* index of destructor */
};
typedef struct s_SIDRE_SHROUD_capsule_data SIDRE_SHROUD_capsule_data;

void SIDRE_SHROUD_memory_destructor(SIDRE_SHROUD_capsule_data *cap);

#ifdef __cplusplus
}
#endif

#endif  // TYPESSIDRE_H
