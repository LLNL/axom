// typesSidre.h
// This is generated code, do not edit
//
// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
// For C users and C++ implementation

#ifndef TYPESSIDRE_H
#define TYPESSIDRE_H


#ifdef __cplusplus
extern "C" {
#endif

struct s_SIDRE_buffer
{
  void* addr;       /* address of C++ memory */
  int idtor;        /* index of destructor */
};
typedef struct s_SIDRE_buffer SIDRE_buffer;

struct s_SIDRE_datastore
{
  void* addr;       /* address of C++ memory */
  int idtor;        /* index of destructor */
};
typedef struct s_SIDRE_datastore SIDRE_datastore;

struct s_SIDRE_group
{
  void* addr;       /* address of C++ memory */
  int idtor;        /* index of destructor */
};
typedef struct s_SIDRE_group SIDRE_group;

struct s_SIDRE_view
{
  void* addr;       /* address of C++ memory */
  int idtor;        /* index of destructor */
};
typedef struct s_SIDRE_view SIDRE_view;

struct s_SID_SHROUD_capsule_data
{
  void* addr;       /* address of C++ memory */
  int idtor;        /* index of destructor */
};
typedef struct s_SID_SHROUD_capsule_data SID_SHROUD_capsule_data;

void SIDRE_SHROUD_memory_destructor(SID_SHROUD_capsule_data* cap);

#ifdef __cplusplus
}
#endif

#endif  // TYPESSIDRE_H
