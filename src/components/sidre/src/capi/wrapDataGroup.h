//
// Copyright (c) 2015, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
//
// All rights reserved.
//
// This source code cannot be distributed without permission and
// further review from Lawrence Livermore National Laboratory.
//
// wrapDataGroup.h
// For C users and C++ implementation

#ifndef WRAPDATAGROUP_H
#define WRAPDATAGROUP_H

#include "wrapDataView.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef EXAMPLE_WRAPPER_IMPL
typedef void ATK_datagroup;
#else
struct s_ATK_datagroup;
typedef struct s_ATK_datagroup ATK_datagroup;
#endif

const char * ATK_datagroup_get_name(const ATK_datagroup * self);

ATK_datagroup * ATK_datagroup_get_parent(ATK_datagroup * self);

ATK_dataview * ATK_datagroup_create_view_and_buffer(ATK_datagroup * self, const char * name);

ATK_datagroup * ATK_datagroup_create_group(ATK_datagroup * self, const char * name);

#ifdef __cplusplus
}
#endif

#endif  // WRAPDATAGROUP_H
