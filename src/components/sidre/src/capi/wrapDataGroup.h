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

#include "sidre/SidreTypes.h"
#include "stdlib.h"

#ifdef __cplusplus
extern "C" {
#endif

// declaration of wrapped types
#ifdef EXAMPLE_WRAPPER_IMPL
typedef void ATK_databuffer;
typedef void ATK_datagroup;
typedef void ATK_datastore;
typedef void ATK_dataview;
#else
struct s_ATK_databuffer;
typedef struct s_ATK_databuffer ATK_databuffer;
struct s_ATK_datagroup;
typedef struct s_ATK_datagroup ATK_datagroup;
struct s_ATK_datastore;
typedef struct s_ATK_datastore ATK_datastore;
struct s_ATK_dataview;
typedef struct s_ATK_dataview ATK_dataview;
#endif

const char * ATK_datagroup_get_name(const ATK_datagroup * self);

const ATK_datagroup * ATK_datagroup_get_parent(const ATK_datagroup * self);

const ATK_datastore * ATK_datagroup_get_data_store(const ATK_datagroup * self);

size_t ATK_datagroup_get_num_views(ATK_datagroup * self);

size_t ATK_datagroup_get_num_groups(ATK_datagroup * self);

bool ATK_datagroup_has_view(ATK_datagroup * self, const char * name);

ATK_dataview * ATK_datagroup_create_view_and_buffer(ATK_datagroup * self, const char * name);

ATK_dataview * ATK_datagroup_create_view_and_buffer_from_type(ATK_datagroup * self, const char * name, int type, ATK_SidreLength len);

ATK_dataview * ATK_datagroup_create_opaque_view(ATK_datagroup * self, const char * name, void * opaque_ptr);

ATK_dataview * ATK_datagroup_create_view(ATK_datagroup * self, const char * name, ATK_databuffer * buff);

ATK_dataview * ATK_datagroup_move_view(ATK_datagroup * self, ATK_dataview * view);

ATK_dataview * ATK_datagroup_copy_view(ATK_datagroup * self, ATK_dataview * view);

void ATK_datagroup_destroy_view_and_buffer(ATK_datagroup * self, const char * name);

ATK_dataview * ATK_datagroup_get_view(ATK_datagroup * self, const char * name);

ATK_IndexType ATK_datagroup_get_view_index(ATK_datagroup * self, const char * name);

const char * ATK_datagroup_get_view_name(const ATK_datagroup * self, ATK_IndexType idx);

bool ATK_datagroup_has_group(ATK_datagroup * self, const char * name);

ATK_datagroup * ATK_datagroup_create_group(ATK_datagroup * self, const char * name);

ATK_datagroup * ATK_datagroup_move_group(ATK_datagroup * self, ATK_datagroup * grp);

void ATK_datagroup_destroy_group(ATK_datagroup * self, const char * name);

ATK_datagroup * ATK_datagroup_get_group(ATK_datagroup * self, const char * name);

ATK_IndexType ATK_datagroup_get_group_index(ATK_datagroup * self, const char * name);

const char * ATK_datagroup_get_group_name(const ATK_datagroup * self, ATK_IndexType idx);

void ATK_datagroup_print(ATK_datagroup * self);

void ATK_datagroup_save(ATK_datagroup * self, const char * obase, const char * protocol);

void ATK_datagroup_load(ATK_datagroup * self, const char * obase, const char * protocol);

#ifdef __cplusplus
}
#endif

#endif  // WRAPDATAGROUP_H
