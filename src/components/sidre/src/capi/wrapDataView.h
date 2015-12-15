// wrapDataView.h
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
// For C users and C++ implementation

#ifndef WRAPDATAVIEW_H
#define WRAPDATAVIEW_H

#include "sidre/SidreTypes.h"
#include "stdlib.h"

#ifdef __cplusplus
extern "C" {
#endif

// declaration of wrapped types
struct s_SIDRE_databuffer;
typedef struct s_SIDRE_databuffer SIDRE_databuffer;
struct s_SIDRE_datagroup;
typedef struct s_SIDRE_datagroup SIDRE_datagroup;
struct s_SIDRE_dataview;
typedef struct s_SIDRE_dataview SIDRE_dataview;

// splicer begin class.DataView.C_definition
// splicer end class.DataView.C_definition

void SIDRE_dataview_allocate_simple(SIDRE_dataview * self);

void SIDRE_dataview_allocate_from_type(SIDRE_dataview * self, int type, SIDRE_SidreLength num_elems);

void SIDRE_dataview_reallocate(SIDRE_dataview * self, SIDRE_SidreLength num_elems);

void SIDRE_dataview_apply_0(SIDRE_dataview * self);

SIDRE_dataview * SIDRE_dataview_attach_buffer(SIDRE_dataview * self, SIDRE_databuffer * buff);

void SIDRE_dataview_apply_nelems(SIDRE_dataview * self, SIDRE_SidreLength num_elems);

void SIDRE_dataview_apply_nelems_offset(SIDRE_dataview * self, SIDRE_SidreLength num_elems, SIDRE_SidreLength offset);

void SIDRE_dataview_apply_nelems_offset_stride(SIDRE_dataview * self, SIDRE_SidreLength num_elems, SIDRE_SidreLength offset, SIDRE_SidreLength stride);

void SIDRE_dataview_apply_type_nelems(SIDRE_dataview * self, int type, SIDRE_SidreLength num_elems);

void SIDRE_dataview_apply_type_nelems_offset(SIDRE_dataview * self, int type, SIDRE_SidreLength num_elems, SIDRE_SidreLength offset);

void SIDRE_dataview_apply_type_nelems_offset_stride(SIDRE_dataview * self, int type, SIDRE_SidreLength num_elems, SIDRE_SidreLength offset, SIDRE_SidreLength stride);

bool SIDRE_dataview_has_buffer(SIDRE_dataview * self);

bool SIDRE_dataview_is_opaque(SIDRE_dataview * self);

const char * SIDRE_dataview_get_name(const SIDRE_dataview * self);

SIDRE_databuffer * SIDRE_dataview_get_buffer(SIDRE_dataview * self);

void * SIDRE_dataview_get_void_ptr(SIDRE_dataview * self);

void SIDRE_dataview_set_scalar_int(SIDRE_dataview * self, int value);

void SIDRE_dataview_set_scalar_long(SIDRE_dataview * self, long value);

void SIDRE_dataview_set_scalar_float(SIDRE_dataview * self, float value);

void SIDRE_dataview_set_scalar_double(SIDRE_dataview * self, double value);

int SIDRE_dataview_get_data_int(SIDRE_dataview * self);

long SIDRE_dataview_get_data_long(SIDRE_dataview * self);

float SIDRE_dataview_get_data_float(SIDRE_dataview * self);

double SIDRE_dataview_get_data_double(SIDRE_dataview * self);

SIDRE_datagroup * SIDRE_dataview_get_owning_group(SIDRE_dataview * self);

int SIDRE_dataview_get_type_id(SIDRE_dataview * self);

size_t SIDRE_dataview_get_total_bytes(SIDRE_dataview * self);

size_t SIDRE_dataview_get_num_elements(SIDRE_dataview * self);

void SIDRE_dataview_print(SIDRE_dataview * self);

#ifdef __cplusplus
}
#endif

#endif  // WRAPDATAVIEW_H
