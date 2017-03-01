// wrapDataGroup.h
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
/**
 * \file wrapDataGroup.h
 * \brief Shroud generated wrapper for DataGroup class
 */
// For C users and C++ implementation

#ifndef WRAPDATAGROUP_H
#define WRAPDATAGROUP_H

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
struct s_SIDRE_datastore;
typedef struct s_SIDRE_datastore SIDRE_datastore;
struct s_SIDRE_dataview;
typedef struct s_SIDRE_dataview SIDRE_dataview;

// splicer begin class.DataGroup.C_definition
// splicer end class.DataGroup.C_definition

const char * SIDRE_datagroup_get_name(const SIDRE_datagroup * self);

void SIDRE_datagroup_get_name_bufferify(const SIDRE_datagroup * self,
                                        char * SH_F_rv, int LSH_F_rv);

const SIDRE_datagroup * SIDRE_datagroup_get_parent(const SIDRE_datagroup * self);

const SIDRE_datastore * SIDRE_datagroup_get_data_store(
  const SIDRE_datagroup * self);

size_t SIDRE_datagroup_get_num_views(const SIDRE_datagroup * self);

size_t SIDRE_datagroup_get_num_groups(const SIDRE_datagroup * self);

bool SIDRE_datagroup_has_view(const SIDRE_datagroup * self, const char * path);

bool SIDRE_datagroup_has_view_bufferify(const SIDRE_datagroup * self,
                                        const char * path, int Lpath);

bool SIDRE_datagroup_has_child_view(const SIDRE_datagroup * self,
                                    const char * name);

bool SIDRE_datagroup_has_child_view_bufferify(const SIDRE_datagroup * self,
                                              const char * name, int Lname);

SIDRE_dataview * SIDRE_datagroup_get_view_from_name(SIDRE_datagroup * self,
                                                    const char * path);

SIDRE_dataview * SIDRE_datagroup_get_view_from_name_bufferify(
  SIDRE_datagroup * self, const char * path, int Lpath);

SIDRE_dataview * SIDRE_datagroup_get_view_from_index(SIDRE_datagroup * self,
                                                     const SIDRE_IndexType idx);

SIDRE_IndexType SIDRE_datagroup_get_view_index(const SIDRE_datagroup * self,
                                               const char * name);

SIDRE_IndexType SIDRE_datagroup_get_view_index_bufferify(
  const SIDRE_datagroup * self, const char * name, int Lname);

const char * SIDRE_datagroup_get_view_name(const SIDRE_datagroup * self,
                                           SIDRE_IndexType idx);

void SIDRE_datagroup_get_view_name_bufferify(const SIDRE_datagroup * self,
                                             SIDRE_IndexType idx,
                                             char * SH_F_rv, int LSH_F_rv);

SIDRE_IndexType SIDRE_datagroup_get_first_valid_view_index(
  const SIDRE_datagroup * self);

SIDRE_IndexType SIDRE_datagroup_get_next_valid_view_index(
  const SIDRE_datagroup * self, SIDRE_IndexType idx);

SIDRE_dataview * SIDRE_datagroup_create_view_and_allocate_nelems(
  SIDRE_datagroup * self, const char * path, int type,
  SIDRE_SidreLength num_elems);

SIDRE_dataview * SIDRE_datagroup_create_view_and_allocate_nelems_bufferify(
  SIDRE_datagroup * self, const char * path, int Lpath, int type,
  SIDRE_SidreLength num_elems);

SIDRE_dataview * SIDRE_datagroup_create_view_and_allocate_shape(
  SIDRE_datagroup * self, const char * path, int type, int ndims,
  SIDRE_SidreLength * shape);

SIDRE_dataview * SIDRE_datagroup_create_view_and_allocate_shape_bufferify(
  SIDRE_datagroup * self, const char * path, int Lpath, int type, int ndims,
  SIDRE_SidreLength * shape);

SIDRE_dataview * SIDRE_datagroup_create_view_scalar_int(SIDRE_datagroup * self,
                                                        const char * path,
                                                        int value);

SIDRE_dataview * SIDRE_datagroup_create_view_scalar_int_bufferify(
  SIDRE_datagroup * self, const char * path, int Lpath, int value);

SIDRE_dataview * SIDRE_datagroup_create_view_scalar_long(SIDRE_datagroup * self,
                                                         const char * path,
                                                         long value);

SIDRE_dataview * SIDRE_datagroup_create_view_scalar_long_bufferify(
  SIDRE_datagroup * self, const char * path, int Lpath, long value);

SIDRE_dataview * SIDRE_datagroup_create_view_scalar_float(
  SIDRE_datagroup * self, const char * path, float value);

SIDRE_dataview * SIDRE_datagroup_create_view_scalar_float_bufferify(
  SIDRE_datagroup * self, const char * path, int Lpath, float value);

SIDRE_dataview * SIDRE_datagroup_create_view_scalar_double(
  SIDRE_datagroup * self, const char * path, double value);

SIDRE_dataview * SIDRE_datagroup_create_view_scalar_double_bufferify(
  SIDRE_datagroup * self, const char * path, int Lpath, double value);

SIDRE_dataview * SIDRE_datagroup_create_view_string(SIDRE_datagroup * self,
                                                    const char * path,
                                                    const char * value);

SIDRE_dataview * SIDRE_datagroup_create_view_string_bufferify(
  SIDRE_datagroup * self, const char * path, int Lpath, const char * value,
  int Lvalue);

SIDRE_dataview * SIDRE_datagroup_create_view_empty(SIDRE_datagroup * self,
                                                   const char * path);

SIDRE_dataview * SIDRE_datagroup_create_view_empty_bufferify(
  SIDRE_datagroup * self, const char * path, int Lpath);

SIDRE_dataview * SIDRE_datagroup_create_view_from_type(SIDRE_datagroup * self,
                                                       const char * path,
                                                       int type,
                                                       SIDRE_SidreLength num_elems);

SIDRE_dataview * SIDRE_datagroup_create_view_from_type_bufferify(
  SIDRE_datagroup * self, const char * path, int Lpath, int type,
  SIDRE_SidreLength num_elems);

SIDRE_dataview * SIDRE_datagroup_create_view_from_type_and_buffer(
  SIDRE_datagroup * self, const char * path, int type,
  SIDRE_SidreLength num_elems, SIDRE_databuffer * buff);

SIDRE_dataview * SIDRE_datagroup_create_view_from_type_and_buffer_bufferify(
  SIDRE_datagroup * self, const char * path, int Lpath, int type,
  SIDRE_SidreLength num_elems, SIDRE_databuffer * buff);

SIDRE_dataview * SIDRE_datagroup_create_view_from_type_external(
  SIDRE_datagroup * self, const char * path, int type,
  SIDRE_SidreLength num_elems, void * external_ptr);

SIDRE_dataview * SIDRE_datagroup_create_view_from_type_external_bufferify(
  SIDRE_datagroup * self, const char * path, int Lpath, int type,
  SIDRE_SidreLength num_elems, void * external_ptr);

SIDRE_dataview * SIDRE_datagroup_create_view_from_shape(SIDRE_datagroup * self,
                                                        const char * path,
                                                        int type, int ndims,
                                                        SIDRE_SidreLength * shape);

SIDRE_dataview * SIDRE_datagroup_create_view_from_shape_bufferify(
  SIDRE_datagroup * self, const char * path, int Lpath, int type, int ndims,
  SIDRE_SidreLength * shape);

SIDRE_dataview * SIDRE_datagroup_create_view_from_shape_and_buffer(
  SIDRE_datagroup * self, const char * path, int type, int ndims,
  SIDRE_SidreLength * shape, SIDRE_databuffer * buff);

SIDRE_dataview * SIDRE_datagroup_create_view_from_shape_and_buffer_bufferify(
  SIDRE_datagroup * self, const char * path, int Lpath, int type, int ndims,
  SIDRE_SidreLength * shape, SIDRE_databuffer * buff);

SIDRE_dataview * SIDRE_datagroup_create_view_from_shape_external(
  SIDRE_datagroup * self, const char * path, int type, int ndims,
  SIDRE_SidreLength * shape, void * external_ptr);

SIDRE_dataview * SIDRE_datagroup_create_view_from_shape_external_bufferify(
  SIDRE_datagroup * self, const char * path, int Lpath, int type, int ndims,
  SIDRE_SidreLength * shape, void * external_ptr);

SIDRE_dataview * SIDRE_datagroup_create_view_into_buffer(SIDRE_datagroup * self,
                                                         const char * path,
                                                         SIDRE_databuffer * buff);

SIDRE_dataview * SIDRE_datagroup_create_view_into_buffer_bufferify(
  SIDRE_datagroup * self, const char * path, int Lpath,
  SIDRE_databuffer * buff);

SIDRE_dataview * SIDRE_datagroup_create_view_external(SIDRE_datagroup * self,
                                                      const char * path,
                                                      void * external_ptr);

SIDRE_dataview * SIDRE_datagroup_create_view_external_bufferify(
  SIDRE_datagroup * self, const char * path, int Lpath, void * external_ptr);

void SIDRE_datagroup_destroy_view(SIDRE_datagroup * self, const char * path);

void SIDRE_datagroup_destroy_view_bufferify(SIDRE_datagroup * self,
                                            const char * path, int Lpath);

void SIDRE_datagroup_destroy_view_and_data_name(SIDRE_datagroup * self,
                                                const char * path);

void SIDRE_datagroup_destroy_view_and_data_name_bufferify(
  SIDRE_datagroup * self, const char * path, int Lpath);

void SIDRE_datagroup_destroy_view_and_data_index(SIDRE_datagroup * self,
                                                 SIDRE_IndexType idx);

SIDRE_dataview * SIDRE_datagroup_move_view(SIDRE_datagroup * self,
                                           SIDRE_dataview * view);

SIDRE_dataview * SIDRE_datagroup_copy_view(SIDRE_datagroup * self,
                                           SIDRE_dataview * view);

bool SIDRE_datagroup_has_group(SIDRE_datagroup * self, const char * path);

bool SIDRE_datagroup_has_group_bufferify(SIDRE_datagroup * self,
                                         const char * path, int Lpath);

bool SIDRE_datagroup_has_child_group(SIDRE_datagroup * self, const char * name);

bool SIDRE_datagroup_has_child_group_bufferify(SIDRE_datagroup * self,
                                               const char * name, int Lname);

SIDRE_datagroup * SIDRE_datagroup_get_group(SIDRE_datagroup * self,
                                            const char * path);

SIDRE_datagroup * SIDRE_datagroup_get_group_bufferify(SIDRE_datagroup * self,
                                                      const char * path,
                                                      int Lpath);

SIDRE_IndexType SIDRE_datagroup_get_group_index(const SIDRE_datagroup * self,
                                                const char * name);

SIDRE_IndexType SIDRE_datagroup_get_group_index_bufferify(
  const SIDRE_datagroup * self, const char * name, int Lname);

const char * SIDRE_datagroup_get_group_name(const SIDRE_datagroup * self,
                                            SIDRE_IndexType idx);

void SIDRE_datagroup_get_group_name_bufferify(const SIDRE_datagroup * self,
                                              SIDRE_IndexType idx,
                                              char * SH_F_rv, int LSH_F_rv);

SIDRE_IndexType SIDRE_datagroup_get_first_valid_group_index(
  const SIDRE_datagroup * self);

SIDRE_IndexType SIDRE_datagroup_get_next_valid_group_index(
  const SIDRE_datagroup * self, SIDRE_IndexType idx);

SIDRE_datagroup * SIDRE_datagroup_create_group(SIDRE_datagroup * self,
                                               const char * path);

SIDRE_datagroup * SIDRE_datagroup_create_group_bufferify(SIDRE_datagroup * self,
                                                         const char * path,
                                                         int Lpath);

void SIDRE_datagroup_destroy_group_name(SIDRE_datagroup * self,
                                        const char * path);

void SIDRE_datagroup_destroy_group_name_bufferify(SIDRE_datagroup * self,
                                                  const char * path, int Lpath);

void SIDRE_datagroup_destroy_group_index(SIDRE_datagroup * self,
                                         SIDRE_IndexType idx);

SIDRE_datagroup * SIDRE_datagroup_move_group(SIDRE_datagroup * self,
                                             SIDRE_datagroup * grp);

void SIDRE_datagroup_print(const SIDRE_datagroup * self);

bool SIDRE_datagroup_is_equivalent_to(const SIDRE_datagroup * self,
                                      const SIDRE_datagroup * other);

void SIDRE_datagroup_save(const SIDRE_datagroup * self, const char * file_path,
                          const char * protocol);

void SIDRE_datagroup_save_bufferify(const SIDRE_datagroup * self,
                                    const char * file_path, int Lfile_path,
                                    const char * protocol, int Lprotocol);

void SIDRE_datagroup_load(SIDRE_datagroup * self, const char * file_path,
                          const char * protocol);

void SIDRE_datagroup_load_bufferify(SIDRE_datagroup * self,
                                    const char * file_path, int Lfile_path,
                                    const char * protocol, int Lprotocol);

void SIDRE_datagroup_load_external_data(SIDRE_datagroup * self,
                                        const char * file_path);

void SIDRE_datagroup_load_external_data_bufferify(SIDRE_datagroup * self,
                                                  const char * file_path,
                                                  int Lfile_path);

#ifdef __cplusplus
}
#endif

#endif  // WRAPDATAGROUP_H
