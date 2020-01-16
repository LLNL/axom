// wrapIOManager.h
// This is generated code, do not edit
//
// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier (BSD-3-Clause)
//
/**
 * \file wrapIOManager.h
 * \brief Shroud generated wrapper for IOManager class
 */
// For C users and C++ implementation

#ifndef WRAPIOMANAGER_H
#define WRAPIOMANAGER_H

#include "axom/sidre/interface/c_fortran/wrapDataStore.h"
#include "axom/sidre/interface/c_fortran/wrapGroup.h"
#include "mpi.h"
#include "typesSPIO.h"

// splicer begin class.IOManager.CXX_declarations
// splicer end class.IOManager.CXX_declarations

#ifdef __cplusplus
extern "C" {
#endif

// splicer begin class.IOManager.C_declarations
// splicer end class.IOManager.C_declarations

SPIO_iomanager* SPIO_iomanager_new_0(MPI_Fint com, SPIO_iomanager* SHC_rv);

SPIO_iomanager* SPIO_iomanager_new_1(MPI_Fint com, bool use_scr,
                                     SPIO_iomanager* SHC_rv);

void SPIO_iomanager_delete(SPIO_iomanager* self);

void SPIO_iomanager_write(SPIO_iomanager* self, SIDRE_group* group,
                          int num_files, const char* file_string,
                          const char* protocol);

void SPIO_iomanager_write_bufferify(SPIO_iomanager* self, SIDRE_group* group,
                                    int num_files, const char* file_string,
                                    int Lfile_string, const char* protocol,
                                    int Lprotocol);

void SPIO_iomanager_write_group_to_root_file(SPIO_iomanager* self,
                                             SIDRE_group* group,
                                             const char* file_name);

void SPIO_iomanager_write_group_to_root_file_bufferify(SPIO_iomanager* self,
                                                       SIDRE_group* group,
                                                       const char* file_name,
                                                       int Lfile_name);

void SPIO_iomanager_write_blueprint_index_to_root_file(SPIO_iomanager* self,
                                                       SIDRE_datastore* datastore, const char* domain_path, const char* file_name,
                                                       const char* mesh_path);

void SPIO_iomanager_write_blueprint_index_to_root_file_bufferify(
  SPIO_iomanager* self, SIDRE_datastore* datastore, const char* domain_path,
  int Ldomain_path, const char* file_name, int Lfile_name,
  const char* mesh_path, int Lmesh_path);

void SPIO_iomanager_read_0(SPIO_iomanager* self, SIDRE_group* group,
                           const char* file_string, const char* protocol);

void SPIO_iomanager_read_0_bufferify(SPIO_iomanager* self, SIDRE_group* group,
                                     const char* file_string, int Lfile_string,
                                     const char* protocol, int Lprotocol);

void SPIO_iomanager_read_1(SPIO_iomanager* self, SIDRE_group* group,
                           const char* file_string, const char* protocol,
                           bool preserve_contents);

void SPIO_iomanager_read_1_bufferify(SPIO_iomanager* self, SIDRE_group* group,
                                     const char* file_string, int Lfile_string,
                                     const char* protocol, int Lprotocol,
                                     bool preserve_contents);

void SPIO_iomanager_read_2(SPIO_iomanager* self, SIDRE_group* group,
                           const char* root_file);

void SPIO_iomanager_read_2_bufferify(SPIO_iomanager* self, SIDRE_group* group,
                                     const char* root_file, int Lroot_file);

void SPIO_iomanager_read_3(SPIO_iomanager* self, SIDRE_group* group,
                           const char* root_file, bool preserve_contents);

void SPIO_iomanager_read_3_bufferify(SPIO_iomanager* self, SIDRE_group* group,
                                     const char* root_file, int Lroot_file,
                                     bool preserve_contents);

void SPIO_iomanager_read_4(SPIO_iomanager* self, SIDRE_group* group,
                           const char* root_file, bool preserve_contents,
                           bool use_scr);

void SPIO_iomanager_read_4_bufferify(SPIO_iomanager* self, SIDRE_group* group,
                                     const char* root_file, int Lroot_file,
                                     bool preserve_contents, bool use_scr);

void SPIO_iomanager_load_external_data(SPIO_iomanager* self, SIDRE_group* group,
                                       const char* root_file);

void SPIO_iomanager_load_external_data_bufferify(SPIO_iomanager* self,
                                                 SIDRE_group* group,
                                                 const char* root_file,
                                                 int Lroot_file);

#ifdef __cplusplus
}
#endif

#endif  // WRAPIOMANAGER_H
