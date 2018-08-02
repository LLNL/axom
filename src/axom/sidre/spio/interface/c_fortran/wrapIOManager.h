// wrapIOManager.h
// This is generated code, do not edit
//
// Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
//
// Produced at the Lawrence Livermore National Laboratory
//
// LLNL-CODE-741217
//
// All rights reserved.
//
// This file is part of Axom.
//
// For details about use and distribution, please read axom/LICENSE.
//
/**
 * \file wrapIOManager.h
 * \brief Shroud generated wrapper for IOManager class
 */
// For C users and C++ implementation

#ifndef WRAPIOMANAGER_H
#define WRAPIOMANAGER_H

#include "mpi.h"
#include "axom/sidre/interface/c_fortran/wrapGroup.h"

// splicer begin class.IOManager.CXX_declarations
// splicer end class.IOManager.CXX_declarations

#ifdef __cplusplus
extern "C" {
#endif

// declaration of shadow types
struct s_SIDRE_group;
typedef struct s_SIDRE_group SIDRE_group;
struct s_SPIO_iomanager;
typedef struct s_SPIO_iomanager SPIO_iomanager;

// splicer begin class.IOManager.C_declarations
// splicer end class.IOManager.C_declarations

SPIO_iomanager* SPIO_iomanager_new_0(MPI_Fint com);

SPIO_iomanager* SPIO_iomanager_new_1(MPI_Fint com, bool use_scr);

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
