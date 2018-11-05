// wrapIOManager.cpp
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
#include "wrapIOManager.h"
#include <stdlib.h>
#include <string>
#include "axom/sidre/core/DataStore.hpp"
#include "axom/sidre/core/Group.hpp"
#include "axom/sidre/spio/IOManager.hpp"

// splicer begin class.IOManager.CXX_definitions
// splicer end class.IOManager.CXX_definitions

extern "C" {

// splicer begin class.IOManager.C_definitions
// splicer end class.IOManager.C_definitions

SPIO_iomanager* SPIO_iomanager_new_0(MPI_Fint com, SPIO_iomanager* SHC_rv)
{
// splicer begin class.IOManager.method.new_0
  MPI_Comm SHCXX_com = MPI_Comm_f2c(com);
  axom::sidre::IOManager* SHCXX_rv = new axom::sidre::IOManager(SHCXX_com);
  SHC_rv->addr = static_cast<void*>(SHCXX_rv);
  SHC_rv->idtor = 0;
  return SHC_rv;
// splicer end class.IOManager.method.new_0
}

SPIO_iomanager* SPIO_iomanager_new_1(MPI_Fint com, bool use_scr,
                                     SPIO_iomanager* SHC_rv)
{
// splicer begin class.IOManager.method.new_1
  MPI_Comm SHCXX_com = MPI_Comm_f2c(com);
  axom::sidre::IOManager* SHCXX_rv = new axom::sidre::IOManager(SHCXX_com,
                                                                use_scr);
  SHC_rv->addr = static_cast<void*>(SHCXX_rv);
  SHC_rv->idtor = 0;
  return SHC_rv;
// splicer end class.IOManager.method.new_1
}

void SPIO_iomanager_delete(SPIO_iomanager* self)
{
// splicer begin class.IOManager.method.delete
  axom::sidre::IOManager* SH_this =
    static_cast<axom::sidre::IOManager*>(self->addr);
  delete SH_this;
  self->addr = NULL;
  return;
// splicer end class.IOManager.method.delete
}

void SPIO_iomanager_write(SPIO_iomanager* self, SIDRE_group* group,
                          int num_files, const char* file_string,
                          const char* protocol)
{
// splicer begin class.IOManager.method.write
  axom::sidre::IOManager* SH_this =
    static_cast<axom::sidre::IOManager*>(self->addr);
  axom::sidre::Group* SHCXX_group =
    static_cast<axom::sidre::Group*>(group->addr);
  const std::string SH_file_string(file_string);
  const std::string SH_protocol(protocol);
  SH_this->write(SHCXX_group, num_files, SH_file_string, SH_protocol);
  return;
// splicer end class.IOManager.method.write
}

void SPIO_iomanager_write_bufferify(SPIO_iomanager* self, SIDRE_group* group,
                                    int num_files, const char* file_string,
                                    int Lfile_string, const char* protocol,
                                    int Lprotocol)
{
// splicer begin class.IOManager.method.write_bufferify
  axom::sidre::IOManager* SH_this =
    static_cast<axom::sidre::IOManager*>(self->addr);
  axom::sidre::Group* SHCXX_group =
    static_cast<axom::sidre::Group*>(group->addr);
  const std::string SH_file_string(file_string, Lfile_string);
  const std::string SH_protocol(protocol, Lprotocol);
  SH_this->write(SHCXX_group, num_files, SH_file_string, SH_protocol);
  return;
// splicer end class.IOManager.method.write_bufferify
}

void SPIO_iomanager_write_group_to_root_file(SPIO_iomanager* self,
                                             SIDRE_group* group,
                                             const char* file_name)
{
// splicer begin class.IOManager.method.write_group_to_root_file
  axom::sidre::IOManager* SH_this =
    static_cast<axom::sidre::IOManager*>(self->addr);
  axom::sidre::Group* SHCXX_group =
    static_cast<axom::sidre::Group*>(group->addr);
  const std::string SH_file_name(file_name);
  SH_this->writeGroupToRootFile(SHCXX_group, SH_file_name);
  return;
// splicer end class.IOManager.method.write_group_to_root_file
}

void SPIO_iomanager_write_group_to_root_file_bufferify(SPIO_iomanager* self,
                                                       SIDRE_group* group,
                                                       const char* file_name,
                                                       int Lfile_name)
{
// splicer begin class.IOManager.method.write_group_to_root_file_bufferify
  axom::sidre::IOManager* SH_this =
    static_cast<axom::sidre::IOManager*>(self->addr);
  axom::sidre::Group* SHCXX_group =
    static_cast<axom::sidre::Group*>(group->addr);
  const std::string SH_file_name(file_name, Lfile_name);
  SH_this->writeGroupToRootFile(SHCXX_group, SH_file_name);
  return;
// splicer end class.IOManager.method.write_group_to_root_file_bufferify
}

void SPIO_iomanager_write_blueprint_index_to_root_file(SPIO_iomanager* self,
                                                       SIDRE_datastore* datastore, const char* domain_path, const char* file_name,
                                                       const char* mesh_name)
{
// splicer begin class.IOManager.method.write_blueprint_index_to_root_file
  axom::sidre::IOManager* SH_this =
    static_cast<axom::sidre::IOManager*>(self->addr);
  axom::sidre::DataStore* SHCXX_datastore =
    static_cast<axom::sidre::DataStore*>(datastore->addr);
  const std::string SH_domain_path(domain_path);
  const std::string SH_file_name(file_name);
  const std::string SH_mesh_name(mesh_name);
  SH_this->writeBlueprintIndexToRootFile(SHCXX_datastore, SH_domain_path,
                                         SH_file_name, SH_mesh_name);
  return;
// splicer end class.IOManager.method.write_blueprint_index_to_root_file
}

void SPIO_iomanager_write_blueprint_index_to_root_file_bufferify(
  SPIO_iomanager* self, SIDRE_datastore* datastore, const char* domain_path,
  int Ldomain_path, const char* file_name, int Lfile_name,
  const char* mesh_name, int Lmesh_name)
{
// splicer begin
// class.IOManager.method.write_blueprint_index_to_root_file_bufferify
  axom::sidre::IOManager* SH_this =
    static_cast<axom::sidre::IOManager*>(self->addr);
  axom::sidre::DataStore* SHCXX_datastore =
    static_cast<axom::sidre::DataStore*>(datastore->addr);
  const std::string SH_domain_path(domain_path, Ldomain_path);
  const std::string SH_file_name(file_name, Lfile_name);
  const std::string SH_mesh_name(mesh_name, Lmesh_name);
  SH_this->writeBlueprintIndexToRootFile(SHCXX_datastore, SH_domain_path,
                                         SH_file_name, SH_mesh_name);
  return;
// splicer end
// class.IOManager.method.write_blueprint_index_to_root_file_bufferify
}

void SPIO_iomanager_read_0(SPIO_iomanager* self, SIDRE_group* group,
                           const char* file_string, const char* protocol)
{
// splicer begin class.IOManager.method.read_0
  axom::sidre::IOManager* SH_this =
    static_cast<axom::sidre::IOManager*>(self->addr);
  axom::sidre::Group* SHCXX_group =
    static_cast<axom::sidre::Group*>(group->addr);
  const std::string SH_file_string(file_string);
  const std::string SH_protocol(protocol);
  SH_this->read(SHCXX_group, SH_file_string, SH_protocol);
  return;
// splicer end class.IOManager.method.read_0
}

void SPIO_iomanager_read_0_bufferify(SPIO_iomanager* self, SIDRE_group* group,
                                     const char* file_string, int Lfile_string,
                                     const char* protocol, int Lprotocol)
{
// splicer begin class.IOManager.method.read_0_bufferify
  axom::sidre::IOManager* SH_this =
    static_cast<axom::sidre::IOManager*>(self->addr);
  axom::sidre::Group* SHCXX_group =
    static_cast<axom::sidre::Group*>(group->addr);
  const std::string SH_file_string(file_string, Lfile_string);
  const std::string SH_protocol(protocol, Lprotocol);
  SH_this->read(SHCXX_group, SH_file_string, SH_protocol);
  return;
// splicer end class.IOManager.method.read_0_bufferify
}

void SPIO_iomanager_read_1(SPIO_iomanager* self, SIDRE_group* group,
                           const char* file_string, const char* protocol,
                           bool preserve_contents)
{
// splicer begin class.IOManager.method.read_1
  axom::sidre::IOManager* SH_this =
    static_cast<axom::sidre::IOManager*>(self->addr);
  axom::sidre::Group* SHCXX_group =
    static_cast<axom::sidre::Group*>(group->addr);
  const std::string SH_file_string(file_string);
  const std::string SH_protocol(protocol);
  SH_this->read(SHCXX_group, SH_file_string, SH_protocol, preserve_contents);
  return;
// splicer end class.IOManager.method.read_1
}

void SPIO_iomanager_read_1_bufferify(SPIO_iomanager* self, SIDRE_group* group,
                                     const char* file_string, int Lfile_string,
                                     const char* protocol, int Lprotocol,
                                     bool preserve_contents)
{
// splicer begin class.IOManager.method.read_1_bufferify
  axom::sidre::IOManager* SH_this =
    static_cast<axom::sidre::IOManager*>(self->addr);
  axom::sidre::Group* SHCXX_group =
    static_cast<axom::sidre::Group*>(group->addr);
  const std::string SH_file_string(file_string, Lfile_string);
  const std::string SH_protocol(protocol, Lprotocol);
  SH_this->read(SHCXX_group, SH_file_string, SH_protocol, preserve_contents);
  return;
// splicer end class.IOManager.method.read_1_bufferify
}

void SPIO_iomanager_read_2(SPIO_iomanager* self, SIDRE_group* group,
                           const char* root_file)
{
// splicer begin class.IOManager.method.read_2
  axom::sidre::IOManager* SH_this =
    static_cast<axom::sidre::IOManager*>(self->addr);
  axom::sidre::Group* SHCXX_group =
    static_cast<axom::sidre::Group*>(group->addr);
  const std::string SH_root_file(root_file);
  SH_this->read(SHCXX_group, SH_root_file);
  return;
// splicer end class.IOManager.method.read_2
}

void SPIO_iomanager_read_2_bufferify(SPIO_iomanager* self, SIDRE_group* group,
                                     const char* root_file, int Lroot_file)
{
// splicer begin class.IOManager.method.read_2_bufferify
  axom::sidre::IOManager* SH_this =
    static_cast<axom::sidre::IOManager*>(self->addr);
  axom::sidre::Group* SHCXX_group =
    static_cast<axom::sidre::Group*>(group->addr);
  const std::string SH_root_file(root_file, Lroot_file);
  SH_this->read(SHCXX_group, SH_root_file);
  return;
// splicer end class.IOManager.method.read_2_bufferify
}

void SPIO_iomanager_read_3(SPIO_iomanager* self, SIDRE_group* group,
                           const char* root_file, bool preserve_contents)
{
// splicer begin class.IOManager.method.read_3
  axom::sidre::IOManager* SH_this =
    static_cast<axom::sidre::IOManager*>(self->addr);
  axom::sidre::Group* SHCXX_group =
    static_cast<axom::sidre::Group*>(group->addr);
  const std::string SH_root_file(root_file);
  SH_this->read(SHCXX_group, SH_root_file, preserve_contents);
  return;
// splicer end class.IOManager.method.read_3
}

void SPIO_iomanager_read_3_bufferify(SPIO_iomanager* self, SIDRE_group* group,
                                     const char* root_file, int Lroot_file,
                                     bool preserve_contents)
{
// splicer begin class.IOManager.method.read_3_bufferify
  axom::sidre::IOManager* SH_this =
    static_cast<axom::sidre::IOManager*>(self->addr);
  axom::sidre::Group* SHCXX_group =
    static_cast<axom::sidre::Group*>(group->addr);
  const std::string SH_root_file(root_file, Lroot_file);
  SH_this->read(SHCXX_group, SH_root_file, preserve_contents);
  return;
// splicer end class.IOManager.method.read_3_bufferify
}

void SPIO_iomanager_read_4(SPIO_iomanager* self, SIDRE_group* group,
                           const char* root_file, bool preserve_contents,
                           bool use_scr)
{
// splicer begin class.IOManager.method.read_4
  axom::sidre::IOManager* SH_this =
    static_cast<axom::sidre::IOManager*>(self->addr);
  axom::sidre::Group* SHCXX_group =
    static_cast<axom::sidre::Group*>(group->addr);
  const std::string SH_root_file(root_file);
  SH_this->read(SHCXX_group, SH_root_file, preserve_contents, use_scr);
  return;
// splicer end class.IOManager.method.read_4
}

void SPIO_iomanager_read_4_bufferify(SPIO_iomanager* self, SIDRE_group* group,
                                     const char* root_file, int Lroot_file,
                                     bool preserve_contents, bool use_scr)
{
// splicer begin class.IOManager.method.read_4_bufferify
  axom::sidre::IOManager* SH_this =
    static_cast<axom::sidre::IOManager*>(self->addr);
  axom::sidre::Group* SHCXX_group =
    static_cast<axom::sidre::Group*>(group->addr);
  const std::string SH_root_file(root_file, Lroot_file);
  SH_this->read(SHCXX_group, SH_root_file, preserve_contents, use_scr);
  return;
// splicer end class.IOManager.method.read_4_bufferify
}

void SPIO_iomanager_load_external_data(SPIO_iomanager* self, SIDRE_group* group,
                                       const char* root_file)
{
// splicer begin class.IOManager.method.load_external_data
  axom::sidre::IOManager* SH_this =
    static_cast<axom::sidre::IOManager*>(self->addr);
  axom::sidre::Group* SHCXX_group =
    static_cast<axom::sidre::Group*>(group->addr);
  const std::string SH_root_file(root_file);
  SH_this->loadExternalData(SHCXX_group, SH_root_file);
  return;
// splicer end class.IOManager.method.load_external_data
}

void SPIO_iomanager_load_external_data_bufferify(SPIO_iomanager* self,
                                                 SIDRE_group* group,
                                                 const char* root_file,
                                                 int Lroot_file)
{
// splicer begin class.IOManager.method.load_external_data_bufferify
  axom::sidre::IOManager* SH_this =
    static_cast<axom::sidre::IOManager*>(self->addr);
  axom::sidre::Group* SHCXX_group =
    static_cast<axom::sidre::Group*>(group->addr);
  const std::string SH_root_file(root_file, Lroot_file);
  SH_this->loadExternalData(SHCXX_group, SH_root_file);
  return;
// splicer end class.IOManager.method.load_external_data_bufferify
}

}  // extern "C"
