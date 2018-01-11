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
// wrapIOManager.cpp
#include "wrapIOManager.h"
#include <string>
#include "sidre/Group.hpp"
#include "sidre/IOManager.hpp"

namespace axom
{
namespace sidre
{

// splicer begin class.IOManager.CXX_definitions
// splicer end class.IOManager.CXX_definitions

extern "C" {

// splicer begin class.IOManager.C_definitions
// splicer end class.IOManager.C_definitions

SPIO_iomanager* SPIO_iomanager_new(MPI_Fint com)
{
// splicer begin class.IOManager.method.new
  IOManager* SH_rv = new IOManager(MPI_Comm_f2c(com));
  return static_cast<SPIO_iomanager*>(static_cast<void*>(SH_rv));
// splicer end class.IOManager.method.new
}

void SPIO_iomanager_delete(SPIO_iomanager* self)
{
// splicer begin class.IOManager.method.delete
  IOManager* SH_this = static_cast<IOManager*>(static_cast<void*>(self));
  delete SH_this;
  return;
// splicer end class.IOManager.method.delete
}

void SPIO_iomanager_write(SPIO_iomanager* self, SIDRE_group* group,
                          int num_files, const char* file_string,
                          const char* protocol)
{
// splicer begin class.IOManager.method.write
  IOManager* SH_this = static_cast<IOManager*>(static_cast<void*>(self));
  axom::sidre::Group* SH_group =
    static_cast<axom::sidre::Group*>(static_cast<void*>(group));
  const std::string SH_file_string(file_string);
  const std::string SH_protocol(protocol);
  SH_this->write(SH_group, num_files, SH_file_string, SH_protocol);
  return;
// splicer end class.IOManager.method.write
}

void SPIO_iomanager_write_bufferify(SPIO_iomanager* self, SIDRE_group* group,
                                    int num_files, const char* file_string,
                                    int Lfile_string, const char* protocol,
                                    int Lprotocol)
{
// splicer begin class.IOManager.method.write_bufferify
  IOManager* SH_this = static_cast<IOManager*>(static_cast<void*>(self));
  axom::sidre::Group* SH_group =
    static_cast<axom::sidre::Group*>(static_cast<void*>(group));
  const std::string SH_file_string(file_string, Lfile_string);
  const std::string SH_protocol(protocol, Lprotocol);
  SH_this->write(SH_group, num_files, SH_file_string, SH_protocol);
  return;
// splicer end class.IOManager.method.write_bufferify
}

void SPIO_iomanager_write_group_to_root_file(SPIO_iomanager* self,
                                             SIDRE_group* group,
                                             const char* file_name)
{
// splicer begin class.IOManager.method.write_group_to_root_file
  IOManager* SH_this = static_cast<IOManager*>(static_cast<void*>(self));
  axom::sidre::Group* SH_group =
    static_cast<axom::sidre::Group*>(static_cast<void*>(group));
  const std::string SH_file_name(file_name);
  SH_this->writeGroupToRootFile(SH_group, SH_file_name);
  return;
// splicer end class.IOManager.method.write_group_to_root_file
}

void SPIO_iomanager_write_group_to_root_file_bufferify(SPIO_iomanager* self,
                                                       SIDRE_group* group,
                                                       const char* file_name,
                                                       int Lfile_name)
{
// splicer begin class.IOManager.method.write_group_to_root_file_bufferify
  IOManager* SH_this = static_cast<IOManager*>(static_cast<void*>(self));
  axom::sidre::Group* SH_group =
    static_cast<axom::sidre::Group*>(static_cast<void*>(group));
  const std::string SH_file_name(file_name, Lfile_name);
  SH_this->writeGroupToRootFile(SH_group, SH_file_name);
  return;
// splicer end class.IOManager.method.write_group_to_root_file_bufferify
}

void SPIO_iomanager_read_0(SPIO_iomanager* self, SIDRE_group* group,
                           const char* file_string, const char* protocol)
{
// splicer begin class.IOManager.method.read_0
  IOManager* SH_this = static_cast<IOManager*>(static_cast<void*>(self));
  axom::sidre::Group* SH_group =
    static_cast<axom::sidre::Group*>(static_cast<void*>(group));
  const std::string SH_file_string(file_string);
  const std::string SH_protocol(protocol);
  SH_this->read(SH_group, SH_file_string, SH_protocol);
  return;
// splicer end class.IOManager.method.read_0
}

void SPIO_iomanager_read_0_bufferify(SPIO_iomanager* self, SIDRE_group* group,
                                     const char* file_string, int Lfile_string,
                                     const char* protocol, int Lprotocol)
{
// splicer begin class.IOManager.method.read_0_bufferify
  IOManager* SH_this = static_cast<IOManager*>(static_cast<void*>(self));
  axom::sidre::Group* SH_group =
    static_cast<axom::sidre::Group*>(static_cast<void*>(group));
  const std::string SH_file_string(file_string, Lfile_string);
  const std::string SH_protocol(protocol, Lprotocol);
  SH_this->read(SH_group, SH_file_string, SH_protocol);
  return;
// splicer end class.IOManager.method.read_0_bufferify
}

void SPIO_iomanager_read_1(SPIO_iomanager* self, SIDRE_group* group,
                           const char* root_file)
{
// splicer begin class.IOManager.method.read_1
  IOManager* SH_this = static_cast<IOManager*>(static_cast<void*>(self));
  axom::sidre::Group* SH_group =
    static_cast<axom::sidre::Group*>(static_cast<void*>(group));
  const std::string SH_root_file(root_file);
  SH_this->read(SH_group, SH_root_file);
  return;
// splicer end class.IOManager.method.read_1
}

void SPIO_iomanager_read_1_bufferify(SPIO_iomanager* self, SIDRE_group* group,
                                     const char* root_file, int Lroot_file)
{
// splicer begin class.IOManager.method.read_1_bufferify
  IOManager* SH_this = static_cast<IOManager*>(static_cast<void*>(self));
  axom::sidre::Group* SH_group =
    static_cast<axom::sidre::Group*>(static_cast<void*>(group));
  const std::string SH_root_file(root_file, Lroot_file);
  SH_this->read(SH_group, SH_root_file);
  return;
// splicer end class.IOManager.method.read_1_bufferify
}

void SPIO_iomanager_read_2(SPIO_iomanager* self, SIDRE_group* group,
                           const char* root_file, bool preserve_contents)
{
// splicer begin class.IOManager.method.read_2
  IOManager* SH_this = static_cast<IOManager*>(static_cast<void*>(self));
  axom::sidre::Group* SH_group =
    static_cast<axom::sidre::Group*>(static_cast<void*>(group));
  const std::string SH_root_file(root_file);
  SH_this->read(SH_group, SH_root_file, preserve_contents);
  return;
// splicer end class.IOManager.method.read_2
}

void SPIO_iomanager_read_2_bufferify(SPIO_iomanager* self, SIDRE_group* group,
                                     const char* root_file, int Lroot_file,
                                     bool preserve_contents)
{
// splicer begin class.IOManager.method.read_2_bufferify
  IOManager* SH_this = static_cast<IOManager*>(static_cast<void*>(self));
  axom::sidre::Group* SH_group =
    static_cast<axom::sidre::Group*>(static_cast<void*>(group));
  const std::string SH_root_file(root_file, Lroot_file);
  SH_this->read(SH_group, SH_root_file, preserve_contents);
  return;
// splicer end class.IOManager.method.read_2_bufferify
}

void SPIO_iomanager_load_external_data(SPIO_iomanager* self,
                                       SIDRE_group* group,
                                       const char* root_file)
{
// splicer begin class.IOManager.method.load_external_data
  IOManager* SH_this = static_cast<IOManager*>(static_cast<void*>(self));
  axom::sidre::Group* SH_group =
    static_cast<axom::sidre::Group*>(static_cast<void*>(group));
  const std::string SH_root_file(root_file);
  SH_this->loadExternalData(SH_group, SH_root_file);
  return;
// splicer end class.IOManager.method.load_external_data
}

void SPIO_iomanager_load_external_data_bufferify(SPIO_iomanager* self,
                                                 SIDRE_group* group,
                                                 const char* root_file,
                                                 int Lroot_file)
{
// splicer begin class.IOManager.method.load_external_data_bufferify
  IOManager* SH_this = static_cast<IOManager*>(static_cast<void*>(self));
  axom::sidre::Group* SH_group =
    static_cast<axom::sidre::Group*>(static_cast<void*>(group));
  const std::string SH_root_file(root_file, Lroot_file);
  SH_this->loadExternalData(SH_group, SH_root_file);
  return;
// splicer end class.IOManager.method.load_external_data_bufferify
}

}  // extern "C"

}  // namespace sidre
}  // namespace axom
