// wrapIOManager.cpp
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
// wrapIOManager.cpp
#include "wrapIOManager.h"
#include <string>
#include "sidre/DataGroup.hpp"
#include "spio/IOManager.hpp"

extern "C" {
namespace asctoolkit {
namespace spio {

SPIO_iomanager * SPIO_iomanager_new(MPI_Fint com)
{

// splicer begin class.IOManager.method.new
IOManager * rv = new IOManager(MPI_Comm_f2c(com));
return static_cast<SPIO_iomanager *>(static_cast<void *>(rv));
// splicer end class.IOManager.method.new
}

void SPIO_iomanager_delete(SPIO_iomanager * self)
{
IOManager *selfobj = static_cast<IOManager *>(static_cast<void *>(self));
// splicer begin class.IOManager.method.delete
delete selfobj;
// splicer end class.IOManager.method.delete
}

void SPIO_iomanager_write(SPIO_iomanager * self, SIDRE_datagroup * group, int num_files, const char * file_string, const char * protocol)
{
IOManager *selfobj = static_cast<IOManager *>(static_cast<void *>(self));
// splicer begin class.IOManager.method.write
asctoolkit::sidre::DataGroup * SH_group = static_cast<asctoolkit::sidre::DataGroup *>(static_cast<void *>(group));
const std::string SH_file_string(file_string);
const std::string SH_protocol(protocol);
selfobj->write(SH_group, num_files, SH_file_string, SH_protocol);
return;
// splicer end class.IOManager.method.write
}

void SPIO_iomanager_write_bufferify(SPIO_iomanager * self, SIDRE_datagroup * group, int num_files, const char * file_string, int Lfile_string, const char * protocol, int Lprotocol)
{
IOManager *selfobj = static_cast<IOManager *>(static_cast<void *>(self));
// splicer begin class.IOManager.method.write_bufferify
asctoolkit::sidre::DataGroup * SH_group = static_cast<asctoolkit::sidre::DataGroup *>(static_cast<void *>(group));
const std::string SH_file_string(file_string, Lfile_string);
const std::string SH_protocol(protocol, Lprotocol);
selfobj->write(SH_group, num_files, SH_file_string, SH_protocol);
return;
// splicer end class.IOManager.method.write_bufferify
}

void SPIO_iomanager_read_0(SPIO_iomanager * self, SIDRE_datagroup * group, const char * file_string, const char * protocol)
{
IOManager *selfobj = static_cast<IOManager *>(static_cast<void *>(self));
// splicer begin class.IOManager.method.read_0
asctoolkit::sidre::DataGroup * SH_group = static_cast<asctoolkit::sidre::DataGroup *>(static_cast<void *>(group));
const std::string SH_file_string(file_string);
const std::string SH_protocol(protocol);
selfobj->read(SH_group, SH_file_string, SH_protocol);
return;
// splicer end class.IOManager.method.read_0
}

void SPIO_iomanager_read_0_bufferify(SPIO_iomanager * self, SIDRE_datagroup * group, const char * file_string, int Lfile_string, const char * protocol, int Lprotocol)
{
IOManager *selfobj = static_cast<IOManager *>(static_cast<void *>(self));
// splicer begin class.IOManager.method.read_0_bufferify
asctoolkit::sidre::DataGroup * SH_group = static_cast<asctoolkit::sidre::DataGroup *>(static_cast<void *>(group));
const std::string SH_file_string(file_string, Lfile_string);
const std::string SH_protocol(protocol, Lprotocol);
selfobj->read(SH_group, SH_file_string, SH_protocol);
return;
// splicer end class.IOManager.method.read_0_bufferify
}

void SPIO_iomanager_read_1(SPIO_iomanager * self, SIDRE_datagroup * group, const char * root_file)
{
IOManager *selfobj = static_cast<IOManager *>(static_cast<void *>(self));
// splicer begin class.IOManager.method.read_1
asctoolkit::sidre::DataGroup * SH_group = static_cast<asctoolkit::sidre::DataGroup *>(static_cast<void *>(group));
const std::string SH_root_file(root_file);
selfobj->read(SH_group, SH_root_file);
return;
// splicer end class.IOManager.method.read_1
}

void SPIO_iomanager_read_1_bufferify(SPIO_iomanager * self, SIDRE_datagroup * group, const char * root_file, int Lroot_file)
{
IOManager *selfobj = static_cast<IOManager *>(static_cast<void *>(self));
// splicer begin class.IOManager.method.read_1_bufferify
asctoolkit::sidre::DataGroup * SH_group = static_cast<asctoolkit::sidre::DataGroup *>(static_cast<void *>(group));
const std::string SH_root_file(root_file, Lroot_file);
selfobj->read(SH_group, SH_root_file);
return;
// splicer end class.IOManager.method.read_1_bufferify
}

// splicer begin class.IOManager.additional_functions
// splicer end class.IOManager.additional_functions

}  // namespace asctoolkit
}  // namespace spio
}  // extern "C"
