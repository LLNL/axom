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
#include "sidre/Group.hpp"
#include "spio/IOManager.hpp"

extern "C" {
namespace axom {
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

void SPIO_iomanager_write(SPIO_iomanager * self, SIDRE_group * grp, int num_files, const char * file_string, const char * protocol)
{
IOManager *selfobj = static_cast<IOManager *>(static_cast<void *>(self));
// splicer begin class.IOManager.method.write
    axom::sidre::Group * SH_grp = static_cast<axom::sidre::Group *>(static_cast<void *>(grp));
    const std::string SH_file_string(file_string);
    const std::string SH_protocol(protocol);
    selfobj->write(SH_grp, num_files, SH_file_string, SH_protocol);
    return;
// splicer end class.IOManager.method.write
}

void SPIO_iomanager_write_bufferify(SPIO_iomanager * self, SIDRE_group * grp, int num_files, const char * file_string, int Lfile_string, const char * protocol, int Lprotocol)
{
IOManager *selfobj = static_cast<IOManager *>(static_cast<void *>(self));
// splicer begin class.IOManager.method.write_bufferify
    axom::sidre::Group * SH_grp = static_cast<axom::sidre::Group *>(static_cast<void *>(grp));
    const std::string SH_file_string(file_string, Lfile_string);
    const std::string SH_protocol(protocol, Lprotocol);
    selfobj->write(SH_grp, num_files, SH_file_string, SH_protocol);
    return;
// splicer end class.IOManager.method.write_bufferify
}

void SPIO_iomanager_write_group_to_root_file(SPIO_iomanager * self, SIDRE_group * grp, const char * file_name)
{
IOManager *selfobj = static_cast<IOManager *>(static_cast<void *>(self));
// splicer begin class.IOManager.method.write_group_to_root_file
    axom::sidre::Group * SH_grp = static_cast<axom::sidre::Group *>(static_cast<void *>(grp));
    const std::string SH_file_name(file_name);
    selfobj->writeGroupToRootFile(SH_grp, SH_file_name);
    return;
// splicer end class.IOManager.method.write_group_to_root_file
}

void SPIO_iomanager_write_group_to_root_file_bufferify(SPIO_iomanager * self, SIDRE_group * grp, const char * file_name, int Lfile_name)
{
IOManager *selfobj = static_cast<IOManager *>(static_cast<void *>(self));
// splicer begin class.IOManager.method.write_group_to_root_file_bufferify
    axom::sidre::Group * SH_grp = static_cast<axom::sidre::Group *>(static_cast<void *>(grp));
    const std::string SH_file_name(file_name, Lfile_name);
    selfobj->writeGroupToRootFile(SH_grp, SH_file_name);
    return;
// splicer end class.IOManager.method.write_group_to_root_file_bufferify
}

void SPIO_iomanager_read_0(SPIO_iomanager * self, SIDRE_group * grp, const char * file_string, const char * protocol)
{
IOManager *selfobj = static_cast<IOManager *>(static_cast<void *>(self));
// splicer begin class.IOManager.method.read_0
    axom::sidre::Group * SH_grp = static_cast<axom::sidre::Group *>(static_cast<void *>(grp));
    const std::string SH_file_string(file_string);
    const std::string SH_protocol(protocol);
    selfobj->read(SH_grp, SH_file_string, SH_protocol);
    return;
// splicer end class.IOManager.method.read_0
}

void SPIO_iomanager_read_0_bufferify(SPIO_iomanager * self, SIDRE_group * grp, const char * file_string, int Lfile_string, const char * protocol, int Lprotocol)
{
IOManager *selfobj = static_cast<IOManager *>(static_cast<void *>(self));
// splicer begin class.IOManager.method.read_0_bufferify
    axom::sidre::Group * SH_grp = static_cast<axom::sidre::Group *>(static_cast<void *>(grp));
    const std::string SH_file_string(file_string, Lfile_string);
    const std::string SH_protocol(protocol, Lprotocol);
    selfobj->read(SH_grp, SH_file_string, SH_protocol);
    return;
// splicer end class.IOManager.method.read_0_bufferify
}

void SPIO_iomanager_read_1(SPIO_iomanager * self, SIDRE_group * grp, const char * root_file)
{
IOManager *selfobj = static_cast<IOManager *>(static_cast<void *>(self));
// splicer begin class.IOManager.method.read_1
    axom::sidre::Group * SH_grp = static_cast<axom::sidre::Group *>(static_cast<void *>(grp));
    const std::string SH_root_file(root_file);
    selfobj->read(SH_grp, SH_root_file);
    return;
// splicer end class.IOManager.method.read_1
}

void SPIO_iomanager_read_1_bufferify(SPIO_iomanager * self, SIDRE_group * grp, const char * root_file, int Lroot_file)
{
IOManager *selfobj = static_cast<IOManager *>(static_cast<void *>(self));
// splicer begin class.IOManager.method.read_1_bufferify
    axom::sidre::Group * SH_grp = static_cast<axom::sidre::Group *>(static_cast<void *>(grp));
    const std::string SH_root_file(root_file, Lroot_file);
    selfobj->read(SH_grp, SH_root_file);
    return;
// splicer end class.IOManager.method.read_1_bufferify
}

void SPIO_iomanager_load_external_data(SPIO_iomanager * self, SIDRE_group * grp, const char * root_file)
{
IOManager *selfobj = static_cast<IOManager *>(static_cast<void *>(self));
// splicer begin class.IOManager.method.load_external_data
    axom::sidre::Group * SH_grp = static_cast<axom::sidre::Group *>(static_cast<void *>(grp));
    const std::string SH_root_file(root_file);
    selfobj->loadExternalData(SH_grp, SH_root_file);
    return;
// splicer end class.IOManager.method.load_external_data
}

void SPIO_iomanager_load_external_data_bufferify(SPIO_iomanager * self, SIDRE_group * grp, const char * root_file, int Lroot_file)
{
IOManager *selfobj = static_cast<IOManager *>(static_cast<void *>(self));
// splicer begin class.IOManager.method.load_external_data_bufferify
    axom::sidre::Group * SH_grp = static_cast<axom::sidre::Group *>(static_cast<void *>(grp));
    const std::string SH_root_file(root_file, Lroot_file);
    selfobj->loadExternalData(SH_grp, SH_root_file);
    return;
// splicer end class.IOManager.method.load_external_data_bufferify
}

// splicer begin class.IOManager.additional_functions
// splicer end class.IOManager.additional_functions

}  // namespace axom
}  // namespace spio
}  // extern "C"
