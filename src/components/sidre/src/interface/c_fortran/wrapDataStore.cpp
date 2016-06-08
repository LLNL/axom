// wrapDataStore.cpp
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
// wrapDataStore.cpp
#include "wrapDataStore.h"
#include <string>
#include "sidre/DataStore.hpp"
#include "sidre/SidreTypes.hpp"

extern "C" {
namespace asctoolkit
{
namespace sidre
{

SIDRE_datastore * SIDRE_datastore_new()
{

// splicer begin class.DataStore.method.new
  DataStore * rv = new DataStore();
  return static_cast<SIDRE_datastore *>(static_cast<void *>(rv));
// splicer end class.DataStore.method.new
}

void SIDRE_datastore_delete(SIDRE_datastore * self)
{
  DataStore * selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.delete
  delete selfobj;
// splicer end class.DataStore.method.delete
}

SIDRE_datagroup * SIDRE_datastore_get_root(SIDRE_datastore * self)
{
  DataStore * selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.get_root
  DataGroup * rv = selfobj->getRoot();
  return static_cast<SIDRE_datagroup *>(static_cast<void *>(rv));
// splicer end class.DataStore.method.get_root
}

SIDRE_databuffer * SIDRE_datastore_get_buffer(SIDRE_datastore * self,
                                              SIDRE_IndexType idx)
{
  DataStore * selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.get_buffer
  DataBuffer * rv = selfobj->getBuffer(idx);
  return static_cast<SIDRE_databuffer *>(static_cast<void *>(rv));
// splicer end class.DataStore.method.get_buffer
}

SIDRE_databuffer * SIDRE_datastore_create_buffer_empty(SIDRE_datastore * self)
{
  DataStore * selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.create_buffer_empty
  DataBuffer * rv = selfobj->createBuffer();
  return static_cast<SIDRE_databuffer *>(static_cast<void *>(rv));
// splicer end class.DataStore.method.create_buffer_empty
}

SIDRE_databuffer * SIDRE_datastore_create_buffer_from_type(
  SIDRE_datastore * self, int type, SIDRE_SidreLength num_elems)
{
  DataStore * selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.create_buffer_from_type
  DataBuffer * rv = selfobj->createBuffer(getTypeID(type), num_elems);
  return static_cast<SIDRE_databuffer *>(static_cast<void *>(rv));
// splicer end class.DataStore.method.create_buffer_from_type
}

void SIDRE_datastore_destroy_buffer(SIDRE_datastore * self, SIDRE_IndexType id)
{
  DataStore * selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.destroy_buffer
  selfobj->destroyBuffer(id);
  return;
// splicer end class.DataStore.method.destroy_buffer
}

size_t SIDRE_datastore_get_num_buffers(SIDRE_datastore * self)
{
  DataStore * selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.get_num_buffers
  size_t rv = selfobj->getNumBuffers();
  return rv;
// splicer end class.DataStore.method.get_num_buffers
}

void SIDRE_datastore_print(SIDRE_datastore * self)
{
  DataStore * selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.print
  selfobj->print();
  return;
// splicer end class.DataStore.method.print
}

void SIDRE_datastore_save_0(SIDRE_datastore * self, const char * file_path,
                            const char * protocol)
{
  DataStore * selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.save_0
  const std::string SH_file_path(file_path);
  const std::string SH_protocol(protocol);
  selfobj->save(SH_file_path, SH_protocol);
  return;
// splicer end class.DataStore.method.save_0
}

void SIDRE_datastore_save_0_bufferify(SIDRE_datastore * self,
                                      const char * file_path, int Lfile_path,
                                      const char * protocol, int Lprotocol)
{
  DataStore * selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.save_0_bufferify
  const std::string SH_file_path(file_path, Lfile_path);
  const std::string SH_protocol(protocol, Lprotocol);
  selfobj->save(SH_file_path, SH_protocol);
  return;
// splicer end class.DataStore.method.save_0_bufferify
}

void SIDRE_datastore_save_1(SIDRE_datastore * self, const char * file_path,
                            const char * protocol,
                            const SIDRE_datagroup * group)
{
  DataStore * selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.save_1
  const std::string SH_file_path(file_path);
  const std::string SH_protocol(protocol);
  selfobj->save(SH_file_path, SH_protocol,
                static_cast<const DataGroup *>(static_cast<const void *>(group)));
  return;
// splicer end class.DataStore.method.save_1
}

void SIDRE_datastore_save_1_bufferify(SIDRE_datastore * self,
                                      const char * file_path, int Lfile_path,
                                      const char * protocol, int Lprotocol,
                                      const SIDRE_datagroup * group)
{
  DataStore * selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.save_1_bufferify
  const std::string SH_file_path(file_path, Lfile_path);
  const std::string SH_protocol(protocol, Lprotocol);
  selfobj->save(SH_file_path, SH_protocol,
                static_cast<const DataGroup *>(static_cast<const void *>(group)));
  return;
// splicer end class.DataStore.method.save_1_bufferify
}

void SIDRE_datastore_load_0(SIDRE_datastore * self, const char * file_path,
                            const char * protocol)
{
  DataStore * selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.load_0
  const std::string SH_file_path(file_path);
  const std::string SH_protocol(protocol);
  selfobj->load(SH_file_path, SH_protocol);
  return;
// splicer end class.DataStore.method.load_0
}

void SIDRE_datastore_load_0_bufferify(SIDRE_datastore * self,
                                      const char * file_path, int Lfile_path,
                                      const char * protocol, int Lprotocol)
{
  DataStore * selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.load_0_bufferify
  const std::string SH_file_path(file_path, Lfile_path);
  const std::string SH_protocol(protocol, Lprotocol);
  selfobj->load(SH_file_path, SH_protocol);
  return;
// splicer end class.DataStore.method.load_0_bufferify
}

void SIDRE_datastore_load_1(SIDRE_datastore * self, const char * file_path,
                            const char * protocol, SIDRE_datagroup * group)
{
  DataStore * selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.load_1
  const std::string SH_file_path(file_path);
  const std::string SH_protocol(protocol);
  selfobj->load(SH_file_path, SH_protocol,
                static_cast<DataGroup *>(static_cast<void *>(group)));
  return;
// splicer end class.DataStore.method.load_1
}

void SIDRE_datastore_load_1_bufferify(SIDRE_datastore * self,
                                      const char * file_path, int Lfile_path,
                                      const char * protocol, int Lprotocol,
                                      SIDRE_datagroup * group)
{
  DataStore * selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.load_1_bufferify
  const std::string SH_file_path(file_path, Lfile_path);
  const std::string SH_protocol(protocol, Lprotocol);
  selfobj->load(SH_file_path, SH_protocol,
                static_cast<DataGroup *>(static_cast<void *>(group)));
  return;
// splicer end class.DataStore.method.load_1_bufferify
}

void SIDRE_datastore_load_external_data_0(SIDRE_datastore * self,
                                          const char * file_path,
                                          const char * protocol)
{
  DataStore * selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.load_external_data_0
  const std::string SH_file_path(file_path);
  const std::string SH_protocol(protocol);
  selfobj->loadExternalData(SH_file_path, SH_protocol);
  return;
// splicer end class.DataStore.method.load_external_data_0
}

void SIDRE_datastore_load_external_data_0_bufferify(SIDRE_datastore * self,
                                                    const char * file_path,
                                                    int Lfile_path,
                                                    const char * protocol,
                                                    int Lprotocol)
{
  DataStore * selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.load_external_data_0_bufferify
  const std::string SH_file_path(file_path, Lfile_path);
  const std::string SH_protocol(protocol, Lprotocol);
  selfobj->loadExternalData(SH_file_path, SH_protocol);
  return;
// splicer end class.DataStore.method.load_external_data_0_bufferify
}

void SIDRE_datastore_load_external_data_1(SIDRE_datastore * self,
                                          const char * file_path,
                                          const char * protocol,
                                          SIDRE_datagroup * group)
{
  DataStore * selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.load_external_data_1
  const std::string SH_file_path(file_path);
  const std::string SH_protocol(protocol);
  selfobj->loadExternalData(SH_file_path, SH_protocol,
                            static_cast<DataGroup *>(static_cast<void *>(group)));
  return;
// splicer end class.DataStore.method.load_external_data_1
}

void SIDRE_datastore_load_external_data_1_bufferify(SIDRE_datastore * self,
                                                    const char * file_path,
                                                    int Lfile_path,
                                                    const char * protocol,
                                                    int Lprotocol,
                                                    SIDRE_datagroup * group)
{
  DataStore * selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.load_external_data_1_bufferify
  const std::string SH_file_path(file_path, Lfile_path);
  const std::string SH_protocol(protocol, Lprotocol);
  selfobj->loadExternalData(SH_file_path, SH_protocol,
                            static_cast<DataGroup *>(static_cast<void *>(group)));
  return;
// splicer end class.DataStore.method.load_external_data_1_bufferify
}

// splicer begin class.DataStore.additional_functions
// splicer end class.DataStore.additional_functions

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"
