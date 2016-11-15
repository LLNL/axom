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

size_t SIDRE_datastore_get_num_buffers(const SIDRE_datastore * self)
{
  const DataStore * selfobj =
    static_cast<const DataStore *>(static_cast<const void *>(self));
// splicer begin class.DataStore.method.get_num_buffers
  size_t rv = selfobj->getNumBuffers();
  return rv;
// splicer end class.DataStore.method.get_num_buffers
}

void SIDRE_datastore_print(const SIDRE_datastore * self)
{
  const DataStore * selfobj =
    static_cast<const DataStore *>(static_cast<const void *>(self));
// splicer begin class.DataStore.method.print
  selfobj->print();
  return;
// splicer end class.DataStore.method.print
}

void SIDRE_datastore_save(const SIDRE_datastore * self, const char * file_path,
                          const char * protocol)
{
  const DataStore * selfobj =
    static_cast<const DataStore *>(static_cast<const void *>(self));
// splicer begin class.DataStore.method.save
  const std::string SH_file_path(file_path);
  const std::string SH_protocol(protocol);
  selfobj->save(SH_file_path, SH_protocol);
  return;
// splicer end class.DataStore.method.save
}

void SIDRE_datastore_save_bufferify(const SIDRE_datastore * self,
                                    const char * file_path, int Lfile_path,
                                    const char * protocol, int Lprotocol)
{
  const DataStore * selfobj =
    static_cast<const DataStore *>(static_cast<const void *>(self));
// splicer begin class.DataStore.method.save_bufferify
  const std::string SH_file_path(file_path, Lfile_path);
  const std::string SH_protocol(protocol, Lprotocol);
  selfobj->save(SH_file_path, SH_protocol);
  return;
// splicer end class.DataStore.method.save_bufferify
}

void SIDRE_datastore_load(SIDRE_datastore * self, const char * file_path,
                          const char * protocol)
{
  DataStore * selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.load
  const std::string SH_file_path(file_path);
  const std::string SH_protocol(protocol);
  selfobj->load(SH_file_path, SH_protocol);
  return;
// splicer end class.DataStore.method.load
}

void SIDRE_datastore_load_bufferify(SIDRE_datastore * self,
                                    const char * file_path, int Lfile_path,
                                    const char * protocol, int Lprotocol)
{
  DataStore * selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.load_bufferify
  const std::string SH_file_path(file_path, Lfile_path);
  const std::string SH_protocol(protocol, Lprotocol);
  selfobj->load(SH_file_path, SH_protocol);
  return;
// splicer end class.DataStore.method.load_bufferify
}

void SIDRE_datastore_load_external_data(SIDRE_datastore * self,
                                        const char * file_path,
                                        const char * protocol)
{
  DataStore * selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.load_external_data
  const std::string SH_file_path(file_path);
  const std::string SH_protocol(protocol);
  selfobj->loadExternalData(SH_file_path, SH_protocol);
  return;
// splicer end class.DataStore.method.load_external_data
}

void SIDRE_datastore_load_external_data_bufferify(SIDRE_datastore * self,
                                                  const char * file_path,
                                                  int Lfile_path,
                                                  const char * protocol,
                                                  int Lprotocol)
{
  DataStore * selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.load_external_data_bufferify
  const std::string SH_file_path(file_path, Lfile_path);
  const std::string SH_protocol(protocol, Lprotocol);
  selfobj->loadExternalData(SH_file_path, SH_protocol);
  return;
// splicer end class.DataStore.method.load_external_data_bufferify
}

// splicer begin class.DataStore.additional_functions
// splicer end class.DataStore.additional_functions

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"
