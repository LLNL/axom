// wrapGroup.cpp
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
#include "wrapGroup.h"
#include <cstring>
#include <string>
#include "sidre/Group.hpp"
#include "sidre/SidreTypes.hpp"

// Copy s into a, blank fill to la characters
// Truncate if a is too short.
static void ShroudStrCopy(char* a, int la, const char* s)
{
  int ls,nm;
  ls = std::strlen(s);
  nm = ls < la ? ls : la;
  std::memcpy(a,s,nm);
  if(la > nm)
    std::memset(a+nm,' ',la-nm);
}

// splicer begin class.Group.CXX_definitions
// splicer end class.Group.CXX_definitions

extern "C" {

// splicer begin class.Group.C_definitions
// splicer end class.Group.C_definitions

SIDRE_IndexType SIDRE_group_get_index(SIDRE_group* self)
{
// splicer begin class.Group.method.get_index
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  axom::sidre::IndexType SHC_rv = SH_this->getIndex();
  return SHC_rv;
// splicer end class.Group.method.get_index
}

const char* SIDRE_group_get_name(const SIDRE_group* self)
{
// splicer begin class.Group.method.get_name
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  const std::string & SHCXX_rv = SH_this->getName();
  const char* SHC_rv = SHCXX_rv.c_str();
  return SHC_rv;
// splicer end class.Group.method.get_name
}

void SIDRE_group_get_name_bufferify(const SIDRE_group* self, char* SHF_rv,
                                    int NSHF_rv)
{
// splicer begin class.Group.method.get_name_bufferify
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  const std::string & SHCXX_rv = SH_this->getName();
  if (SHCXX_rv.empty())
  {
    std::memset(SHF_rv, ' ', NSHF_rv);
  }
  else
  {
    ShroudStrCopy(SHF_rv, NSHF_rv, SHCXX_rv.c_str());
  }
  return;
// splicer end class.Group.method.get_name_bufferify
}

void SIDRE_group_get_path_bufferify(const SIDRE_group* self, char* SHF_rv,
                                    int NSHF_rv)
{
// splicer begin class.Group.method.get_path_bufferify
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  std::string SHCXX_rv = SH_this->getPath();
  if (SHCXX_rv.empty())
  {
    std::memset(SHF_rv, ' ', NSHF_rv);
  }
  else
  {
    ShroudStrCopy(SHF_rv, NSHF_rv, SHCXX_rv.c_str());
  }
  return;
// splicer end class.Group.method.get_path_bufferify
}

void SIDRE_group_get_path_name_bufferify(const SIDRE_group* self, char* SHF_rv,
                                         int NSHF_rv)
{
// splicer begin class.Group.method.get_path_name_bufferify
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  std::string SHCXX_rv = SH_this->getPathName();
  if (SHCXX_rv.empty())
  {
    std::memset(SHF_rv, ' ', NSHF_rv);
  }
  else
  {
    ShroudStrCopy(SHF_rv, NSHF_rv, SHCXX_rv.c_str());
  }
  return;
// splicer end class.Group.method.get_path_name_bufferify
}

const SIDRE_group* SIDRE_group_get_parent(const SIDRE_group* self)
{
// splicer begin class.Group.method.get_parent
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  const axom::sidre::Group* SHCXX_rv = SH_this->getParent();
  const SIDRE_group* SHC_rv =
    static_cast<const SIDRE_group*>(static_cast<const void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.get_parent
}

size_t SIDRE_group_get_num_groups(const SIDRE_group* self)
{
// splicer begin class.Group.method.get_num_groups
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  size_t SHC_rv = SH_this->getNumGroups();
  return SHC_rv;
// splicer end class.Group.method.get_num_groups
}

size_t SIDRE_group_get_num_views(const SIDRE_group* self)
{
// splicer begin class.Group.method.get_num_views
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  size_t SHC_rv = SH_this->getNumViews();
  return SHC_rv;
// splicer end class.Group.method.get_num_views
}

const SIDRE_datastore* SIDRE_group_get_data_store(const SIDRE_group* self)
{
// splicer begin class.Group.method.get_data_store
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  const axom::sidre::DataStore* SHCXX_rv = SH_this->getDataStore();
  const SIDRE_datastore* SHC_rv =
    static_cast<const SIDRE_datastore*>(static_cast<const void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.get_data_store
}

bool SIDRE_group_has_view(const SIDRE_group* self, const char* path)
{
// splicer begin class.Group.method.has_view
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  const std::string SH_path(path);
  bool SHC_rv = SH_this->hasView(SH_path);
  return SHC_rv;
// splicer end class.Group.method.has_view
}

bool SIDRE_group_has_view_bufferify(const SIDRE_group* self, const char* path,
                                    int Lpath)
{
// splicer begin class.Group.method.has_view_bufferify
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  const std::string SH_path(path, Lpath);
  bool SHC_rv = SH_this->hasView(SH_path);
  return SHC_rv;
// splicer end class.Group.method.has_view_bufferify
}

bool SIDRE_group_has_child_view(const SIDRE_group* self, const char* name)
{
// splicer begin class.Group.method.has_child_view
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  const std::string SH_name(name);
  bool SHC_rv = SH_this->hasChildView(SH_name);
  return SHC_rv;
// splicer end class.Group.method.has_child_view
}

bool SIDRE_group_has_child_view_bufferify(const SIDRE_group* self,
                                          const char* name, int Lname)
{
// splicer begin class.Group.method.has_child_view_bufferify
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  const std::string SH_name(name, Lname);
  bool SHC_rv = SH_this->hasChildView(SH_name);
  return SHC_rv;
// splicer end class.Group.method.has_child_view_bufferify
}

SIDRE_IndexType SIDRE_group_get_view_index(const SIDRE_group* self,
                                           const char* name)
{
// splicer begin class.Group.method.get_view_index
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  const std::string SH_name(name);
  axom::sidre::IndexType SHC_rv = SH_this->getViewIndex(SH_name);
  return SHC_rv;
// splicer end class.Group.method.get_view_index
}

SIDRE_IndexType SIDRE_group_get_view_index_bufferify(const SIDRE_group* self,
                                                     const char* name,
                                                     int Lname)
{
// splicer begin class.Group.method.get_view_index_bufferify
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  const std::string SH_name(name, Lname);
  axom::sidre::IndexType SHC_rv = SH_this->getViewIndex(SH_name);
  return SHC_rv;
// splicer end class.Group.method.get_view_index_bufferify
}

const char* SIDRE_group_get_view_name(const SIDRE_group* self,
                                      SIDRE_IndexType idx)
{
// splicer begin class.Group.method.get_view_name
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  const std::string & SHCXX_rv = SH_this->getViewName(idx);
  // C_error_pattern
  if (!axom::sidre::nameIsValid(SHCXX_rv))
  {
    return SIDRE_InvalidName;
  }

  const char* SHC_rv = SHCXX_rv.c_str();
  return SHC_rv;
// splicer end class.Group.method.get_view_name
}

void SIDRE_group_get_view_name_bufferify(const SIDRE_group* self,
                                         SIDRE_IndexType idx, char* SHF_rv,
                                         int NSHF_rv)
{
// splicer begin class.Group.method.get_view_name_bufferify
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  const std::string & SHCXX_rv = SH_this->getViewName(idx);
  // C_error_pattern
  if (!axom::sidre::nameIsValid(SHCXX_rv))
  {
    std::memset(SHF_rv, ' ', NSHF_rv);
    return;
  }

  if (SHCXX_rv.empty())
  {
    std::memset(SHF_rv, ' ', NSHF_rv);
  }
  else
  {
    ShroudStrCopy(SHF_rv, NSHF_rv, SHCXX_rv.c_str());
  }
  return;
// splicer end class.Group.method.get_view_name_bufferify
}

SIDRE_view* SIDRE_group_get_view_from_name(SIDRE_group* self, const char* path)
{
// splicer begin class.Group.method.get_view_from_name
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path);
  axom::sidre::View* SHCXX_rv = SH_this->getView(SH_path);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.get_view_from_name
}

SIDRE_view* SIDRE_group_get_view_from_name_bufferify(SIDRE_group* self,
                                                     const char* path,
                                                     int Lpath)
{
// splicer begin class.Group.method.get_view_from_name_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path, Lpath);
  axom::sidre::View* SHCXX_rv = SH_this->getView(SH_path);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.get_view_from_name_bufferify
}

SIDRE_view* SIDRE_group_get_view_from_index(SIDRE_group* self,
                                            const SIDRE_IndexType idx)
{
// splicer begin class.Group.method.get_view_from_index
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  axom::sidre::View* SHCXX_rv = SH_this->getView(idx);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.get_view_from_index
}

SIDRE_IndexType SIDRE_group_get_first_valid_view_index(const SIDRE_group* self)
{
// splicer begin class.Group.method.get_first_valid_view_index
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  axom::sidre::IndexType SHC_rv = SH_this->getFirstValidViewIndex();
  return SHC_rv;
// splicer end class.Group.method.get_first_valid_view_index
}

SIDRE_IndexType SIDRE_group_get_next_valid_view_index(const SIDRE_group* self,
                                                      SIDRE_IndexType idx)
{
// splicer begin class.Group.method.get_next_valid_view_index
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  axom::sidre::IndexType SHC_rv = SH_this->getNextValidViewIndex(idx);
  return SHC_rv;
// splicer end class.Group.method.get_next_valid_view_index
}

SIDRE_view* SIDRE_group_create_view_empty(SIDRE_group* self, const char* path)
{
// splicer begin class.Group.method.create_view_empty
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path);
  axom::sidre::View* SHCXX_rv = SH_this->createView(SH_path);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_empty
}

SIDRE_view* SIDRE_group_create_view_empty_bufferify(SIDRE_group* self,
                                                    const char* path, int Lpath)
{
// splicer begin class.Group.method.create_view_empty_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path, Lpath);
  axom::sidre::View* SHCXX_rv = SH_this->createView(SH_path);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_empty_bufferify
}

SIDRE_view* SIDRE_group_create_view_from_type(SIDRE_group* self,
                                              const char* path, int type,
                                              SIDRE_SidreLength num_elems)
{
// splicer begin class.Group.method.create_view_from_type
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  axom::sidre::View* SHCXX_rv = SH_this->createView(SH_path, SHCXX_type,
                                                    num_elems);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_from_type
}

SIDRE_view* SIDRE_group_create_view_from_type_bufferify(SIDRE_group* self,
                                                        const char* path,
                                                        int Lpath, int type,
                                                        SIDRE_SidreLength num_elems)
{
// splicer begin class.Group.method.create_view_from_type_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path, Lpath);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  axom::sidre::View* SHCXX_rv = SH_this->createView(SH_path, SHCXX_type,
                                                    num_elems);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_from_type_bufferify
}

SIDRE_view* SIDRE_group_create_view_from_shape(SIDRE_group* self,
                                               const char* path, int type,
                                               int ndims,
                                               SIDRE_SidreLength* shape)
{
// splicer begin class.Group.method.create_view_from_shape
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  axom::sidre::View* SHCXX_rv = SH_this->createView(SH_path, SHCXX_type, ndims,
                                                    shape);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_from_shape
}

SIDRE_view* SIDRE_group_create_view_from_shape_bufferify(SIDRE_group* self,
                                                         const char* path,
                                                         int Lpath, int type,
                                                         int ndims,
                                                         SIDRE_SidreLength* shape)
{
// splicer begin class.Group.method.create_view_from_shape_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path, Lpath);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  axom::sidre::View* SHCXX_rv = SH_this->createView(SH_path, SHCXX_type, ndims,
                                                    shape);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_from_shape_bufferify
}

SIDRE_view* SIDRE_group_create_view_into_buffer(SIDRE_group* self,
                                                const char* path,
                                                SIDRE_buffer* buff)
{
// splicer begin class.Group.method.create_view_into_buffer
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path);
  axom::sidre::Buffer* SHCXX_buff =
    static_cast<axom::sidre::Buffer*>(static_cast<void*>(buff));
  axom::sidre::View* SHCXX_rv = SH_this->createView(SH_path, SHCXX_buff);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_into_buffer
}

SIDRE_view* SIDRE_group_create_view_into_buffer_bufferify(SIDRE_group* self,
                                                          const char* path,
                                                          int Lpath,
                                                          SIDRE_buffer* buff)
{
// splicer begin class.Group.method.create_view_into_buffer_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path, Lpath);
  axom::sidre::Buffer* SHCXX_buff =
    static_cast<axom::sidre::Buffer*>(static_cast<void*>(buff));
  axom::sidre::View* SHCXX_rv = SH_this->createView(SH_path, SHCXX_buff);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_into_buffer_bufferify
}

SIDRE_view* SIDRE_group_create_view_from_type_and_buffer(SIDRE_group* self,
                                                         const char* path,
                                                         int type,
                                                         SIDRE_SidreLength num_elems,
                                                         SIDRE_buffer* buff)
{
// splicer begin class.Group.method.create_view_from_type_and_buffer
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  axom::sidre::Buffer* SHCXX_buff =
    static_cast<axom::sidre::Buffer*>(static_cast<void*>(buff));
  axom::sidre::View* SHCXX_rv = SH_this->createView(SH_path, SHCXX_type,
                                                    num_elems, SHCXX_buff);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_from_type_and_buffer
}

SIDRE_view* SIDRE_group_create_view_from_type_and_buffer_bufferify(
  SIDRE_group* self, const char* path, int Lpath, int type,
  SIDRE_SidreLength num_elems, SIDRE_buffer* buff)
{
// splicer begin class.Group.method.create_view_from_type_and_buffer_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path, Lpath);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  axom::sidre::Buffer* SHCXX_buff =
    static_cast<axom::sidre::Buffer*>(static_cast<void*>(buff));
  axom::sidre::View* SHCXX_rv = SH_this->createView(SH_path, SHCXX_type,
                                                    num_elems, SHCXX_buff);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_from_type_and_buffer_bufferify
}

SIDRE_view* SIDRE_group_create_view_from_shape_and_buffer(SIDRE_group* self,
                                                          const char* path,
                                                          int type, int ndims,
                                                          SIDRE_SidreLength* shape,
                                                          SIDRE_buffer* buff)
{
// splicer begin class.Group.method.create_view_from_shape_and_buffer
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  axom::sidre::Buffer* SHCXX_buff =
    static_cast<axom::sidre::Buffer*>(static_cast<void*>(buff));
  axom::sidre::View* SHCXX_rv = SH_this->createView(SH_path, SHCXX_type, ndims,
                                                    shape, SHCXX_buff);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_from_shape_and_buffer
}

SIDRE_view* SIDRE_group_create_view_from_shape_and_buffer_bufferify(
  SIDRE_group* self, const char* path, int Lpath, int type, int ndims,
  SIDRE_SidreLength* shape, SIDRE_buffer* buff)
{
// splicer begin class.Group.method.create_view_from_shape_and_buffer_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path, Lpath);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  axom::sidre::Buffer* SHCXX_buff =
    static_cast<axom::sidre::Buffer*>(static_cast<void*>(buff));
  axom::sidre::View* SHCXX_rv = SH_this->createView(SH_path, SHCXX_type, ndims,
                                                    shape, SHCXX_buff);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_from_shape_and_buffer_bufferify
}

SIDRE_view* SIDRE_group_create_view_external(SIDRE_group* self,
                                             const char* path,
                                             void* external_ptr)
{
// splicer begin class.Group.method.create_view_external
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path);
  axom::sidre::View* SHCXX_rv = SH_this->createView(SH_path, external_ptr);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_external
}

SIDRE_view* SIDRE_group_create_view_external_bufferify(SIDRE_group* self,
                                                       const char* path,
                                                       int Lpath,
                                                       void* external_ptr)
{
// splicer begin class.Group.method.create_view_external_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path, Lpath);
  axom::sidre::View* SHCXX_rv = SH_this->createView(SH_path, external_ptr);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_external_bufferify
}

SIDRE_view* SIDRE_group_create_view_from_type_external(SIDRE_group* self,
                                                       const char* path,
                                                       int type,
                                                       SIDRE_SidreLength num_elems,
                                                       void* external_ptr)
{
// splicer begin class.Group.method.create_view_from_type_external
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  axom::sidre::View* SHCXX_rv = SH_this->createView(SH_path, SHCXX_type,
                                                    num_elems, external_ptr);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_from_type_external
}

SIDRE_view* SIDRE_group_create_view_from_type_external_bufferify(
  SIDRE_group* self, const char* path, int Lpath, int type,
  SIDRE_SidreLength num_elems, void* external_ptr)
{
// splicer begin class.Group.method.create_view_from_type_external_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path, Lpath);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  axom::sidre::View* SHCXX_rv = SH_this->createView(SH_path, SHCXX_type,
                                                    num_elems, external_ptr);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_from_type_external_bufferify
}

SIDRE_view* SIDRE_group_create_view_from_shape_external(SIDRE_group* self,
                                                        const char* path,
                                                        int type, int ndims,
                                                        SIDRE_SidreLength* shape,
                                                        void* external_ptr)
{
// splicer begin class.Group.method.create_view_from_shape_external
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  axom::sidre::View* SHCXX_rv = SH_this->createView(SH_path, SHCXX_type, ndims,
                                                    shape, external_ptr);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_from_shape_external
}

SIDRE_view* SIDRE_group_create_view_from_shape_external_bufferify(
  SIDRE_group* self, const char* path, int Lpath, int type, int ndims,
  SIDRE_SidreLength* shape, void* external_ptr)
{
// splicer begin class.Group.method.create_view_from_shape_external_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path, Lpath);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  axom::sidre::View* SHCXX_rv = SH_this->createView(SH_path, SHCXX_type, ndims,
                                                    shape, external_ptr);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_from_shape_external_bufferify
}

SIDRE_view* SIDRE_group_create_view_and_allocate_nelems(SIDRE_group* self,
                                                        const char* path,
                                                        int type,
                                                        SIDRE_SidreLength num_elems)
{
// splicer begin class.Group.method.create_view_and_allocate_nelems
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  axom::sidre::View* SHCXX_rv = SH_this->createViewAndAllocate(SH_path,
                                                               SHCXX_type,
                                                               num_elems);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_and_allocate_nelems
}

SIDRE_view* SIDRE_group_create_view_and_allocate_nelems_bufferify(
  SIDRE_group* self, const char* path, int Lpath, int type,
  SIDRE_SidreLength num_elems)
{
// splicer begin class.Group.method.create_view_and_allocate_nelems_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path, Lpath);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  axom::sidre::View* SHCXX_rv = SH_this->createViewAndAllocate(SH_path,
                                                               SHCXX_type,
                                                               num_elems);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_and_allocate_nelems_bufferify
}

SIDRE_view* SIDRE_group_create_view_and_allocate_shape(SIDRE_group* self,
                                                       const char* path,
                                                       int type, int ndims,
                                                       SIDRE_SidreLength* shape)
{
// splicer begin class.Group.method.create_view_and_allocate_shape
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  axom::sidre::View* SHCXX_rv = SH_this->createViewAndAllocate(SH_path,
                                                               SHCXX_type,
                                                               ndims, shape);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_and_allocate_shape
}

SIDRE_view* SIDRE_group_create_view_and_allocate_shape_bufferify(
  SIDRE_group* self, const char* path, int Lpath, int type, int ndims,
  SIDRE_SidreLength* shape)
{
// splicer begin class.Group.method.create_view_and_allocate_shape_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path, Lpath);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  axom::sidre::View* SHCXX_rv = SH_this->createViewAndAllocate(SH_path,
                                                               SHCXX_type,
                                                               ndims, shape);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_and_allocate_shape_bufferify
}

SIDRE_view* SIDRE_group_create_view_scalar_int(SIDRE_group* self,
                                               const char* path, int value)
{
// splicer begin class.Group.method.create_view_scalar_int
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path);
  axom::sidre::View* SHCXX_rv = SH_this->createViewScalar<int>(SH_path, value);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_scalar_int
}

SIDRE_view* SIDRE_group_create_view_scalar_int_bufferify(SIDRE_group* self,
                                                         const char* path,
                                                         int Lpath, int value)
{
// splicer begin class.Group.method.create_view_scalar_int_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path, Lpath);
  axom::sidre::View* SHCXX_rv = SH_this->createViewScalar<int>(SH_path, value);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_scalar_int_bufferify
}

SIDRE_view* SIDRE_group_create_view_scalar_long(SIDRE_group* self,
                                                const char* path, long value)
{
// splicer begin class.Group.method.create_view_scalar_long
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path);
  axom::sidre::View* SHCXX_rv = SH_this->createViewScalar<long>(SH_path, value);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_scalar_long
}

SIDRE_view* SIDRE_group_create_view_scalar_long_bufferify(SIDRE_group* self,
                                                          const char* path,
                                                          int Lpath, long value)
{
// splicer begin class.Group.method.create_view_scalar_long_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path, Lpath);
  axom::sidre::View* SHCXX_rv = SH_this->createViewScalar<long>(SH_path, value);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_scalar_long_bufferify
}

SIDRE_view* SIDRE_group_create_view_scalar_float(SIDRE_group* self,
                                                 const char* path, float value)
{
// splicer begin class.Group.method.create_view_scalar_float
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path);
  axom::sidre::View* SHCXX_rv =
    SH_this->createViewScalar<float>(SH_path, value);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_scalar_float
}

SIDRE_view* SIDRE_group_create_view_scalar_float_bufferify(SIDRE_group* self,
                                                           const char* path,
                                                           int Lpath,
                                                           float value)
{
// splicer begin class.Group.method.create_view_scalar_float_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path, Lpath);
  axom::sidre::View* SHCXX_rv =
    SH_this->createViewScalar<float>(SH_path, value);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_scalar_float_bufferify
}

SIDRE_view* SIDRE_group_create_view_scalar_double(SIDRE_group* self,
                                                  const char* path,
                                                  double value)
{
// splicer begin class.Group.method.create_view_scalar_double
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path);
  axom::sidre::View* SHCXX_rv =
    SH_this->createViewScalar<double>(SH_path, value);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_scalar_double
}

SIDRE_view* SIDRE_group_create_view_scalar_double_bufferify(SIDRE_group* self,
                                                            const char* path,
                                                            int Lpath,
                                                            double value)
{
// splicer begin class.Group.method.create_view_scalar_double_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path, Lpath);
  axom::sidre::View* SHCXX_rv =
    SH_this->createViewScalar<double>(SH_path, value);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_scalar_double_bufferify
}

SIDRE_view* SIDRE_group_create_view_string(SIDRE_group* self, const char* path,
                                           const char* value)
{
// splicer begin class.Group.method.create_view_string
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path);
  const std::string SH_value(value);
  axom::sidre::View* SHCXX_rv = SH_this->createViewString(SH_path, SH_value);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_string
}

SIDRE_view* SIDRE_group_create_view_string_bufferify(SIDRE_group* self,
                                                     const char* path,
                                                     int Lpath,
                                                     const char* value,
                                                     int Lvalue)
{
// splicer begin class.Group.method.create_view_string_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path, Lpath);
  const std::string SH_value(value, Lvalue);
  axom::sidre::View* SHCXX_rv = SH_this->createViewString(SH_path, SH_value);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_view_string_bufferify
}

void SIDRE_group_destroy_view(SIDRE_group* self, const char* path)
{
// splicer begin class.Group.method.destroy_view
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path);
  SH_this->destroyView(SH_path);
  return;
// splicer end class.Group.method.destroy_view
}

void SIDRE_group_destroy_view_bufferify(SIDRE_group* self, const char* path,
                                        int Lpath)
{
// splicer begin class.Group.method.destroy_view_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path, Lpath);
  SH_this->destroyView(SH_path);
  return;
// splicer end class.Group.method.destroy_view_bufferify
}

void SIDRE_group_destroy_view_and_data_name(SIDRE_group* self, const char* path)
{
// splicer begin class.Group.method.destroy_view_and_data_name
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path);
  SH_this->destroyViewAndData(SH_path);
  return;
// splicer end class.Group.method.destroy_view_and_data_name
}

void SIDRE_group_destroy_view_and_data_name_bufferify(SIDRE_group* self,
                                                      const char* path,
                                                      int Lpath)
{
// splicer begin class.Group.method.destroy_view_and_data_name_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path, Lpath);
  SH_this->destroyViewAndData(SH_path);
  return;
// splicer end class.Group.method.destroy_view_and_data_name_bufferify
}

void SIDRE_group_destroy_view_and_data_index(SIDRE_group* self,
                                             SIDRE_IndexType idx)
{
// splicer begin class.Group.method.destroy_view_and_data_index
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  SH_this->destroyViewAndData(idx);
  return;
// splicer end class.Group.method.destroy_view_and_data_index
}

SIDRE_view* SIDRE_group_move_view(SIDRE_group* self, SIDRE_view* view)
{
// splicer begin class.Group.method.move_view
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  axom::sidre::View* SHCXX_view =
    static_cast<axom::sidre::View*>(static_cast<void*>(view));
  axom::sidre::View* SHCXX_rv = SH_this->moveView(SHCXX_view);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.move_view
}

SIDRE_view* SIDRE_group_copy_view(SIDRE_group* self, SIDRE_view* view)
{
// splicer begin class.Group.method.copy_view
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  axom::sidre::View* SHCXX_view =
    static_cast<axom::sidre::View*>(static_cast<void*>(view));
  axom::sidre::View* SHCXX_rv = SH_this->copyView(SHCXX_view);
  SIDRE_view* SHC_rv = static_cast<SIDRE_view*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.copy_view
}

bool SIDRE_group_has_group(SIDRE_group* self, const char* path)
{
// splicer begin class.Group.method.has_group
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path);
  bool SHC_rv = SH_this->hasGroup(SH_path);
  return SHC_rv;
// splicer end class.Group.method.has_group
}

bool SIDRE_group_has_group_bufferify(SIDRE_group* self, const char* path,
                                     int Lpath)
{
// splicer begin class.Group.method.has_group_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path, Lpath);
  bool SHC_rv = SH_this->hasGroup(SH_path);
  return SHC_rv;
// splicer end class.Group.method.has_group_bufferify
}

bool SIDRE_group_has_child_group(SIDRE_group* self, const char* name)
{
// splicer begin class.Group.method.has_child_group
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_name(name);
  bool SHC_rv = SH_this->hasChildGroup(SH_name);
  return SHC_rv;
// splicer end class.Group.method.has_child_group
}

bool SIDRE_group_has_child_group_bufferify(SIDRE_group* self, const char* name,
                                           int Lname)
{
// splicer begin class.Group.method.has_child_group_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_name(name, Lname);
  bool SHC_rv = SH_this->hasChildGroup(SH_name);
  return SHC_rv;
// splicer end class.Group.method.has_child_group_bufferify
}

SIDRE_IndexType SIDRE_group_get_group_index(const SIDRE_group* self,
                                            const char* name)
{
// splicer begin class.Group.method.get_group_index
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  const std::string SH_name(name);
  axom::sidre::IndexType SHC_rv = SH_this->getGroupIndex(SH_name);
  return SHC_rv;
// splicer end class.Group.method.get_group_index
}

SIDRE_IndexType SIDRE_group_get_group_index_bufferify(const SIDRE_group* self,
                                                      const char* name,
                                                      int Lname)
{
// splicer begin class.Group.method.get_group_index_bufferify
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  const std::string SH_name(name, Lname);
  axom::sidre::IndexType SHC_rv = SH_this->getGroupIndex(SH_name);
  return SHC_rv;
// splicer end class.Group.method.get_group_index_bufferify
}

const char* SIDRE_group_get_group_name(const SIDRE_group* self,
                                       SIDRE_IndexType idx)
{
// splicer begin class.Group.method.get_group_name
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  const std::string & SHCXX_rv = SH_this->getGroupName(idx);
  // C_error_pattern
  if (!axom::sidre::nameIsValid(SHCXX_rv))
  {
    return SIDRE_InvalidName;
  }

  const char* SHC_rv = SHCXX_rv.c_str();
  return SHC_rv;
// splicer end class.Group.method.get_group_name
}

void SIDRE_group_get_group_name_bufferify(const SIDRE_group* self,
                                          SIDRE_IndexType idx, char* SHF_rv,
                                          int NSHF_rv)
{
// splicer begin class.Group.method.get_group_name_bufferify
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  const std::string & SHCXX_rv = SH_this->getGroupName(idx);
  // C_error_pattern
  if (!axom::sidre::nameIsValid(SHCXX_rv))
  {
    std::memset(SHF_rv, ' ', NSHF_rv);
    return;
  }

  if (SHCXX_rv.empty())
  {
    std::memset(SHF_rv, ' ', NSHF_rv);
  }
  else
  {
    ShroudStrCopy(SHF_rv, NSHF_rv, SHCXX_rv.c_str());
  }
  return;
// splicer end class.Group.method.get_group_name_bufferify
}

SIDRE_group* SIDRE_group_get_group_from_name(SIDRE_group* self,
                                             const char* path)
{
// splicer begin class.Group.method.get_group_from_name
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path);
  axom::sidre::Group* SHCXX_rv = SH_this->getGroup(SH_path);
  SIDRE_group* SHC_rv = static_cast<SIDRE_group*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.get_group_from_name
}

SIDRE_group* SIDRE_group_get_group_from_name_bufferify(SIDRE_group* self,
                                                       const char* path,
                                                       int Lpath)
{
// splicer begin class.Group.method.get_group_from_name_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path, Lpath);
  axom::sidre::Group* SHCXX_rv = SH_this->getGroup(SH_path);
  SIDRE_group* SHC_rv = static_cast<SIDRE_group*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.get_group_from_name_bufferify
}

SIDRE_group* SIDRE_group_get_group_from_index(SIDRE_group* self,
                                              SIDRE_IndexType idx)
{
// splicer begin class.Group.method.get_group_from_index
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  axom::sidre::Group* SHCXX_rv = SH_this->getGroup(idx);
  SIDRE_group* SHC_rv = static_cast<SIDRE_group*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.get_group_from_index
}

SIDRE_IndexType SIDRE_group_get_first_valid_group_index(const SIDRE_group* self)
{
// splicer begin class.Group.method.get_first_valid_group_index
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  axom::sidre::IndexType SHC_rv = SH_this->getFirstValidGroupIndex();
  return SHC_rv;
// splicer end class.Group.method.get_first_valid_group_index
}

SIDRE_IndexType SIDRE_group_get_next_valid_group_index(const SIDRE_group* self,
                                                       SIDRE_IndexType idx)
{
// splicer begin class.Group.method.get_next_valid_group_index
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  axom::sidre::IndexType SHC_rv = SH_this->getNextValidGroupIndex(idx);
  return SHC_rv;
// splicer end class.Group.method.get_next_valid_group_index
}

SIDRE_group* SIDRE_group_create_group(SIDRE_group* self, const char* path)
{
// splicer begin class.Group.method.create_group
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path);
  axom::sidre::Group* SHCXX_rv = SH_this->createGroup(SH_path);
  SIDRE_group* SHC_rv = static_cast<SIDRE_group*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_group
}

SIDRE_group* SIDRE_group_create_group_bufferify(SIDRE_group* self,
                                                const char* path, int Lpath)
{
// splicer begin class.Group.method.create_group_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path, Lpath);
  axom::sidre::Group* SHCXX_rv = SH_this->createGroup(SH_path);
  SIDRE_group* SHC_rv = static_cast<SIDRE_group*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.create_group_bufferify
}

void SIDRE_group_destroy_group_name(SIDRE_group* self, const char* path)
{
// splicer begin class.Group.method.destroy_group_name
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path);
  SH_this->destroyGroup(SH_path);
  return;
// splicer end class.Group.method.destroy_group_name
}

void SIDRE_group_destroy_group_name_bufferify(SIDRE_group* self,
                                              const char* path, int Lpath)
{
// splicer begin class.Group.method.destroy_group_name_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_path(path, Lpath);
  SH_this->destroyGroup(SH_path);
  return;
// splicer end class.Group.method.destroy_group_name_bufferify
}

void SIDRE_group_destroy_group_index(SIDRE_group* self, SIDRE_IndexType idx)
{
// splicer begin class.Group.method.destroy_group_index
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  SH_this->destroyGroup(idx);
  return;
// splicer end class.Group.method.destroy_group_index
}

SIDRE_group* SIDRE_group_move_group(SIDRE_group* self, SIDRE_group* grp)
{
// splicer begin class.Group.method.move_group
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  axom::sidre::Group* SHCXX_grp =
    static_cast<axom::sidre::Group*>(static_cast<void*>(grp));
  axom::sidre::Group* SHCXX_rv = SH_this->moveGroup(SHCXX_grp);
  SIDRE_group* SHC_rv = static_cast<SIDRE_group*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.Group.method.move_group
}

void SIDRE_group_print(const SIDRE_group* self)
{
// splicer begin class.Group.method.print
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  SH_this->print();
  return;
// splicer end class.Group.method.print
}

bool SIDRE_group_is_equivalent_to(const SIDRE_group* self,
                                  const SIDRE_group* other)
{
// splicer begin class.Group.method.is_equivalent_to
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  const axom::sidre::Group* SHCXX_other =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(other));
  bool SHC_rv = SH_this->isEquivalentTo(SHCXX_other);
  return SHC_rv;
// splicer end class.Group.method.is_equivalent_to
}

void SIDRE_group_save(const SIDRE_group* self, const char* file_path,
                      const char* protocol)
{
// splicer begin class.Group.method.save
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  const std::string SH_file_path(file_path);
  const std::string SH_protocol(protocol);
  SH_this->save(SH_file_path, SH_protocol);
  return;
// splicer end class.Group.method.save
}

void SIDRE_group_save_bufferify(const SIDRE_group* self, const char* file_path,
                                int Lfile_path, const char* protocol,
                                int Lprotocol)
{
// splicer begin class.Group.method.save_bufferify
  const axom::sidre::Group* SH_this =
    static_cast<const axom::sidre::Group*>(static_cast<const void*>(self));
  const std::string SH_file_path(file_path, Lfile_path);
  const std::string SH_protocol(protocol, Lprotocol);
  SH_this->save(SH_file_path, SH_protocol);
  return;
// splicer end class.Group.method.save_bufferify
}

void SIDRE_group_load_0(SIDRE_group* self, const char* file_path,
                        const char* protocol)
{
// splicer begin class.Group.method.load_0
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_file_path(file_path);
  const std::string SH_protocol(protocol);
  SH_this->load(SH_file_path, SH_protocol);
  return;
// splicer end class.Group.method.load_0
}

void SIDRE_group_load_0_bufferify(SIDRE_group* self, const char* file_path,
                                  int Lfile_path, const char* protocol,
                                  int Lprotocol)
{
// splicer begin class.Group.method.load_0_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_file_path(file_path, Lfile_path);
  const std::string SH_protocol(protocol, Lprotocol);
  SH_this->load(SH_file_path, SH_protocol);
  return;
// splicer end class.Group.method.load_0_bufferify
}

void SIDRE_group_load_1(SIDRE_group* self, const char* file_path,
                        const char* protocol, bool preserve_contents)
{
// splicer begin class.Group.method.load_1
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_file_path(file_path);
  const std::string SH_protocol(protocol);
  SH_this->load(SH_file_path, SH_protocol, preserve_contents);
  return;
// splicer end class.Group.method.load_1
}

void SIDRE_group_load_1_bufferify(SIDRE_group* self, const char* file_path,
                                  int Lfile_path, const char* protocol,
                                  int Lprotocol, bool preserve_contents)
{
// splicer begin class.Group.method.load_1_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_file_path(file_path, Lfile_path);
  const std::string SH_protocol(protocol, Lprotocol);
  SH_this->load(SH_file_path, SH_protocol, preserve_contents);
  return;
// splicer end class.Group.method.load_1_bufferify
}

void SIDRE_group_load_external_data(SIDRE_group* self, const char* file_path)
{
// splicer begin class.Group.method.load_external_data
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_file_path(file_path);
  SH_this->loadExternalData(SH_file_path);
  return;
// splicer end class.Group.method.load_external_data
}

void SIDRE_group_load_external_data_bufferify(SIDRE_group* self,
                                              const char* file_path,
                                              int Lfile_path)
{
// splicer begin class.Group.method.load_external_data_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_file_path(file_path, Lfile_path);
  SH_this->loadExternalData(SH_file_path);
  return;
// splicer end class.Group.method.load_external_data_bufferify
}

bool SIDRE_group_rename(SIDRE_group* self, const char* new_name)
{
// splicer begin class.Group.method.rename
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_new_name(new_name);
  bool SHC_rv = SH_this->rename(SH_new_name);
  return SHC_rv;
// splicer end class.Group.method.rename
}

bool SIDRE_group_rename_bufferify(SIDRE_group* self, const char* new_name,
                                  int Lnew_name)
{
// splicer begin class.Group.method.rename_bufferify
  axom::sidre::Group* SH_this =
    static_cast<axom::sidre::Group*>(static_cast<void*>(self));
  const std::string SH_new_name(new_name, Lnew_name);
  bool SHC_rv = SH_this->rename(SH_new_name);
  return SHC_rv;
// splicer end class.Group.method.rename_bufferify
}

}  // extern "C"
