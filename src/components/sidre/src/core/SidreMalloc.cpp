/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

/*!
 ******************************************************************************
 *
 * \file    SidreMalloc.cpp
 *
 * \brief   Header file containing definition to a malloc based MetaBuffer.
 *
 ******************************************************************************
 */

// Standard C++ headers
//#include <string>
#include <cstdlib>

// Other CS Toolkit headers
//#include "common/CommonTypes.hpp"

// SiDRe project headers
#include "sidre/SidreMalloc.hpp"

namespace asctoolkit
{
namespace sidre
{

class MallocContext {
public:
    MallocContext () :
	m_type(CONDUIT_EMPTY_T),
	m_nitems(0),
	m_number_of_bytes(0),
	m_data_pointer(ATK_NULLPTR)
    {};
    // m_bytes_per_item?
    TypeID m_type;
    SidreLength m_nitems;
    size_t m_number_of_bytes;
    void *m_data_pointer;
};

//---------- Static
class StaticMetaBuffer : public MetaBuffer
{
public:
  virtual void *getDataPointer(void *context) const
    {
	MallocContext * mcontext = static_cast<MallocContext *>(context);
	return mcontext->m_data_pointer;
    }

  virtual size_t getNumberOfElements(void *context) const
    {
	MallocContext * mcontext = static_cast<MallocContext *>(context);
        return mcontext->m_nitems;
    }

    virtual TypeID getTypeID(void *context) const
    {
	MallocContext * mcontext = static_cast<MallocContext *>(context);
	return mcontext->m_type;
    }

  virtual void allocate(void *context, TypeID type, SidreLength nitems) const
    {
	SLIC_ASSERT_MSG( false,
			 "Cannot call allocate on a static buffer");
	return;
    }

private:
};

//----------------------------------------------------------------------

// XXX - temp code
StaticMetaBuffer *STATIC_JUNK_DUMMY = NULL;
void RegisterStaticMetaBuffers(void);

//----------------------------------------------------------------------

DataView *registerStaticNode(DataGroup * group,
			     const std::string& name,
			     void *addr,
			     TypeID type, SidreLength len )
{
  if (STATIC_JUNK_DUMMY == NULL) {
      RegisterStaticMetaBuffers();
  }
  MallocContext * context = new MallocContext;

  context->m_type = type;
  context->m_nitems = len;
  context->m_number_of_bytes = 0;  // XXX
  context->m_data_pointer = addr;

  DataView * view  = group->createViewWithMetaBuffer(name,
						     context,
						     STATIC_JUNK_DUMMY);
  return view;
}

// XXX - create overloads for each native type
DataView *registerStaticNode(DataGroup * group,
			     const std::string& name,
			     int *addr,
			     SidreLength len)
{
    return  registerStaticNode(group, name,
			       static_cast<void *>(addr),
			       CONDUIT_NATIVE_INT_DATATYPE_ID, len);
}

void RegisterStaticMetaBuffers(void)
{
    STATIC_JUNK_DUMMY = new StaticMetaBuffer;
}

//----------------------------------------------------------------------

class MallocMetaBuffer : public MetaBuffer
{
public:
  virtual void *getDataPointer(void *context) const
    {
	MallocContext * mcontext = static_cast<MallocContext *>(context);
	return mcontext->m_data_pointer;
    }

  virtual size_t getNumberOfElements(void *context) const
    {
	MallocContext * mcontext = static_cast<MallocContext *>(context);
        return mcontext->m_nitems;
    }

    virtual TypeID getTypeID(void *context) const
    {
	MallocContext * mcontext = static_cast<MallocContext *>(context);
	return mcontext->m_type;
    }

  virtual void allocate(void *context, TypeID type, SidreLength nitems) const
    {
	MallocContext * mcontext = static_cast<MallocContext *>(context);
	// XXX - replace sizeof
	mcontext->m_nitems = nitems;
	mcontext->m_number_of_bytes = sizeof(int) * nitems;
	mcontext->m_data_pointer = malloc( mcontext->m_number_of_bytes); 
	mcontext->m_type = type;
	return;
    }

private:
};

//----------------------------------------------------------------------

// XXX - temp code
MallocMetaBuffer *MALLOC_JUNK_DUMMY = NULL;
void RegisterMallocMetaBuffers(void);

//----------------------------------------------------------------------

    // XXX TODO - release context

DataView *registerMallocNode(DataGroup * group,
			     const std::string& name)
{
  if (MALLOC_JUNK_DUMMY == NULL) {
      RegisterMallocMetaBuffers();
  }
  MallocContext * context = new MallocContext;
  DataView *  view  = group->createViewWithMetaBuffer(name,
						      context,
						      MALLOC_JUNK_DUMMY);
  return view;
}

DataView *registerMallocNode(DataGroup * group,
			     const std::string& name,
			     TypeID type, SidreLength len )
{
  if (MALLOC_JUNK_DUMMY == NULL) {
      RegisterMallocMetaBuffers();
  }
  MallocContext * context = new MallocContext;
  DataView *  view  = group->createViewWithMetaBuffer(name,
						      context,
						      MALLOC_JUNK_DUMMY);
  view->m_type = type;
  view->m_nitems = len;
  return view;
}


void RegisterMallocMetaBuffers(void)
{
    MALLOC_JUNK_DUMMY = new MallocMetaBuffer;
}

} /* end namespace sidre */
} /* end namespace asctoolkit */
