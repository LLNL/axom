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
#include "sidre/SidreVector.hpp"

namespace asctoolkit
{
namespace sidre
{

class VectorContext {
public:
    VectorContext () :
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
class VectorMetaBuffer : public MetaBuffer
{
public:
  virtual void *getDataPointer(void *context) const
    {
	std::vector<int> * mcontext = static_cast<std::vector<int> *>(context);
#ifdef USE_CXX11
	return mcontext->data();
#else
	return  mcontext->empty() ? ATK_NULLPTR : &(*mcontext)[0];
#endif
    }

  virtual size_t getNumberOfElements(void *context) const
    {
	std::vector<int> * mcontext = static_cast<std::vector<int> *>(context);
        return mcontext->size();
    }

  virtual void allocate(void *context, TypeID type, SidreLength nitems) const
    {
	//	std::vector<int> * mcontext = static_cast<std::vector<int> *>(context);
	return;
    }
};

//----------------------------------------------------------------------

// XXX - temp code
VectorMetaBuffer *JUNK_DUMMY = NULL;
void RegisterVectorMetaBuffers(void);

//----------------------------------------------------------------------

DataView *registerVectorNode(DataGroup * group,
			     const std::string& name,
			     std::vector<int> *vect)
{
  if (JUNK_DUMMY == NULL) {
      RegisterVectorMetaBuffers();
  }
  DataView * view  = group->createViewWithMetaBuffer(name,
						     static_cast<void *>(vect),
						     JUNK_DUMMY);
  return view;
}

void RegisterVectorMetaBuffers(void)
{
    JUNK_DUMMY = new VectorMetaBuffer;
}

} /* end namespace sidre */
} /* end namespace asctoolkit */
