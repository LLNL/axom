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

class VectorMetaBuffer : public MetaBuffer
{
public:
  virtual void *getDataPointer() const
    {
#ifdef USE_CXX11
	return m_context->data();
#else
	return  m_context->empty() ? ATK_NULLPTR : &(*m_context)[0];
#endif
    }

  virtual size_t getNumberOfElements() const
    {
        return m_context->size();
    }

  virtual TypeID getTypeID() const
  {
      // XXX fix
      return CONDUIT_NATIVE_INT_DATATYPE_ID;
  }

  virtual void * allocate(TypeID type, SidreLength nitems) const
    {
      m_context->reserve(nitems);
      return getDataPointer();
    }

    virtual void release() const
    {
      m_context->clear();
    }

    virtual void * reallocate(TypeID type, SidreLength nitems) const
    {
      m_context->resize(nitems);
      return getDataPointer();
    }


  VectorMetaBuffer(std::vector<int> * context) :
    m_context(context)
  { }

private:
  std::vector<int> * m_context;
};

//----------------------------------------------------------------------

DataView *registerVectorNode(DataGroup * group,
			     const std::string& name,
			     std::vector<int> * vect)
{
  VectorMetaBuffer * metabuffer = new VectorMetaBuffer(vect);
  return group->createMetaBufferView(name, metabuffer);
}

} /* end namespace sidre */
} /* end namespace asctoolkit */
