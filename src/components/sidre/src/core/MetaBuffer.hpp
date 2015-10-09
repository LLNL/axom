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
 * \file    MetaBuffer.hpp
 *
 * \brief   Header file containing definition of MetaBuffer class.
 *
 ******************************************************************************
 */

#ifndef METABUFFER_HPP_
#define METABUFFER_HPP_

// SiDRe project headers
#include <sidre/SidreTypes.h>

namespace asctoolkit
{
namespace sidre
{

/*!
 * \class MetaBuffer
 *
 * \brief Abstract class used to control behavior to access meta-buffer.
 *
 */
class MetaBuffer
{
public:

  /*!
   * \brief Return void-pointer to data associated with MetaBuffer.
   */
  virtual void *getDataPointer() const = 0;

  /*!
   * \brief Return total number of elements allocated by this MetaBuffer object.
   */
  virtual size_t getNumberOfElements() const = 0;

  /*!
   * \brief Return total number of bytes allocated by this DataView object.
   */
  virtual TypeID getTypeID() const = 0;

  /*!
   * \brief Allocate buffer to requested size.
   */
  virtual void * allocate(TypeID type, SidreLength nitems) const = 0;

};


} /* end namespace sidre */
} /* end namespace asctoolkit */

#endif /* METABUFFER_HPP_ */
