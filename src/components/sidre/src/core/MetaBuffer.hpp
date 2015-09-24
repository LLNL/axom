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

// Standard C++ headers
#include <memory>
#include <map>
#include <string>
#include <vector>

#if defined(USE_CXX11)
#include <unordered_map>
#else
#include "boost/unordered_map.hpp"
#endif

// Other CS Toolkit headers
#include "slic/slic.hpp"

// SiDRe project headers
#include "Collections.hpp"
#include "DataView.hpp"
#include "SidreTypes.hpp"



namespace asctoolkit
{
namespace sidre
{

// using directives to make Conduit usage easier and less visible
using conduit::Node;
using conduit::Schema;

class DataBuffer;
class DataGroup;
class DataStore;

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
  virtual void *getDataPointer(void *context) const = 0;

  /*!
   * \brief Return total number of elements allocated by this MetaBuffer object.
   */
  virtual size_t getNumberOfElements(void *context) const = 0;

};


} /* end namespace sidre */
} /* end namespace asctoolkit */

#endif /* METABUFFER_HPP_ */
