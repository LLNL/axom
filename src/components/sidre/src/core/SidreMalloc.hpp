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
 * \file    SidreMalloc.hpp
 *
 * \brief   Header file containing definition to a malloc based MetaBuffer.
 *
 ******************************************************************************
 */

#ifndef SIDREMALLOC_HPP_
#define SIDREMALLOC_HPP_

// Standard C++ headers
#include <string>

// Other CS Toolkit headers
#include "common/CommonTypes.hpp"
#include "conduit/conduit.hpp"

// SiDRe project headers
#include "sidre/DataGroup.hpp"
#include "sidre/DataView.hpp"
//#include "sidre/MetaBuffer.hpp"

namespace asctoolkit
{
namespace sidre
{
class DataView;

DataView *registerStaticNode(DataGroup * group,
			     const std::string& name,
			     void *addr,
			     TypeID type, SidreLength len );

// XXX - create overloads for each native type
DataView *registerStaticNode(DataGroup * group,
			     const std::string& name,
			     int *addr,
			     SidreLength len = 1 );

DataView *registerMallocNode(DataGroup * group,
			     const std::string& name);

DataView *registerMallocNode(DataGroup * group,
			     const std::string& name,
			     TypeID type, SidreLength len );


} /* end namespace sidre */
} /* end namespace asctoolkit */

#endif /* SIDREMALLOC_HPP_ */

