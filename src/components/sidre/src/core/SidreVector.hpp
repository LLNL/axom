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

#ifndef SIDREVECTOR_HPP_
#define SIDREVECTOR_HPP_

// Standard C++ headers
#include <string>
#include <vector>

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

DataView *registerVectorNode(DataGroup * group,
			     const std::string& name,
			     std::vector<int> *vect);

} /* end namespace sidre */
} /* end namespace asctoolkit */

#endif /* SIDREVECTOR_HPP_ */

