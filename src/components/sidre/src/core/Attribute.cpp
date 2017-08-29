/*
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
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
 * \file Attribute.cpp
 *
 * \brief   Implementation file for Attribute class.
 *
 ******************************************************************************
 */

// Associated header file
#include "Attribute.hpp"

// Other axom headers
//#include "axom/Types.hpp"
//#include "slic/slic.hpp"

// Sidre component headers
//#include "Buffer.hpp"
//#include "Group.hpp"
//#include "DataStore.hpp"

namespace axom
{
namespace sidre
{

/*
 *************************************************************************
 *
 * PRIVATE ctor for Attribute not associated with any data.
 *
 *************************************************************************
 */
Attribute::Attribute( const std::string& name) :
  m_name(name),
  m_index(InvalidIndex)
{}

/*
 *************************************************************************
 *
 * PRIVATE dtor.
 *
 *************************************************************************
 */
Attribute::~Attribute()
{}

} /* end namespace sidre */
} /* end namespace axom */
