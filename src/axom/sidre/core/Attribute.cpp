/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
//#include "axom/core/Types.hpp"
//#include "axom/slic/interface/slic.hpp"

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
