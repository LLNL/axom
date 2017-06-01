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
 * \file AttrValues.cpp
 *
 * \brief   Implementation file for AttrValues class.
 *
 ******************************************************************************
 */

// Associated header file
#include "AttrValues.hpp"

// Other axom headers
#include "axom/Types.hpp"
#include "slic/slic.hpp"

// Sidre component headers

namespace axom
{
namespace sidre
{

/*
 *************************************************************************
 *
 *
 *************************************************************************
 */
bool AttrValues::setAttrValue(const Attribute * attr,
			      const std::string & value)
{

  if (m_values == AXOM_NULLPTR)
  {
    m_values = new(std::nothrow) Values( );
  }

  IndexType iattr = attr->getIndex();

  if ((size_t) iattr < m_values->size())
  {
    // replace existing attribute
    (*m_values)[iattr] = value;
  }
  else
  {
    // Create all attributes up to iattr
    m_values->reserve(iattr + 1);
    for(int n=m_values->size(); n < iattr; ++n)
    {
      m_values->push_back("");
    }
    m_values->push_back(value);
  } 


#if 0
std::vector<int>::size_type sz = myvector.size();
 std::vector<std::string> * avec = static_cast(std::vector<std::string> *) m_attributes[iattr];
  // Make sure there are enough elements
#endif

  return true;
}


/*
 *************************************************************************
 *
 *
 *************************************************************************
 */
const std::string & AttrValues::getAttribute( const Attribute * attr ) const
{
  static std::string none("NONE");

  if (m_values == AXOM_NULLPTR)
  {
    // No attributes have been set;
    return attr->getDefault();
  }

  IndexType iattr = attr->getIndex();

  if ((size_t) iattr > m_values->size())
  {
    // This attribute has not been set for this View
    return attr->getDefault();
  }

  return (*m_values)[iattr];
}


/*
 *************************************************************************
 *
 * PRIVATE ctor for AttrValues
 *
 *************************************************************************
 */
AttrValues::AttrValues() :
  m_values(AXOM_NULLPTR)
{}

/*
 *************************************************************************
 *
 * PRIVATE dtor.
 *
 *************************************************************************
 */
AttrValues::~AttrValues()
{
}

} /* end namespace sidre */
} /* end namespace axom */
