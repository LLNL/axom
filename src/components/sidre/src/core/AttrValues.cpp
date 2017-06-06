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
 * Return true if the attribute has been explicitly set; else false.
 *
 *************************************************************************
 */
bool AttrValues::hasAttrValue( const Attribute * attr ) const
{
  if (m_values == AXOM_NULLPTR)
  {
    // No attributes have been set in this View;
    return false;
  }

  IndexType iattr = attr->getIndex();

  if ((size_t) iattr > m_values->size())
  {
    // This attribute has not been set for this View
    return false;
  }

  Node & value = (*m_values)[iattr];
  
  if (value.schema().dtype().is_empty())
  {
    return false;
  }

  return true;
}

/*
 *************************************************************************
 *
 * Set Attribute.
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
    // Create all attributes up to iattr, push back empty Nodes
    m_values->reserve(iattr + 1);
    for(int n=m_values->size(); n < iattr + 1; ++n)
    {
      m_values->push_back(Node());
    }
    (*m_values)[iattr] = value;
  } 

  return true;
}


/*
 *************************************************************************
 *
 * Return attribute.
 *
 *************************************************************************
 */
const char * AttrValues::getAttribute( const Attribute * attr ) const
{
  if (m_values == AXOM_NULLPTR)
  {
    // No attributes have been set in this View;
    return attr->getDefault();
  }

  IndexType iattr = attr->getIndex();

  if ((size_t) iattr > m_values->size())
  {
    // This attribute has not been set for this View
    return attr->getDefault();
  }

  Node & value = (*m_values)[iattr];
  
  if (value.schema().dtype().is_empty())
  {
    return attr->getDefault();
  }

  return value.as_char8_str();
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
