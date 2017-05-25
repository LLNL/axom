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
bool AttrValues::setAttrValue(const Attribute * attr, IndexType idx,
			      const std::string & value)
{

  if (m_attributes == AXOM_NULLPTR)
  {
    m_attributes = new(std::nothrow) Attributes( );
  }

  IndexType iattr = attr->getIndex();

  // If the vector for this attribute does not exist yet, create it.
  // Set new attributes to AXOM_NULLPTR.
  if ((size_t) iattr >= m_attributes->size()) {
    m_attributes->reserve(iattr + 1);
    for(int n=m_attributes->size(); n < iattr; ++n)
    {
      m_attributes->push_back(AXOM_NULLPTR);
    }
    std::vector<std::string> * avec = new std::vector<std::string>();
    m_attributes->push_back(avec);
      //    (*m_attributes)[iattr] = avec;
  }

  {
    //    std::vector<std::string> * avec = static_cast<std::vector<std::string> *>(m_attributes[iattr]);
    std::vector<std::string> * avec = (*m_attributes)[iattr];
    //    std::vector<std::string> * avec = m_attributes->at(iattr);

    if ((size_t) idx < avec->size())
    {
      // replace existing attribute
      (*avec)[idx] = value;
    }
    else
    {
      // Create attributes for all Views upto idx
      avec->reserve(idx + 1);
      for(int n=avec->size(); n < idx; ++n)
      {
	avec->push_back("");
      }
      avec->push_back(value);
    }
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
const std::string & AttrValues::getAttribute( const Attribute * attr, IndexType idx ) const
{
  static std::string none("NONE");

  if (m_attributes == AXOM_NULLPTR)
  {
    // No attributes have been set;
    return none;
  }

  IndexType iattr = attr->getIndex();

  if ((size_t) iattr > m_attributes->size()) {
    // This attribute has not been set for any Views in this Group.
    return none;
  }

  std::vector<std::string> * avec = (*m_attributes)[iattr];
  if ((size_t) idx > avec->size())
  {
    // This attribute has not been set for this View
    return none;
  }

  return (*avec)[iattr];
}


/*
 *************************************************************************
 *
 * PRIVATE ctor for AttrValues
 *
 *************************************************************************
 */
AttrValues::AttrValues() :
  m_attributes(AXOM_NULLPTR)
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
