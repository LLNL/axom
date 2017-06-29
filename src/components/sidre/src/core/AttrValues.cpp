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
bool AttrValues::hasValue( const Attribute * attr ) const
{
  if (attr == AXOM_NULLPTR)
  {
    SLIC_CHECK_MSG(attr != AXOM_NULLPTR,
		   "hasValue: called without an Attribute");
    return false;
  }

  if (m_values == AXOM_NULLPTR)
  {
    // No attributes have been set in this View.
    return false;
  }

  IndexType iattr = attr->getIndex();

  if ((size_t) iattr >= m_values->size())
  {
    // This attribute has not been set for this View.
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
 * Set Attribute to its default value. If the Attribute has not been
 * set yet, then it is already the default value.  Otherwise, reset
 * the Node to empty.  Future calls to hasAttributeValue will return
 * false for the attribute.
 *
 *************************************************************************
 */
bool AttrValues::setDefault( const Attribute * attr )
{
  if (attr == AXOM_NULLPTR)
  {
    SLIC_CHECK_MSG(attr != AXOM_NULLPTR,
		   "setDefault: called without an Attribute");
    return false;
  }

  if (m_values == AXOM_NULLPTR)
  {
    // No attributes have been set in this View, already default.
    return true;
  }

  IndexType iattr = attr->getIndex();

  if ((size_t) iattr >= m_values->size())
  {
    // This attribute has not been set for this View, already default.
    return true;
  }

  Node & value = (*m_values)[iattr];
  value.reset();

  return true;
}

/*
 *************************************************************************
 *
 * PRIVATE Create a Node for the Attribute.
 *
 * Create vector of Nodes and push empty nodes up to attr's index.
 *
 *************************************************************************
 */
bool AttrValues::createNode(const Attribute * attr)
{
  if (attr == AXOM_NULLPTR)
  {
    SLIC_CHECK_MSG(attr != AXOM_NULLPTR,
		   "createNode: called without an Attribute");
    return false;
  }

  if (m_values == AXOM_NULLPTR)
  {
    m_values = new(std::nothrow) Values( );
  }

  IndexType iattr = attr->getIndex();

  if ((size_t) iattr >= m_values->size())
  {
    // Create all attributes up to iattr, push back empty Nodes
    m_values->reserve(iattr + 1);
    for(int n=m_values->size(); n < iattr + 1; ++n)
    {
      m_values->push_back(Node());
    }
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
Node::ConstValue AttrValues::getScalar( const Attribute * attr ) const
{
  if (attr == AXOM_NULLPTR)
  {
    SLIC_CHECK_MSG(attr != AXOM_NULLPTR,
		   "getScalar: called without an Attribute");
    return getEmptyNodeRef().value();
  }

  const Node & node = getValueNodeRef(attr);
  return node.value();
}

/*
 *************************************************************************
 *
 * Return attribute.
 *
 *************************************************************************
 */
const char * AttrValues::getString( const Attribute * attr ) const
{
  if (attr == AXOM_NULLPTR)
  {
    SLIC_CHECK_MSG(attr != AXOM_NULLPTR,
		   "getString: called without an Attribute");
    return AXOM_NULLPTR;
  }

  if (attr->getTypeID() != CHAR8_STR_ID)
  {
    SLIC_CHECK_MSG(attr->getTypeID() == CHAR8_STR_ID,
		   "getString: Called on attribute '"
		   << attr->getName() 
		   << "' which is type "
		   << DataType::id_to_name(attr->getTypeID())
                   << ".");
    return AXOM_NULLPTR;
  }

  const Node & node = getValueNodeRef(attr);
  return node.as_char8_str();
}

/*
 *************************************************************************
 *
 * Return Reference to Node.
 *
 *************************************************************************
 */
const Node & AttrValues::getValueNodeRef( const Attribute * attr ) const
{
  if (attr == AXOM_NULLPTR)
  {
    SLIC_CHECK_MSG(attr != AXOM_NULLPTR,
		   "getValueNodeRef: called without an Attribute");
    return getEmptyNodeRef();
  }

  if (m_values == AXOM_NULLPTR)
  {
    // No attributes have been set in this View;
    return attr->getDefaultNodeRef();
  }

  IndexType iattr = attr->getIndex();

  if ((size_t) iattr >= m_values->size())
  {
    // This attribute has not been set for this View
    return attr->getDefaultNodeRef();
  }

  Node & value = (*m_values)[iattr];
  
  if (value.schema().dtype().is_empty())
  {
    return attr->getDefaultNodeRef();
  }

  return value;
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
