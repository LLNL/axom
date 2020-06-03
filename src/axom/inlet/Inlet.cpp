// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Inlet.cpp
 *
 * \brief This file contains the class implementation of Inlet, the main class
 *        for the Inlet component.
 *******************************************************************************
 */

#include "axom/inlet/Inlet.hpp"

#include "fmt/fmt.hpp"
#include "axom/slic.hpp"

namespace axom
{
namespace inlet
{

bool Inlet::addGroup(const std::string& name,
                     const std::string& description)
{
  SLIC_ASSERT_MSG(m_reader != nullptr, "Inlet's Reader class not set");
  SLIC_ASSERT_MSG(m_ds != nullptr, "Inlet's Sidre Datastore class not set");

  axom::sidre::Group* root = m_ds->getRoot();
  if (root->hasGroup(name))
  {
    SLIC_WARNING("Inlet: Cannot create Group that already exists: " + name);
    return false;
  }
  root->createGroup(name);
  return true;
}

bool Inlet::addInt(const std::string& name,
                   const std::string& description,
                   int defaultValue)
{
  SLIC_ASSERT_MSG(m_reader != nullptr, "Inlet's Reader class not set");
  SLIC_ASSERT_MSG(m_ds != nullptr, "Inlet's Sidre Datastore class not set");

  int value;
  if(!m_reader->getInt(name, value))
  {
    value = defaultValue;
  }

  axom::sidre::Group* root = m_ds->getRoot();
  if (root->hasView(name))
  {
    SLIC_WARNING("Inlet: Cannot create Integer that already exists: " + name);
    return false;
  }
  root->createViewScalar(name, value);

  return true;
}

bool Inlet::addInt(const std::string& name,
                   const std::string& description,
                   bool required)
{
  SLIC_ASSERT_MSG(m_reader != nullptr, "Inlet's Reader class not set");
  SLIC_ASSERT_MSG(m_ds != nullptr, "Inlet's Sidre Datastore class not set");

  int value;
  if(m_reader->getInt(name, value))
  {
    axom::sidre::Group* root = m_ds->getRoot();
    if (root->hasView(name))
    {
      SLIC_WARNING("Inlet: Cannot create Integer that already exists: " + name);
      return false;
    }
    root->createViewScalar(name, value);
  }
  else if(required)
  {
    SLIC_ERROR("Inlet: Required integer was not found in input deck: " + name);
  }
  else
  {
    // not found in input deck but not required
    return false;
  }
  return true;
}

bool Inlet::get(const std::string& name, int& value)
{
  SLIC_ASSERT_MSG(m_ds != nullptr, "Inlet's Sidre Datastore class not set");

  axom::sidre::Group* root = m_ds->getRoot();
  
  if (!root->hasView(name))
  {
    return false;
  }

  axom::sidre::View* view = root->getView(name);
  if (view->getTypeID() != axom::sidre::INT_ID)
  {
    std::string error_msg = fmt::format("Integer Field named '{0}' was asked for"
                                        " but recieved type {1}",
                                        name, view->getTypeID());
    SLIC_ERROR(error_msg);
  }

  value = view->getScalar();
  return true;
}

} // end namespace inlet
} // end namespace axom
