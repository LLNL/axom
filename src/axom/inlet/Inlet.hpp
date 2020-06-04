// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Inlet.hpp
 *
 * \brief This file contains the class definition of Inlet, the main class
 *        for the Inlet component.
 *******************************************************************************
 */

#ifndef INLET_INLET_HPP
#define INLET_INLET_HPP

#include <string>
#include <vector>

#include "axom/inlet/Reader.hpp"

#include "axom/sidre.hpp"

namespace axom
{
namespace inlet
{

class Inlet
{
public:
  void reader(Reader* reader) { m_reader = reader; };
  Reader* reader() { return m_reader; };

  //TODO change to datastore group
  void datastore(axom::sidre::DataStore* ds) { m_ds = ds; };
  axom::sidre::DataStore* datastore() { return m_ds; };

  // Functions that define the input deck schema

  //TODO rename group to not conflict with sidre group
  bool addGroup(const std::string& name, const std::string& description);

  bool addInt(const std::string& name,
              const std::string& description,
              int defaultValue);
  bool addInt(const std::string& name,
              const std::string& description,
              bool required=false);

  // Functions that get the values out of the datastore
  
  bool get(const std::string& name, int& value);

  // TODO add update value functions
private:
  Reader* m_reader = nullptr;
  axom::sidre::DataStore* m_ds = nullptr;
};

} // end namespace inlet
} // end namespace axom

#endif
