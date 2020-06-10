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

  void sidreGroup(axom::sidre::Group* group) { m_sidreGroup = group; };
  axom::sidre::Group* sidreGroup() { return m_sidreGroup; };

  // Functions that define the input deck schema
  axom::sidre::Group* addGroup(const std::string& name,
                               const std::string& description);

  axom::sidre::Group* addBool(const std::string& name,
                              const std::string& description);
  axom::sidre::Group* addDouble(const std::string& name,
                                const std::string& description);
  axom::sidre::Group* addInt(const std::string& name,
                             const std::string& description);
  axom::sidre::Group* addString(const std::string& name,
                                const std::string& description);

  // Functions that get the values out of the datastore
  bool get(const std::string& name, bool& value);
  bool get(const std::string& name, double& value);
  bool get(const std::string& name, int& value);
  bool get(const std::string& name, std::string& value);

  // TODO add update value functions
private:
  axom::sidre::Group* add(const std::string& name,
                          const std::string& description);
  axom::sidre::View* get(const std::string& name);

  Reader* m_reader = nullptr;
  axom::sidre::Group* m_sidreGroup = nullptr;
};

} // end namespace inlet
} // end namespace axom

#endif
