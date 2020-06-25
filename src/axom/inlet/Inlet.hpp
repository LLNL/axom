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

#include <memory>
#include <string>
#include <vector>

//#include "axom/inlet/SchemaCreator.hpp"
#include "axom/inlet/Group.hpp"
#include "axom/inlet/Field.hpp"
#include "axom/inlet/Reader.hpp"

#include "axom/sidre.hpp"

namespace axom
{
namespace inlet
{

class Inlet : public SchemaCreator
{
public:
  Inlet(std::shared_ptr<Reader> reader,
        axom::sidre::Group* sidreRootGroup) :
    m_reader(reader),
    m_sidreRootGroup(sidreRootGroup),
    m_group(std::make_shared<Group>("", "", m_reader, m_sidreRootGroup)) {}

  virtual ~Inlet() = default;

  std::shared_ptr<Reader> reader() { return m_reader; };
  axom::sidre::Group* sidreGroup() { return m_sidreRootGroup; };

  // Functions that define the input deck schema
  std::shared_ptr<Group> addGroup(const std::string& name,
                                  const std::string& description);

  std::shared_ptr<Field> addBool(const std::string& name,
                                 const std::string& description);
  std::shared_ptr<Field> addDouble(const std::string& name,
                                   const std::string& description);
  std::shared_ptr<Field> addInt(const std::string& name,
                                const std::string& description);
  std::shared_ptr<Field> addString(const std::string& name,
                                   const std::string& description);

  // Functions that get the values out of the datastore
  bool get(const std::string& name, bool& value);
  bool get(const std::string& name, double& value);
  bool get(const std::string& name, int& value);
  bool get(const std::string& name, std::string& value);

  // TODO add update value functions
private:
  axom::sidre::View* baseGet(const std::string& name);

  std::shared_ptr<Reader> m_reader;
  axom::sidre::Group* m_sidreRootGroup = nullptr;
  std::shared_ptr<Group> m_group;
};

} // end namespace inlet
} // end namespace axom

#endif
