// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef INLET_GROUP_HPP
#define INLET_GROUP_HPP

#include <memory>
#include <string>

#include "axom/inlet/Field.hpp"
#include "axom/inlet/Reader.hpp"
#include "axom/inlet/SchemaCreator.hpp"

#include "axom/sidre.hpp"

namespace axom
{
namespace inlet
{


//TODO: possibly rename Group not to clash with sidre::Group
class Group : public SchemaCreator
{
public:
  Group(const std::string& name,
        const std::string& description,
        std::shared_ptr<Reader> reader,
        axom::sidre::Group* sidreRootGroup) : 
    m_name(name),
    m_reader(reader),
    m_sidreRootGroup(sidreRootGroup)
    {
      SLIC_ASSERT_MSG(m_reader, "Inlet's Reader class not valid");
      SLIC_ASSERT_MSG(m_sidreRootGroup != nullptr, "Inlet's Sidre Datastore class not set");

      axom::sidre::Group* sidreGroup = nullptr;
      if (m_name == "")
      {
        sidreGroup = m_sidreRootGroup;
      }
      else
      {
        if (!m_sidreRootGroup->hasGroup(name))
        {
          sidreGroup = m_sidreRootGroup->createGroup(name);
        }
        else
        {
          sidreGroup = m_sidreRootGroup->getGroup(name);
        }
      }

      if(description == "")
      {
        if (sidreGroup->hasView("description"))
        {
          //TODO: warn user?
          sidreGroup->destroyViewAndData("description");
        }
        sidreGroup->createViewString("description", description);
      }
    }

  virtual ~Group() = default;

  axom::sidre::Group* sidreGroup() { return m_sidreRootGroup->getGroup(m_name); };

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
private:
  axom::sidre::Group* baseFieldAdd(const std::string& name,
                                   const std::string& description);

  std::string m_name;
  std::shared_ptr<Reader> m_reader;
  axom::sidre::Group* m_sidreRootGroup;
};

} // end namespace inlet
} // end namespace axom

#endif
