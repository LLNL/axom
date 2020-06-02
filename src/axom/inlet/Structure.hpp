// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Structure.hpp
 *
 * \brief This file contains the class definition of Structure.
 *******************************************************************************
 */

#ifndef INLET_STRUCTURE_HPP
#define INLET_STRUCTURE_HPP

#include <string>
#include <vector>

#include "axom/inlet/Field.hpp"
#include "axom/inlet/IntField.hpp"
#include "axom/inlet/GroupField.hpp"
#include "axom/inlet/Map.hpp"
#include "axom/inlet/Backend.hpp"

namespace axom
{
namespace inlet
{


/*!
 *******************************************************************************
 * \class Structure
 *
 * \brief This class is used to define the structure of your input deck.
 *
 *******************************************************************************
 */
class Structure
{
public:
  void map(Map* map) { m_map = map; };
  Map* map() { return m_map; };

  void backend(Backend* backend) { m_backend = backend; }
  std::vector<std::string> names();

  GroupField* addGroup(const std::string& name, const std::string& description);
  GroupField* addGroup(std::string&& rname, std::string&& rdescription);

  IntField* addIntField(const std::string& name,
                        const std::string& description,
                        int defaultValue);
  IntField* addIntField(const std::string& name,
                        const std::string& description,
                        bool required=false);

  IntField* getIntField(const std::string& name);
private:
  Map* m_map = nullptr;
  Backend* m_backend = nullptr;
};

} // end namespace inlet
} // end namespace axom

#endif
