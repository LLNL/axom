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

#ifndef SLIM_STRUCTURE_HPP
#define SLIM_STRUCTURE_HPP

#include <string>
#include <vector>

#include "axom/slim/Field.hpp"
#include "axom/slim/IntField.hpp"
#include "axom/slim/GroupField.hpp"
#include "axom/slim/Map.hpp"

namespace axom
{
namespace slim
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

    GroupField* addGroup(const std::string& name, const std::string& description);
    GroupField* addGroup(std::string&& rname, std::string&& rdescription);

    IntField* addIntField(const std::string& name,
                         const std::string& description,
                         int defaultValue);
    IntField* addIntField(const std::string& name,
                         const std::string& description);

    const std::vector<Field*>& fields() { return m_fields; };
private:
    Map* m_map = nullptr;
    std::vector<Field*> m_fields;
};

} // end namespace slim
} // end namespace axom

#endif
