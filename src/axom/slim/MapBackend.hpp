// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file MapBackend.hpp
 *
 * \brief This file contains the class definition of MapBackend.
 *******************************************************************************
 */

#ifndef SLIM_MAPBACKEND_HPP
#define SLIM_MAPBACKEND_HPP

#include <string>
#include <map>

#include "axom/slim/Field.hpp"
#include "axom/slim/Backend.hpp"

namespace axom
{
namespace slim
{

class MapBackend : public Backend
{
public:
    void add(Field* field) { m_fields[field->name()] = field; }
    Field* get(const std::string& name)
    {
        auto p = m_fields.find(name);
        if(p == m_fields.end()) {
            return nullptr;
        }
        return p->second;
    }
private:
    std::map<std::string, Field*> m_fields;
};

} // end namespace slim
} // end namespace axom

#endif
