// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Group.hpp
 *
 * \brief This file contains the class definition of Group.
 *******************************************************************************
 */

#ifndef SLIM_GROUP_HPP
#define SLIM_GROUP_HPP

#include <string>

namespace axom
{
namespace slim
{


/*!
 *******************************************************************************
 * \class Group
 *
 * \brief This class is used to define a group for your input deck.
 *
 *******************************************************************************
 */
class Group
{
public:
    Group(const std::string& name, const std::string& description)
    : m_name(name)
    , m_description(description) {};

    Group(std::string&& rname, std::string&& rdescription)
    : m_name(std::move(rname))
    , m_description(std::move(rdescription)) {};

private:
    std::string m_name;
    std::string m_description;
};

} // end namespace slim
} // end namespace axom

#endif
