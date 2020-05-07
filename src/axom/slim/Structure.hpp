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
    axom::slim::Group CreateGroup(const std::string& name, const std::string& description);
    axom::slim::Group CreateGroup(std::string&& rname, std::string&& rdescription);

    const std::vector<axom::slim::Group>& Groups() { return m_groups; };
private:
    std::vector<axom::slim::Group> m_groups;
};

} // end namespace slim
} // end namespace axom

#endif
