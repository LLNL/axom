// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Structure.cpp
 *
 * \brief This file contains the class implementation of Structure.
 *******************************************************************************
 */

#include "axom/slim/Structure.hpp"

#include "axom/slim/Group.hpp"

namespace axom
{
namespace slim
{

axom::slim::Group Structure::CreateGroup(const std::string& name,
                                         const std::string& description)
{
    axom::slim::Group group = axom::slim::Group(name, description);
    m_groups(group);
    return group;
}

axom::slim::Group Structure::CreateGroup(std::string&& rname,
                                         std::string&& rdescription)
{
    axom::slim::Group group = axom::slim::Group(rname, rdescription);
    m_groups(group);
    return group;
}

} // end namespace slim
} // end namespace axom

#endif
