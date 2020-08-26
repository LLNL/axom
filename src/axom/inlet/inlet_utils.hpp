// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/sidre.hpp"
#include "fmt/fmt.hpp"
#include "axom/core/utilities/StringUtilities.hpp"

#ifndef INLET_UTILS_HPP
#define INLET_UTILS_HPP

namespace axom
{
namespace inlet
{


/*!
*****************************************************************************
* \brief This function is used to mark if anything went wrong during the 
* defining phase of inlet so verify() will properly fail.
*
* \param [in] root Pointer to the Sidre Root Group where the warning flag 
* will be set.
*****************************************************************************
*/
void setWarningFlag(axom::sidre::Group* root);

/*!
*****************************************************************************
* \brief This function appends the prefix name to the ending name.
*
* \param [in] The prefix string name.
* \param [in] The ending string name.
*
* \return The appended string.
*****************************************************************************
*/
std::string appendPrefix(const std::string& prefix, const std::string& name);

/*!
*****************************************************************************
* \brief This function extracts the Table name from the full name.
*
* \param [in] The prefix of the name, to be removed.
* \param [in] The full name.
*
* \return The extracted string.
*****************************************************************************
*/
std::string removePrefix(const std::string& prefix, const std::string& name);

}
}

#endif
