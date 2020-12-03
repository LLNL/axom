// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <memory>
#include <utility>

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
* \brief This function is used to configure the Inlet object corresponding
* to the provided Sidre group as required
*
* \param [in] target Reference to the Sidre group to set the required 
* status of
* \param [in] root Reference to the Sidre Root Group where the warning flag 
* will be set on failure
* \param [in] required Whether the Inlet object is required
*****************************************************************************
*/
void setRequired(axom::sidre::Group& target,
                 axom::sidre::Group& root,
                 bool required);

/*!
*****************************************************************************
* \brief This function is used to determine if the Inlet object corresponding
* to the provided Sidre group is required
*
* \param [in] target Reference to the Sidre group to check the required 
* status of
* \param [in] root Reference to the Sidre Root Group where the warning flag 
* will be set on failure
* \return Whether the Inlet object is required
*****************************************************************************
*/
bool checkIfRequired(const axom::sidre::Group& target, axom::sidre::Group& root);

/*!
*****************************************************************************
* \brief This function is used to verify the required-ness of the Inlet object
* corresponding to the provided Sidre group
*
* \param [in] target Reference to the Sidre group to verify the required-ness of
* \param [in] condition The condition that must be true if the object is required
* \param [in] type The type of the object as a string, for use in the warning message
* 
* \return False if the object was required but \p condition was false, True otherwise
* \post If the function returns False, a warning message will be emitted
*****************************************************************************
*/
bool verifyRequired(const axom::sidre::Group& target,
                    const bool condition,
                    const std::string& type);

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

namespace cpp11_compat
{
/*!
*****************************************************************************
* \brief This function provides backwards compatibility for std::make_unique,
* which is not implemented until C++14.  It should be removed when either
* Axom or the Inlet component is no longer required to support C++11
*
* \tparam T The type to construct
* \tparam Args The variadic argument list to forward to T's constructor
*
* \return A unique ptr constructed with the given arguments
*****************************************************************************
*/
template <typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

}  // namespace cpp11_compat

}  // namespace inlet
}  // namespace axom

#endif
