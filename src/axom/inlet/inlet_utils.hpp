#include "axom/sidre.hpp"

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
* \brief This function appends the prefix path to the ending path.
*
* \param [in] The prefix string path.
* \param [in] The ending string path.
*
* \return The appended string.
*****************************************************************************
*/
std::string concatenatePaths(const std::string& prefix, const std::string& name);

/*!
*****************************************************************************
* \brief This function extracts the path from the full path.
*
* \param [in] The prefix of the path.
* \param [in] The full string path.
*
* \return The extracted string.
*****************************************************************************
*/
std::string getName(const std::string& prefix, const std::string& pathName);

}
}

#endif
