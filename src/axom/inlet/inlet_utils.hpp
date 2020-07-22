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

}
}

#endif
