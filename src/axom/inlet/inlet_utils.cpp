#include "axom/inlet/inlet_utils.hpp"

namespace axom
{
namespace inlet
{


void setWarningFlag(axom::sidre::Group* root) {
  if (!root->hasView("warningFlag")) {
    root->createViewScalar("warningFlag", 1);
  }
}

}
}
