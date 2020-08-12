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

std::string getFullName(const std::string& prefix, const std::string& name)
{
  if (prefix == "")
  {
    return name;
  }
  else
  {
    return prefix + "/" + name;
  }
}

std::string getPath(const std::string& prefix, const std::string& pathName) {
  return pathName.substr(0, prefix.size());
}

}
}
