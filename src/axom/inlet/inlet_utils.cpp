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

std::string appendPrefix(const std::string& prefix, const std::string& name)
{
  return (prefix == "") ? name : prefix + "/" + name;
}

std::string removePrefix(const std::string& prefix, const std::string& name) {
  SLIC_ASSERT_MSG(name.find(prefix) != std::string::npos, 
                  fmt::format("Provided name {0} does not contain prefix{1}", 
                  name, prefix));
  return name.substr(prefix.size()-1);
}

}
}
