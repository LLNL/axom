// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "Function.hpp"

namespace axom
{
namespace inlet
{
std::string Function::name() const
{
  return removePrefix(m_sidreRootGroup->getPathName(),
                      m_sidreGroup->getPathName());
}

Function& Function::required(bool isRequired)
{
  SLIC_ASSERT_MSG(m_sidreGroup != nullptr,
                  "[Inlet] Function specific Sidre Datastore Group not set");
  setRequired(*m_sidreGroup, *m_sidreRootGroup, isRequired);
  return *this;
}

bool Function::isRequired() const
{
  SLIC_ASSERT_MSG(m_sidreGroup != nullptr,
                  "[Inlet] Function specific Sidre Datastore Group not set");
  return checkRequired(*m_sidreGroup, *m_sidreRootGroup);
}

Function& Function::registerVerifier(std::function<bool(const Function&)> lambda)
{
  SLIC_WARNING_IF(m_verifier,
                  fmt::format("[Inlet] Verifier for Function "
                              "already set: {0}",
                              name()));
  m_verifier = lambda;
  return *this;
}

bool Function::verify() const
{
  bool verified = true;
  // If this function was required, make sure something was defined in it
  if(m_sidreGroup->hasView("required"))
  {
    int8 required = m_sidreGroup->getView("required")->getData();
    if(required && !m_func)
    {
      const std::string msg = fmt::format(
        "[Inlet] Required Function not "
        "specified: {0}",
        m_sidreGroup->getPathName());
      SLIC_WARNING(msg);
      verified = false;
    }
  }
  // Verify this Function if a lambda was configured
  if(m_verifier && !m_verifier(*this))
  {
    verified = false;
    SLIC_WARNING(fmt::format("[Inlet] Function failed verification: {0}", name()));
  }

  return verified;
}

bool AggregateFunction::verify() const
{
  return std::all_of(
    m_funcs.begin(),
    m_funcs.end(),
    [](const Verifiable<Function>& func) { return func.verify(); });
}

AggregateFunction& AggregateFunction::required(bool isRequired)
{
  for(auto& func : m_funcs)
  {
    func.get().required(isRequired);
  }
  return *this;
}

bool AggregateFunction::isRequired() const
{
  return std::any_of(
    m_funcs.begin(),
    m_funcs.end(),
    [](const Verifiable<Function>& func) { return func.isRequired(); });
}

AggregateFunction& AggregateFunction::registerVerifier(
  std::function<bool(const Function&)> lambda)
{
  for(auto& func : m_funcs)
  {
    func.get().registerVerifier(lambda);
  }
  return *this;
}

}  // end namespace inlet
}  // end namespace axom
