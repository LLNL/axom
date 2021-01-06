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
  return checkIfRequired(*m_sidreGroup, *m_sidreRootGroup);
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
  verified &=
    verifyRequired(*m_sidreGroup, static_cast<bool>(m_func), "Function");
  // Verify this Function if a lambda was configured
  if(m_verifier && !m_verifier(*this))
  {
    verified = false;
    SLIC_WARNING(fmt::format("[Inlet] Function failed verification: {0}", name()));
  }

  return verified;
}

}  // end namespace inlet
}  // end namespace axom
