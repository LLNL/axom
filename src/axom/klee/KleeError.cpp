// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/KleeError.hpp"

namespace axom
{
namespace klee
{
KleeError::KleeError(const inlet::VerificationError &error) : m_errors {{error}}
{ }

KleeError::KleeError(const std::vector<inlet::VerificationError> &errors)
  : m_errors {errors}
{ }

const char *KleeError::what() const noexcept
{
  return m_errors[0].message.data();
}

}  // namespace klee
}  // namespace axom
