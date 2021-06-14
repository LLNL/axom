// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_KLEEERROR_HPP
#define AXOM_KLEEERROR_HPP

#include <exception>
#include <vector>

#include "axom/inlet/inlet_utils.hpp"

namespace axom
{
namespace klee
{
class KleeError : public std::exception
{
public:
  explicit KleeError(const inlet::VerificationError &error);

  explicit KleeError(const std::vector<inlet::VerificationError> &errors);

  const char *what() const noexcept override;

  const std::vector<inlet::VerificationError> &getErrors() const
  {
    return m_errors;
  }

private:
  std::vector<inlet::VerificationError> m_errors;
};

}  // namespace klee
}  // namespace axom

#endif  //AXOM_KLEEERROR_HPP
