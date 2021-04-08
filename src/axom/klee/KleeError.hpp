// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_KLEEERROR_HPP
#define AXOM_KLEEERROR_HPP

#include <exception>
#include <string>
#include <utility>

namespace axom
{
namespace klee
{
class KleeError : public std::exception
{
public:
  explicit KleeError(std::string message) : m_message {std::move(message)} { }

  const char *what() const noexcept override { return m_message.data(); }

private:
  std::string m_message;
};

}  // namespace klee
}  // namespace axom

#endif  //AXOM_KLEEERROR_HPP
