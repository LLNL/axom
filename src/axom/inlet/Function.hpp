// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Function.hpp
 *
 * \brief This file contains the class definition of Inlet's Function class.
 *******************************************************************************
 */

#ifndef INLET_FUNCTION_HPP
#define INLET_FUNCTION_HPP

// #include <memory>
// #include <string>
#include <functional>
// #include <unordered_map>
#include <tuple>
// #include <type_traits>

#include "fmt/fmt.hpp"

#include "axom/sidre.hpp"

#include "axom/inlet/Field.hpp"
#include "axom/inlet/Reader.hpp"
#include "axom/inlet/Verifiable.hpp"

#include "axom/sidre.hpp"

namespace axom
{
namespace inlet
{

namespace detail
{
  using Vec2RetDouble = double(const primal::Vector3D&);

  template <typename... FunctionTypes>
  class FunctionWrapper
  {
    public:
    template <typename Ret, typename... Args>
    Ret operator()(Args&&... args) const
    {
      const auto& func = std::get<std::function<Ret(Args...)>>(m_funcs);
      SLIC_ERROR_IF(!func, "[Inlet] Function with requested type does not exist");
      return func(std::forward<Args>(args)...);
    }
    private:
      std::tuple<std::function<FunctionTypes>...> m_funcs;
  };

  using InletFunctionWrapper = FunctionWrapper<Vec2RetDouble>;
} // end namespace detail

class Function : public Verifiable<Function>
{
  public:
  Function(axom::sidre::Group* sidreGroup,
        axom::sidre::Group* root,
        bool docEnabled = true)
    : m_sidreGroup(sidreGroup)
    , m_sidreRootGroup(root)
    , m_docEnabled(docEnabled)
  { }

  // InletType type() const { return InletType::Function; }

  /*!
   *****************************************************************************
   * \brief Returns a function of requested type
   * 
   * \return The value
   * \tparam T The type to retrieve
   *****************************************************************************
   */
  template <typename FuncType>
  std::function<FuncType> get() const;

  private:
  // This function's sidre group
  axom::sidre::Group* m_sidreGroup = nullptr;
  axom::sidre::Group* m_sidreRootGroup = nullptr;
  bool m_docEnabled;
  std::function<bool(const Function&)> m_verifier;
  detail::InletFunctionWrapper m_functions;
};

}  // end namespace inlet
}  // end namespace axom

#endif // INLET_FUNCTION_HPP
