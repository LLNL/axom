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
#include <type_traits>

#include "fmt/fmt.hpp"

#include "axom/sidre.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include "axom/inlet/Field.hpp"
#include "axom/inlet/Verifiable.hpp"

#include "axom/sidre.hpp"

namespace axom
{
namespace inlet
{
enum class InletFunctionType
{
  Vec2D,
  Vec3D,
  Double
};

namespace detail
{
template <InletFunctionType>
struct inlet_function_type;

template <>
struct inlet_function_type<InletFunctionType::Vec2D>
{
  using type = primal::Vector2D;
};

template <>
struct inlet_function_type<InletFunctionType::Vec3D>
{
  using type = primal::Vector3D;
};

template <>
struct inlet_function_type<InletFunctionType::Double>
{
  using type = double;
};

// Makes non-primitive parameters (not bool/int/double) const references
template <typename Arg>
struct inlet_function_arg_type
{
  using type = typename std::conditional<
    std::is_arithmetic<Arg>::value,
    Arg,
    typename std::add_lvalue_reference<typename std::add_const<Arg>::type>::type>::type;
};

template <InletFunctionType Result, InletFunctionType Arg>
struct inlet_function_signature
{
  using type = typename inlet_function_type<Result>::type(
    typename inlet_function_arg_type<typename inlet_function_type<Arg>::type>::type);
};

using Vec3_Double =
  inlet_function_signature<InletFunctionType::Double, InletFunctionType::Vec3D>::type;

template <typename... FunctionTypes>
class FunctionWrapper
{
public:
  template <typename FuncType>
  FunctionWrapper(std::function<FuncType>&& func)
  {
    m_function_valid = static_cast<bool>(func);
    std::get<std::function<FuncType>>(m_funcs) = std::move(func);
  }

  FunctionWrapper(FunctionWrapper&& other) = default;

  template <typename Ret,
            typename... Args,
            typename SFINAE = decltype(std::get<std::function<Ret(Args...)>>(
              std::declval<std::tuple<std::function<FunctionTypes>...>>()))>
  Ret call(Args const&... args) const
  {
    const auto& func = std::get<std::function<Ret(Args const&...)>>(m_funcs);
    SLIC_ERROR_IF(!func, "[Inlet] Function with requested type does not exist");
    return func(std::forward<Args const&>(args)...);
  }

  operator bool() const { return m_function_valid; }

private:
  std::tuple<std::function<FunctionTypes>...> m_funcs;
  bool m_function_valid = false;
};

}  // end namespace detail

using InletFunctionWrapper = detail::FunctionWrapper<detail::Vec3_Double>;

class Function : public Verifiable<Function>
{
public:
  Function(axom::sidre::Group* sidreGroup,
           axom::sidre::Group* root,
           InletFunctionWrapper&& func,
           bool docEnabled = true)
    : m_sidreGroup(sidreGroup)
    , m_sidreRootGroup(root)
    , m_docEnabled(docEnabled)
    , m_functions(std::move(func))
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
  InletFunctionWrapper m_functions;
};

}  // end namespace inlet
}  // end namespace axom

#endif  // INLET_FUNCTION_HPP
