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

#include <memory>
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
#include "axom/inlet/inlet_utils.hpp"

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
/*!
 *******************************************************************************
 * \class inlet_function_type
 *
 * \brief A type trait for mapping values of the InletFunctionType enum to types
 * \tparam InletFunctionType The value of the enum, a non-type template parameter
 * 
 * Specializations must be defined for each member of InletFunctionType
 *******************************************************************************
 */
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

/*!
 *******************************************************************************
 * \class inlet_function_arg_type
 *
 * \brief A type trait for modifying function argument types to enforce const
 * correctness and to avoid copies
 * \tparam Arg The argument type
 * 
 * Maps a type "Arg" to "const Arg&" if that type does not satisfy std::is_arithmetic,
 * which matches Inlet primitive types bool, int, and double, but not vectors or
 * std::strings - though only vectors and doubles are currently supported as arg
 * types
 *******************************************************************************
 */
template <typename Arg>
struct inlet_function_arg_type
{
  using type = typename std::conditional<
    std::is_arithmetic<Arg>::value,
    Arg,
    typename std::add_lvalue_reference<typename std::add_const<Arg>::type>::type>::type;
};

/*!
 *******************************************************************************
 * \class inlet_function_signature
 *
 * \brief A type trait for building a function signature type
 * \tparam Result The function's return type
 * \tparam Args The parameter pack for the function's arguments - currently only
 * single-argument functions are supported
 * 
 * Creates a function signature usable as the template parameter to std::function
 * given a set of InletFunctionTypes
 *******************************************************************************
 */
template <InletFunctionType Result, InletFunctionType... Args>
struct inlet_function_signature
{
  using type = typename inlet_function_type<Result>::type(
    typename inlet_function_arg_type<typename inlet_function_type<Args...>::type>::type);
};

}  // end namespace detail

/*!
 *******************************************************************************
 * \class FunctionWrapper
 *
 * \brief A sum type for callables with arbitrary signature
 * \tparam FunctionTypes The types of supported functions - can be generated
 * with detail::inlet_function_signature
 * 
 * This provides an interface not templated on a specific function signature
 * for uniform retrieval through the Reader interface
 *******************************************************************************
 */
template <typename... FunctionTypes>
class FunctionWrapper
{
public:
  /*!
   *******************************************************************************
   * \brief Primary constructor, initializes a member of the "variant"
   * 
   * \param [in] func The function to initialize with
   * \tparam FuncType The function's signature
   * 
   * \note "Empty" functions are allowable and result in the object comparing
   * false when converted to bool
   *******************************************************************************
   */
  template <typename FuncType>
  FunctionWrapper(std::function<FuncType>&& func)
  {
    m_function_valid = static_cast<bool>(func);
    std::get<std::unique_ptr<std::function<FuncType>>>(m_funcs) =
      cpp11_compat::make_unique<std::function<FuncType>>(std::move(func));
  }

  /*!
   *******************************************************************************
   * \brief Default constructor
   * 
   * Should never be called as it is only required when a return occurs after a 
   * SLIC_ERROR which does not indicate termination to the compiler.  If exceptions
   * are ever switched to this function and its calls should be removed.
   *******************************************************************************
   */
  FunctionWrapper() = default;
  FunctionWrapper(FunctionWrapper&& other) = default;

  /*!
   *******************************************************************************
   * \brief Calls the function
   * 
   * \param [in] args The parameter pack for the function's arguments
   * \tparam Ret The user-specified return type, needed to fully disambiguate the
   * function to call
   * \tparam Args The types of the user-specified arguments, deduced
   * 
   * \return The function's result
   *******************************************************************************
   */
  template <typename Ret, typename... Args>
  Ret call(Args&&... args) const
  {
    using ArgTypes = typename detail::inlet_function_arg_type<Args...>::type;
    const auto& func =
      *std::get<std::unique_ptr<std::function<Ret(ArgTypes)>>>(m_funcs);
    SLIC_ERROR_IF(!func || !m_function_valid,
                  "[Inlet] Function with requested type does not exist");
    return func(std::forward<ArgTypes>(args)...);
  }

  /*!
   *******************************************************************************
   * \brief Checks whether the function exists
   *******************************************************************************
   */
  operator bool() const { return m_function_valid; }

private:
  // This is on the heap to reduce size - each pointer is only 8 bytes vs 32 bytes
  // for a std::function, and it is guaranteed that only one of the pointers will
  // actually point to something
  std::tuple<std::unique_ptr<std::function<FunctionTypes>>...> m_funcs;
  bool m_function_valid = false;
};

// The set of supported type signatures to be used in the "variant" type
namespace detail
{
// This could definitely be done with more metaprogramming, which may be the only option
// as the number of argumetns increases

using Vec2_Vec3 =
  inlet_function_signature<InletFunctionType::Vec2D, InletFunctionType::Vec3D>::type;

using Vec2_Vec2 =
  inlet_function_signature<InletFunctionType::Vec2D, InletFunctionType::Vec2D>::type;

using Vec2_Double =
  inlet_function_signature<InletFunctionType::Vec2D, InletFunctionType::Double>::type;

using Vec3_Vec3 =
  inlet_function_signature<InletFunctionType::Vec3D, InletFunctionType::Vec3D>::type;

using Vec3_Vec2 =
  inlet_function_signature<InletFunctionType::Vec3D, InletFunctionType::Vec2D>::type;

using Vec3_Double =
  inlet_function_signature<InletFunctionType::Vec3D, InletFunctionType::Double>::type;

using Double_Vec3 =
  inlet_function_signature<InletFunctionType::Double, InletFunctionType::Vec3D>::type;

using Double_Vec2 =
  inlet_function_signature<InletFunctionType::Double, InletFunctionType::Vec2D>::type;

using Double_Double =
  inlet_function_signature<InletFunctionType::Double, InletFunctionType::Double>::type;

using BasicFunctionWrapper = FunctionWrapper<Vec2_Vec3,
                                             Vec2_Vec2,
                                             Vec2_Double,
                                             Vec3_Vec3,
                                             Vec3_Vec2,
                                             Vec3_Double,
                                             Double_Vec3,
                                             Double_Vec2,
                                             Double_Double>;

}  // end namespace detail

using InletFunctionWrapper = detail::BasicFunctionWrapper;

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
