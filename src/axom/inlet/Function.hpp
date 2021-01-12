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
#include <functional>
#include <tuple>
#include <type_traits>

#include "fmt/fmt.hpp"

#include "axom/sidre.hpp"

#include "axom/inlet/Field.hpp"
#include "axom/inlet/InletVector.hpp"
#include "axom/inlet/Verifiable.hpp"
#include "axom/inlet/inlet_utils.hpp"

#include "axom/sidre.hpp"

namespace axom
{
namespace inlet
{
/*!
 *******************************************************************************
 * \brief The tags used to describe function signatures in the input file
 * 
 * \note Additions to this enumeration should be propagated to FunctionType, LuaReader
 * and func_signature_tuples defined below (the mapping from enum to types)
 * 
 * \note Vec3D corresponds to a three-dimensional vector, Double corresponds to
 * a floating-point scalar
 * 
 * \note A two-dimensional vector was intentionally excluded for simplicity as
 * a 3D vector can be used in its place (with the third component being empty/zero)
 *******************************************************************************
 */
enum class FunctionTag
{
  Vec3D,
  Double
};

/*!
 *******************************************************************************
 * \brief The types used to describe function signatures in the input file
 * 
 * \note These are aliases intended to improve readability
 *******************************************************************************
 */
struct FunctionType
{
  using Vec3D = InletVector;
  using Double = double;
};

namespace detail
{
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
 * \class cleanup_function_signature
 *
 * \brief Takes a Ret(Args...) signature and applies cvref qualifiers to applicable
 * arguments
 * 
 * \tparam Ret The function's return type
 * \tparam Args... The function's arguments
 * 
 * Designed to be used with user-facing retrieval functions so the user is not 
 * required to specify the qualifiers in something like a get<std::function<Ret(Args...)>
 *******************************************************************************
 */
template <typename>
struct cleanup_function_signature;

template <typename Ret, typename... Args>
struct cleanup_function_signature<Ret(Args...)>
{
  using type = Ret(typename inlet_function_arg_type<Args>::type...);
};

// The permutation metaprogramming was assisted by the following:
// https://nolte.dev/jekyll/update/2019/06/06/fun-with-tuples.html
// https://stackoverflow.com/questions/55005063/create-a-type-list-combination-of-types-in-c

/*!
 *******************************************************************************
 * \brief Adds a single type to a single tuple type via prepending
 *
 * \tparam T The type to add to the tuple, e.g., <A>
 * \tparam Ts... The types in the existing tuple, e.g., <B, C, D>
 * 
 * In the above example the result type is std::tuple<A, B, C, D>
 *******************************************************************************
 */
template <typename T, typename... Ts>
struct tuple_concat;

template <typename T, typename... Ts>
struct tuple_concat<T, std::tuple<Ts...>>
{
  // Expand the parameter pack within the single tuple alias decl
  using type = std::tuple<T, Ts...>;
};

/*!
 *******************************************************************************
 * \brief Adds a single type to a parameter pack of tuples
 *
 * \tparam T The type to add to the tuple, e.g., <A>
 * \tparam Ts... The types of the tuples to add to, e.g., <std::tuple<B, C>, std::tuple<D, E>>
 * 
 * In the above example the result type is std::tuple<std::tuple<A, B, C>, std::tuple<A, D, E>>
 * 
 * The Python-like pseudocode for this function is roughly:
 * \code{.py}
 * def multi_tuple_concat(T, Ts):
 *   return [T + T_ for T_ in Ts]
 * \endcode
 *******************************************************************************
 */
template <typename T, typename... Ts>
struct multi_tuple_concat
{
  // Expand the parameter pack across calls to tuple_concat so T is prepended
  // to each tuple
  using type = std::tuple<typename tuple_concat<T, Ts>::type...>;
};

/*!
 *******************************************************************************
 * \brief Recursive helper struct for enumerating permutations of a list of types
 *
 * \tparam N The zero-indexed length of the permutations to produce ("stack height")
 * \tparam Ts... The original (1-dimensional) set of types to permute, e.g. <A, B, C>
 * \tparam Us... The current "result" tuple of tuples, each member U in Us will
 * be a tuple whose size is the original requested permutation length minus the
 * current stack height N
 * 
 * \note This function works by modifying the "result" type.  The multi_tuple_concat
 * is called to add each original type T to each U in Us.  The resulting tuple
 * of tuples from the concatenation of each T is then concatenated with std::tuple_cat.
 * For example, with Ts = <A, B, C> and Us = std::tuple<std::tuple<A>, std::tuple<B>, std::tuple<C>>,
 * the full tuple_cat call is std::tuple_cat(std::tuple<std::tuple<A, A>, std::tuple<A, B>, 
 * std::tuple<A, C>>, std::tuple<std::tuple<B, A>, std::tuple<B, B>, std::tuple<B, C>>, ...)
 * 
 * The Us pack is expanded in the multi_tuple_concat "call" whereas the Ts pack
 * is expanded as part of the std::tuple_cat "call"
 * 
 * The Python-like pseudocode for this function is roughly:
 * \code{.py}
 * def permutation_helper(Ts, Us, N):
 *   if N == 0:
 *     return Us
 *   else:
 *     result = []
 *     for T in Ts:
 *       result.append(permute(Ts, [T + U for U in Us], N - 1))
 *     return result
 * \endcode
 *******************************************************************************
 */
template <std::size_t N>
struct permutation_helper
{
  // Implemented as a function to allow for two parameter packs
  template <typename... Ts, typename... Us>
  static auto get(std::tuple<Ts...> t, std::tuple<Us...> u)
    -> decltype(permutation_helper<N - 1>::get(
      t,
      std::tuple_cat(typename multi_tuple_concat<Ts, Us...>::type()...)));
};

// Base case -> N = 0 -> Us... is the result, so just "return" it
template <>
struct permutation_helper<0u>
{
  template <typename... Ts, typename... Us>
  static auto get(std::tuple<Ts...> t, std::tuple<Us...> u) -> decltype(u);
};

/*!
 *******************************************************************************
 * \brief Entry point for retrieving a permutation of types
 *
 * \tparam N The (1-indexed) permutation length
 * \tparam Ts... The types to permute
 * 
 * \note This initializes the result by creating a tuple of tuples of each type
 * in the list, e.g., for Ts = <A, B, C>, the second declval will be 
 * std::tuple<std::tuple<A>, std::tuple<B>, std::tuple<C>> (the first will be just
 * std::tuple<A, B, C>)
 * 
 * The Python-like pseudocode for this function is roughly:
 * \code{.py}
 * def type_permutations(Ts, N):
 *   return permutation_helper(Ts, [(T) for T in Ts])
 * \endcode
 *******************************************************************************
 */
template <std::size_t N, typename... Ts>
using type_permutations = decltype(permutation_helper<N - 1>::get(
  std::declval<std::tuple<Ts...>>(),
  std::declval<std::tuple<std::tuple<Ts>...>>()));

/*!
 *******************************************************************************
 * \brief Converts a tuple into a function signature
 *
 * \tparam Ret The function return type
 * \tparam Args... The function's argument types
 * 
 * \note This would be called as tuple_to_inlet_signature<Ret, Arg1, Arg2, ..., etc>
 * to convert to Ret(Arg1, Arg2, ...) after applying the cvref transformation
 * by inlet_function_arg_type<T>
 *******************************************************************************
 */
template <typename... Ts>
struct tuple_to_inlet_signature;

template <typename Ret, typename... Args>
struct tuple_to_inlet_signature<std::tuple<Ret, Args...>>
{
  using type = Ret(typename inlet_function_arg_type<Args>::type...);
};

}  // end namespace detail

/*!
 *******************************************************************************
 * \class FunctionWrapper
 *
 * \brief A sum type for callables with arbitrary signature
 * \tparam FunctionTags The types of supported functions
 * 
 * This provides an interface not templated on a specific function signature
 * for uniform retrieval through the Reader interface
 *******************************************************************************
 */
template <typename... FunctionTags>
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

  FunctionWrapper() = default;
  FunctionWrapper(FunctionWrapper&& other) = default;

  /*!
   *******************************************************************************
   * \brief Calls the function
   * 
   * \param [in] args The parameter pack for the function's arguments
   * \tparam Ret The user-specified return type, needed to fully disambiguate the
   * function to call
   * \tparam Args The types of the user-specified arguments, deduced automatically
   * 
   * \return The function's result
   *******************************************************************************
   */
  template <typename Ret, typename... Args>
  Ret call(Args&&... args) const
  {
    const auto& ptr = std::get<std::unique_ptr<
      std::function<Ret(typename detail::inlet_function_arg_type<Args>::type...)>>>(
      m_funcs);
    SLIC_ERROR_IF(
      !m_function_valid || !ptr || !(*ptr),
      fmt::format("[Inlet] Function '{0}' with requested type does not exist",
                  m_name));

    const auto& func = *ptr;
    return func(
      std::forward<typename detail::inlet_function_arg_type<Args>::type>(args)...);
  }

  template <typename FuncType>
  std::function<FuncType> get() const
  {
    const auto& ptr = std::get<std::unique_ptr<
      std::function<typename detail::cleanup_function_signature<FuncType>::type>>>(
      m_funcs);
    SLIC_ERROR_IF(
      !ptr,
      fmt::format("[Inlet] Function '{0}' with requested type does not exist",
                  m_name));
    return *ptr;
  }

  /*!
   *******************************************************************************
   * \brief Checks whether the function exists
   *******************************************************************************
   */
  explicit operator bool() const { return m_function_valid; }

  /*!
   *******************************************************************************
   * \brief Sets the function's name
   * 
   * \note Needs to be separate from constructor to allow the compiler to deduce
   * template arguments correctly
   *******************************************************************************
   */
  void setName(std::string&& name) { m_name = std::move(name); }

private:
  // This is on the heap to reduce size - each pointer is only 8 bytes vs 32 bytes
  // for a std::function, and it is guaranteed that only one of the pointers will
  // actually point to something
  std::tuple<std::unique_ptr<std::function<FunctionTags>>...> m_funcs;
  bool m_function_valid = false;
  std::string m_name;
};

namespace detail
{
/*!
 *******************************************************************************
 * \brief "Unwraps" a tuple of tuples, converts the inner tuples to function
 * signatures, and defines a FunctionWrapper templated on those signatures
 *
 * \tparam Tuples The tuple of tuples to expand, e.g.,
 * std::tuple<std::tuple<A, B>, std::tuple<A, B, C>> to use the signatures
 * A(B) and A(B, C) (subject to cvref qualifiers added by inlet_function_arg_type)
 *******************************************************************************
 */
template <typename... Tuples>
struct tuples_to_wrapper;

template <typename... Tuples>
struct tuples_to_wrapper<std::tuple<Tuples...>>
{
  using type =
    FunctionWrapper<typename tuple_to_inlet_signature<Tuples>::type...>;
};

/*!
 *****************************************************************************
 * \brief Generates the permutations of a parameter pack up to a certain
 * length, including all shorter permutations
 * 
 * \tparam N The length of the longest permutation, will also generate
 * permutations of length N-1, N-2, ..., 1
 * \tparam Ts... The parameter pack of types to enumerate
 *****************************************************************************
 */
template <std::size_t N, typename... Ts>
struct arg_tuples
{
  using type = decltype(
    std::tuple_cat(std::declval<type_permutations<N, Ts...>>(),
                   std::declval<typename arg_tuples<N - 1, Ts...>::type>()));
};

// Base case - empty tuple
template <typename... Ts>
struct arg_tuples<0u, Ts...>
{
  using type = std::tuple<>;
};

/*!
 *****************************************************************************
 * \brief The maximum number of user-specified arguments to a function
 *
 * \note Be very cautious when increasing this, as it will result in exponential
 * function generation - specifically, m^n where n is this MAX_NUM_ARGS
 * and m is the number of elements in the FunctionTag enumeration
 *****************************************************************************
 */
static constexpr std::size_t MAX_NUM_ARGS = 2u;

// Get the permutations of all possible signatures
// Add one as return types also need to be permuted
using func_signature_tuples =
  arg_tuples<MAX_NUM_ARGS + 1, FunctionType::Vec3D, FunctionType::Double>::type;

using BasicFunctionWrapper = tuples_to_wrapper<func_signature_tuples>::type;

}  // end namespace detail

using FunctionVariant = detail::BasicFunctionWrapper;

class Function : public Verifiable<Function>
{
public:
  Function(axom::sidre::Group* sidreGroup,
           axom::sidre::Group* root,
           FunctionVariant&& func,
           bool docEnabled = true)
    : m_sidreGroup(sidreGroup)
    , m_sidreRootGroup(root)
    , m_docEnabled(docEnabled)
    , m_func(std::move(func))
  {
    m_func.setName(name());
  }

  /*!
   *****************************************************************************
   * \brief Returns a function of requested type
   * 
   * \return The value
   * \tparam T The type to retrieve
   *****************************************************************************
   */
  template <typename FuncType>
  std::function<FuncType> get() const
  {
    return m_func.get<FuncType>();
  }

  template <typename Ret, typename... Args>
  Ret call(Args&&... args) const
  {
    return m_func.call<Ret>(std::forward<Args>(args)...);
  }

  /*!
   *******************************************************************************
   * \brief Checks whether the function exists
   *******************************************************************************
   */
  operator bool() const { return static_cast<bool>(m_func); }

  /*!
   *****************************************************************************
   * \return The full name for this Function.
   *****************************************************************************
  */
  std::string name() const;

  /*!
   *****************************************************************************
   * \brief This will be called by Inlet::verify to verify the contents of this
   *  Function
   *****************************************************************************
   */
  bool verify() const;

  /*!
   *****************************************************************************
   * \brief Set the required status of this Table.
   *
   * Set whether this Table is required, or not, to be in the input file.
   * The default behavior is to not be required.
   *
   * \param [in] isRequired Boolean value of whether Table is required
   *
   * \return Reference to this instance of Table
   *****************************************************************************
   */
  Function& required(bool isRequired = true);

  /*!
   *****************************************************************************
   * \brief Return the required status of this Table.
   *
   * Return that this Function is required, or not, to be in the input file.
   * The default behavior is to not be required.
   *
   * \return Boolean value of whether this Function is required
   *****************************************************************************
   */
  bool isRequired() const;

  /*!
   *****************************************************************************
   * \brief Registers the function object that will verify this function
   * during the verification stage.
   * 
   * \param [in] The function object that will be called by Table::verify().
   *****************************************************************************
  */
  Function& registerVerifier(std::function<bool(const Function&)> lambda);

private:
  // This function's sidre group
  axom::sidre::Group* m_sidreGroup = nullptr;
  axom::sidre::Group* m_sidreRootGroup = nullptr;
  bool m_docEnabled;
  std::function<bool(const Function&)> m_verifier;
  FunctionVariant m_func;
};

}  // end namespace inlet
}  // end namespace axom

#endif  // INLET_FUNCTION_HPP
