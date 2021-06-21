// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
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
 * and func_signature_lists defined below (the mapping from enum to types)
 * 
 * \note Vector corresponds to a vector with max dimension of three, Double corresponds to
 * a floating-point scalar
 * 
 * \note A two-dimensional vector was intentionally excluded for simplicity as
 * a 3D vector can be used in its place (with the third component being empty/zero)
 *******************************************************************************
 */
enum class FunctionTag
{
  Vector,
  Double,
  Void,
  String
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
  using Vector = InletVector;
  using Double = double;
  using Void = void;
  using String = std::string;
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

/*!
 *******************************************************************************
 * \brief A lightweight type used to store a list of types
 * \tparam Ts The parameter pack containing the types to store
 *******************************************************************************
 */
template <typename... Ts>
struct TypeList
{ };

/*!
 *******************************************************************************
 * \brief Helper structure used to concatenate parameter packs (lists) of TypeLists
 * \tparam TypeLists The TypeLists to concatenate (in order)
 *******************************************************************************
 */
template <typename... TypeLists>
struct Concat;

template <>
struct Concat<>
{
  using type = TypeList<>;
};

template <typename... Types>
struct Concat<TypeList<Types...>>
{
  using type = TypeList<Types...>;
};

template <typename... First, typename... Second>
struct Concat<TypeList<First...>, TypeList<Second...>>
{
  using type = TypeList<First..., Second...>;
};

template <typename... First, typename... Second, typename... Tail>
struct Concat<TypeList<First...>, TypeList<Second...>, Tail...>
{
  // Recursive call used for n > 2 TypeLists - recursive calls will always
  // be to n - 2 and thus always terminate in one of the above base cases
  using type = typename Concat<TypeList<First..., Second...>, Tail...>::type;
};

/*!
 *******************************************************************************
 * \brief Implements concatenation functionality similar to std::tuple_cat
 * \tparam TypeLists The TypeLists to concatenate (in order)
 * \note This function should only be used in unevaluated contexts (e.g., within
 * a \p decltype expression)
 *******************************************************************************
 */
template <typename... TypeLists>
typename Concat<TypeLists...>::type type_list_cat(TypeLists&&...)
{
  return {};
}

// The permutation metaprogramming was assisted by the following:
// https://nolte.dev/jekyll/update/2019/06/06/fun-with-tuples.html
// https://stackoverflow.com/questions/55005063/create-a-type-list-combination-of-types-in-c

/*!
 *******************************************************************************
 * \brief Adds a single type to a single TypeList via prepending
 *
 * \tparam T The type to add to the TypeList, e.g., <A>
 * \tparam Ts... The types in the existing TypeList, e.g., <B, C, D>
 * 
 * In the above example the result type is TypeList<A, B, C, D>
 *******************************************************************************
 */
template <typename T, typename... Ts>
struct list_prepend;

template <typename T, typename... Ts>
struct list_prepend<T, TypeList<Ts...>>
{
  // Expand the parameter pack within the single list alias decl
  using type = TypeList<T, Ts...>;
};

/*!
 *******************************************************************************
 * \brief Adds a single type to a parameter pack of TypeLists
 *
 * \tparam T The type to add to the TypeLists, e.g., <A>
 * \tparam Ts... The types of the TypeLists to add to, e.g., <TypeList<B, C>, TypeList<D, E>>
 * 
 * In the above example the result type is TypeList<TypeList<A, B, C>, TypeList<A, D, E>>
 * 
 * The Python-like pseudocode for this function is roughly:
 * \code{.py}
 * def multi_list_prepend(T, Ts):
 *   return [T + T_ for T_ in Ts]
 * \endcode
 *******************************************************************************
 */
template <typename T, typename... Ts>
struct multi_list_prepend
{
  // Expand the parameter pack across calls to list_prepend so T is prepended
  // to each TypeList
  using type = TypeList<typename list_prepend<T, Ts>::type...>;
};

/*!
 *******************************************************************************
 * \brief Recursive helper struct for enumerating permutations of a list of types
 *
 * \tparam N The zero-indexed length of the permutations to produce ("stack height")
 * \tparam Ts... The original (1-dimensional) set of types to permute, e.g. <A, B, C>
 * \tparam Us... The current "result" list of lists, each member U in Us will
 * be a list whose size is the original requested permutation length minus the
 * current stack height N
 * 
 * \note This function works by modifying the "result" type.  The multi_list_prepend
 * is called to add each original type T to each U in Us.  The resulting list
 * of lists from the concatenation of each T is then concatenated with type_list_cat.
 * For example, with Ts = <A, B, C> and Us = TypeList<TypeList<A>, TypeList<B>, TypeList<C>>,
 * the full type_list_cat call is type_list_cat(TypeList<TypeList<A, A>, TypeList<A, B>, 
 * TypeList<A, C>>, TypeList<TypeList<B, A>, TypeList<B, B>, TypeList<B, C>>, ...)
 * 
 * The Us pack is expanded in the multi_list_prepend "call" whereas the Ts pack
 * is expanded as part of the type_list_cat "call"
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
  static auto get(TypeList<Ts...> t, TypeList<Us...> u)
    -> decltype(permutation_helper<N - 1>::get(
      t,
      type_list_cat(typename multi_list_prepend<Ts, Us...>::type()...)));
};

// Base case -> N = 0 -> Us... is the result, so just "return" it
template <>
struct permutation_helper<0u>
{
  template <typename... Ts, typename... Us>
  static auto get(TypeList<Ts...> t, TypeList<Us...> u) -> decltype(u);
};

/*!
 *******************************************************************************
 * \brief Entry point for retrieving a permutation of types
 *
 * \tparam N The (1-indexed) permutation length
 * \tparam Ts... The types to permute
 * 
 * \note This initializes the result by creating a list of lists of each type
 * in the list, e.g., for Ts = <A, B, C>, the second declval will be 
 * TypeList<TypeList<A>, TypeList<B>, TypeList<C>> (the first will be just
 * TypeList<A, B, C>)
 * 
 * The Python-like pseudocode for this function is roughly:
 * \code{.py}
 * def type_permutations(Ts, N):
 *   return permutation_helper(Ts, [(T) for T in Ts])
 * \endcode
 *******************************************************************************
 */
template <std::size_t N, typename... Ts>
using type_permutations = decltype(
  permutation_helper<N - 1>::get(std::declval<TypeList<Ts...>>(),
                                 std::declval<TypeList<TypeList<Ts>...>>()));

/*!
 *******************************************************************************
 * \brief Converts a list into a function signature
 *
 * \tparam Ret The function return type
 * \tparam Args... The function's argument types
 * 
 * \note This would be called as list_to_inlet_signature<Ret, Arg1, Arg2, ..., etc>
 * to convert to Ret(Arg1, Arg2, ...) after applying the cvref transformation
 * by inlet_function_arg_type<T>
 *******************************************************************************
 */
template <typename... Ts>
struct list_to_inlet_signature;

template <typename Ret, typename... Args>
struct list_to_inlet_signature<TypeList<Ret, Args...>>
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
 * \brief "Unwraps" a list of lists, converts the inner lists to function
 * signatures, and defines a FunctionWrapper templated on those signatures
 *
 * \tparam TypeLists The list of lists to expand, e.g.,
 * TypeList<TypeList<A, B>, TypeList<A, B, C>> to use the signatures
 * A(B) and A(B, C) (subject to cvref qualifiers added by inlet_function_arg_type)
 *******************************************************************************
 */
template <typename... TypeLists>
struct lists_to_wrapper;

template <typename... TypeLists>
struct lists_to_wrapper<TypeList<TypeLists...>>
{
  using type =
    FunctionWrapper<typename list_to_inlet_signature<TypeLists>::type...>;
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
struct arg_lists
{
  using type = decltype(
    type_list_cat(std::declval<type_permutations<N, Ts...>>(),
                  std::declval<typename arg_lists<N - 1, Ts...>::type>()));
};

// Base case - empty list
template <typename... Ts>
struct arg_lists<0u, Ts...>
{
  using type = TypeList<>;
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

// Permissible return types
using ret_list =
  TypeList<FunctionType::Void, FunctionType::Vector, FunctionType::Double, FunctionType::String>;

// First, permissible argument types are permuted
using arg_permutations =
  arg_lists<MAX_NUM_ARGS, FunctionType::Vector, FunctionType::Double, FunctionType::String>::type;

// Then return types are prepended to the lists to create lists representing a
// full function signature
using arg_permutations_with_returns = decltype(
  permutation_helper<1u>::get(ret_list {}, std::declval<arg_permutations>()));

// Then function signatures are created representing functions that take no arguments
// but return something
using no_arguments_with_returns =
  decltype(permutation_helper<1u>::get(ret_list {}, TypeList<TypeList<>>()));

// The two sets of full function signatures (with args and without args) are concatenated
using func_signature_lists =
  decltype(type_list_cat(std::declval<arg_permutations_with_returns>(),
                         std::declval<no_arguments_with_returns>()));

using BasicFunctionWrapper = lists_to_wrapper<func_signature_lists>::type;

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
   *****************************************************************************
   * \brief Returns pointer to the Sidre Group class for this Function.
   *
   * Provides access to the Sidre Group class that holds all the stored
   * information for this Function instance.
   *****************************************************************************
   */
  const axom::sidre::Group* sidreGroup() const { return m_sidreGroup; };

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

  bool verify(std::vector<VerificationError>* errors = nullptr) const override;

  Function& required(bool isRequired = true) override;

  bool isRequired() const override;

  using Verifiable<Function>::registerVerifier;

  Function& registerVerifier(Verifier lambda) override;

private:
  // This function's sidre group
  axom::sidre::Group* m_sidreGroup = nullptr;
  axom::sidre::Group* m_sidreRootGroup = nullptr;
  bool m_docEnabled;
  Verifier m_verifier;
  FunctionVariant m_func;
};

}  // end namespace inlet
}  // end namespace axom

#endif  // INLET_FUNCTION_HPP
