#ifndef SINA_CPPBRIDGE_HPP
#define SINA_CPPBRIDGE_HPP

/**
 * @file
 *
 * This file contains functions which can go away once this library is
 * migrated to C++14 or later.
 */

#include <memory>
#include <utility>

namespace sina { namespace internal {

/**
 * Make a unique pointer pointing to a new object.
 *
 * Delete when switch to C++14 in favor of std::make_unique.
 *
 * @tparam T the type of the object
 * @tparam ParamTypes the types of the parameters
 * @param params the parameters to the constructor
 * @return a unique_ptr to the newly-constructed object
 */
template<typename T, typename... ParamTypes>
std::unique_ptr<T> make_unique(ParamTypes &&... params) {
    return std::unique_ptr<T>(new T(std::forward<ParamTypes>(params)...));
};

}} // namespaces

#endif //SINA_CPPBRIDGE_HPP
