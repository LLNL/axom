// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SLAM_POLYVALUE_HPP
#define SLAM_POLYVALUE_HPP

#include <memory>

namespace axom
{
namespace slam
{
namespace detail
{
template <typename Base>
inline std::unique_ptr<Base> construct_nullptr(const Base&)
{
  return std::unique_ptr<Base> {};
}

template <typename T, typename Base>
struct copy_construct_impl
{
  static std::unique_ptr<Base> construct(const Base& other)
  {
    return std::unique_ptr<Base>(new T(static_cast<const T&>(other)));
  }
};

}  // namespace detail

template <typename T>
struct PolyValue
{
public:
  template <typename U>
  PolyValue(std::remove_pointer_t<U> value)
    : m_poly(new U(std::move(value)))
    , m_copier(&detail::copy_construct_impl<U, T>::construct)
  {
    static_assert(!std::is_abstract<U>::value,
                  "Cannot construct with abstract type.");
  }

  template <typename U>
  [[deprecated]] PolyValue(const U* value) : PolyValue()
  {
    static_assert(!std::is_abstract<U>::value,
                  "Cannot construct with abstract type.");
    if(value)
    {
      m_poly.reset(new U(*value));
      m_copier = &detail::copy_construct_impl<U, T>::construct;
    }
  }

  PolyValue(std::nullptr_t = nullptr) : m_copier(&detail::construct_nullptr<T>)
  { }

  PolyValue(const PolyValue& other)
    : m_poly((other.m_copier)(*(other.m_poly)))
    , m_copier(other.m_copier)
  { }

  PolyValue& operator=(const PolyValue& other)
  {
    if(this != &other)
    {
      m_poly = (other.m_copier)(*(other.m_poly));
      m_copier = other.m_copier;
    }
    return *this;
  }

  PolyValue(PolyValue&&) = default;
  PolyValue& operator=(PolyValue&&) = default;

  ~PolyValue() = default;

  const T* operator->() const { return m_poly.get(); }
  T* operator->() { return m_poly.get(); }

  const T& operator*() const { return *m_poly; }
  T& operator*() { return *m_poly; }

  const T* get() const { return m_poly.get(); }
  T* get() { return m_poly.get(); }

private:
  std::unique_ptr<T> m_poly;
  std::unique_ptr<T> (*m_copier)(const T&);
};

}  // namespace slam
}  // namespace axom

#endif  // SLAM_POLYVALUE_HPP
