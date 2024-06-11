// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_VIEWS_NODE_ARRAY_VIEW_HPP_
#define AXOM_MIR_VIEWS_NODE_ARRAY_VIEW_HPP_

#include "axom/slic/slic.hpp"

#include <conduit/conduit.hpp>
#include <iostream>

namespace axom
{
namespace mir
{
namespace views
{
/*
using Node = ::conduit::Node;
using int8 = ::conduit::int8;
using int16 = ::conduit::int16;
using int32 = ::conduit::int32;
using int64 = ::conduit::int64;
using uint8 = ::conduit::uint8;
using uint16 = ::conduit::uint16;
using uint32 = ::conduit::uint32;
using uint64 = ::conduit::uint64;
using float32 = ::conduit::float32;
using float64 = ::conduit::float64;
using index_t = ::conduit::index_t;
using DataType = ::conduit::DataType;
*/

namespace detail
{

struct Delimiter {};

/// Used to separate arguments.
constexpr Delimiter ArgumentDelimiter;

template <typename ... Args>
constexpr int encode_types(Args... args)
{
  return (... | args);
}

template <typename ... Args>
constexpr int select_types(Args... args)
{
  return encode_types((1 << args)...);
}

constexpr bool type_selected(int flag, int bit)
{
  return flag & (1 << bit);
}

constexpr int select_all_types()
{
  return -1;
}

constexpr int select_index_types()
{
  return select_types(conduit::DataType::INT32_ID, conduit::DataType::INT64_ID, conduit::DataType::UINT32_ID, conduit::DataType::UINT64_ID);
}

constexpr int select_float_types()
{
  return select_types(conduit::DataType::FLOAT32_ID, conduit::DataType::FLOAT64_ID);
}

//------------------------------------------------------------------------------
// General Node to ArrayView. Handle all types.
//------------------------------------------------------------------------------
template <bool Enabled, typename FuncType> std::enable_if_t<Enabled, void>
Node_to_ArrayView_single_int8(const conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::int8> view(n.as_int8_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_single_int8(const conduit::Node &n, FuncType &&func)
{
  std::cout << "Unsupported int8 node." << std::endl;
}

template <bool Enabled, typename FuncType> std::enable_if_t<Enabled, void>
Node_to_ArrayView_single_int8(conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::int8> view(n.as_int8_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_single_int8(conduit::Node &n, FuncType &&func)
{
  std::cout << "Unsupported int8 node." << std::endl;
}

template <bool Enabled, typename FuncType> std::enable_if_t<Enabled, void>
Node_to_ArrayView_single_int16(const conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::int16> view(n.as_int16_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_single_int16(const conduit::Node &n, FuncType &&func)
{
  std::cout << "Unsupported int16 node." << std::endl;
}

template <bool Enabled, typename FuncType> std::enable_if_t<Enabled, void>
Node_to_ArrayView_single_int16(conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::int16> view(n.as_int16_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_single_int16(conduit::Node &n, FuncType &&func)
{
  std::cout << "Unsupported int16 node." << std::endl;
}

template <bool Enabled, typename FuncType> std::enable_if_t<Enabled, void>
Node_to_ArrayView_single_int32(const conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::int32> view(n.as_int32_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_single_int32(const conduit::Node &n, FuncType &&func)
{
  std::cout << "Unsupported int32 node." << std::endl;
}

template <bool Enabled, typename FuncType> std::enable_if_t<Enabled, void>
Node_to_ArrayView_single_int32(conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::int32> view(n.as_int32_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_single_int32(conduit::Node &n, FuncType &&func)
{
  std::cout << "Unsupported int32 node." << std::endl;
}

template <bool Enabled, typename FuncType> std::enable_if_t<Enabled, void>
Node_to_ArrayView_single_int64(const conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::int64> view(n.as_int64_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_single_int64(const conduit::Node &n, FuncType &&func)
{
  std::cout << "Unsupported int64 node." << std::endl;
}

template <bool Enabled, typename FuncType> std::enable_if_t<Enabled, void>
Node_to_ArrayView_single_int64(conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::int64> view(n.as_int64_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_single_int64(conduit::Node &n, FuncType &&func)
{
  std::cout << "Unsupported int64 node." << std::endl;
}

template <bool Enabled, typename FuncType> std::enable_if_t<Enabled, void>
Node_to_ArrayView_single_uint8(const conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::uint8> view(n.as_uint8_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_single_uint8(const conduit::Node &n, FuncType &&func)
{
  std::cout << "Unsupported uint8 node." << std::endl;
}

template <bool Enabled, typename FuncType> std::enable_if_t<Enabled, void>
Node_to_ArrayView_single_uint8(conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::uint8> view(n.as_uint8_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_single_uint8(conduit::Node &n, FuncType &&func)
{
  std::cout << "Unsupported uint8 node." << std::endl;
}

template <bool Enabled, typename FuncType> std::enable_if_t<Enabled, void>
Node_to_ArrayView_single_uint16(const conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::uint16> view(n.as_uint16_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_single_uint16(const conduit::Node &n, FuncType &&func)
{
  std::cout << "Unsupported uint16 node." << std::endl;
}

template <bool Enabled, typename FuncType> std::enable_if_t<Enabled, void>
Node_to_ArrayView_single_uint16(conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::uint16> view(n.as_uint16_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_single_uint16(conduit::Node &n, FuncType &&func)
{
  std::cout << "Unsupported uint16 node." << std::endl;
}

template <bool Enabled, typename FuncType> std::enable_if_t<Enabled, void>
Node_to_ArrayView_single_uint32(const conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::uint32> view(n.as_uint32_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_single_uint32(const conduit::Node &n, FuncType &&func)
{
  std::cout << "Unsupported uint32 node." << std::endl;
}

template <bool Enabled, typename FuncType> std::enable_if_t<Enabled, void>
Node_to_ArrayView_single_uint32(conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::uint32> view(n.as_uint32_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_single_uint32(conduit::Node &n, FuncType &&func)
{
  std::cout << "Unsupported uint32 node." << std::endl;
}

template <bool Enabled, typename FuncType> std::enable_if_t<Enabled, void>
Node_to_ArrayView_single_uint64(const conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::uint64> view(n.as_uint64_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_single_uint64(const conduit::Node &n, FuncType &&func)
{
  std::cout << "Unsupported uint64 node." << std::endl;
}

template <bool Enabled, typename FuncType> std::enable_if_t<Enabled, void>
Node_to_ArrayView_single_uint64(conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::uint64> view(n.as_uint64_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_single_uint64(conduit::Node &n, FuncType &&func)
{
  std::cout << "Unsupported uint64 node." << std::endl;
}

template <bool Enabled, typename FuncType> std::enable_if_t<Enabled, void>
Node_to_ArrayView_single_float32(const conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::float32> view(n.as_float32_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_single_float32(const conduit::Node &n, FuncType &&func)
{
  std::cout << "Unsupported float32 node." << std::endl;
}

template <bool Enabled, typename FuncType> std::enable_if_t<Enabled, void>
Node_to_ArrayView_single_float32(conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::float32> view(n.as_float32_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_single_float32(conduit::Node &n, FuncType &&func)
{
  std::cout << "Unsupported float32 node." << std::endl;
}

template <bool Enabled, typename FuncType> std::enable_if_t<Enabled, void>
Node_to_ArrayView_single_float64(const conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::float64> view(n.as_float64_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_single_float64(const conduit::Node &n, FuncType &&func)
{
  std::cout << "Unsupported float64 node." << std::endl;
}

template <bool Enabled, typename FuncType> std::enable_if_t<Enabled, void>
Node_to_ArrayView_single_float64(conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::float64> view(n.as_float64_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_single_float64(conduit::Node &n, FuncType &&func)
{
  std::cout << "Unsupported float64 node." << std::endl;
}

template <int Types = select_all_types(), typename FuncType>
void Node_to_ArrayView_single(const conduit::Node &n, FuncType &&func)
{
  /* Later, with C++17, we can do this:
    if constexpr (type_selected(Types, conduit::DataType::INT8_ID))
    {
      if(n.dtype().is_int8())
      {
        axom::ArrayView<conduit::int8> view(n.as_int8_ptr(), size);
        func(view);
      }
    }
  */

  if(n.dtype().is_int8())
  {
    Node_to_ArrayView_single_int8<type_selected(Types, conduit::DataType::INT8_ID)>(n, func);
  }
  else if(n.dtype().is_int16())
  {
    Node_to_ArrayView_single_int16<type_selected(Types, conduit::DataType::INT16_ID)>(n, func);
  }
  else if(n.dtype().is_int32())
  {
    Node_to_ArrayView_single_int32<type_selected(Types, conduit::DataType::INT32_ID)>(n, func);
  }
  else if(n.dtype().is_int64())
  {
    Node_to_ArrayView_single_int64<type_selected(Types, conduit::DataType::INT64_ID)>(n, func);
  }
  else if(n.dtype().is_uint8())
  {
    Node_to_ArrayView_single_uint8<type_selected(Types, conduit::DataType::UINT8_ID)>(n, func);
  }
  else if(n.dtype().is_uint16())
  {
    Node_to_ArrayView_single_uint16<type_selected(Types, conduit::DataType::UINT16_ID)>(n, func);
  }
  else if(n.dtype().is_uint32())
  {
    Node_to_ArrayView_single_uint32<type_selected(Types, conduit::DataType::UINT32_ID)>(n, func);
  }
  else if(n.dtype().is_uint64())
  {
    Node_to_ArrayView_single_uint64<type_selected(Types, conduit::DataType::UINT64_ID)>(n, func);
  }
  else if(n.dtype().is_float32())
  {
    Node_to_ArrayView_single_float32<type_selected(Types, conduit::DataType::FLOAT32_ID)>(n, func);
  }
  else if(n.dtype().is_float64())
  {
    Node_to_ArrayView_single_float64<type_selected(Types, conduit::DataType::FLOAT64_ID)>(n, func);
  }
}

template <int Types = select_all_types(), typename FuncType>
void Node_to_ArrayView_single(conduit::Node &n, FuncType &&func)
{
  if(n.dtype().is_int8())
  {
    Node_to_ArrayView_single_int8<type_selected(Types, conduit::DataType::INT8_ID)>(n, func);
  }
  else if(n.dtype().is_int16())
  {
    Node_to_ArrayView_single_int16<type_selected(Types, conduit::DataType::INT16_ID)>(n, func);
  }
  else if(n.dtype().is_int32())
  {
    Node_to_ArrayView_single_int32<type_selected(Types, conduit::DataType::INT32_ID)>(n, func);
  }
  else if(n.dtype().is_int64())
  {
    Node_to_ArrayView_single_int64<type_selected(Types, conduit::DataType::INT64_ID)>(n, func);
  }
  else if(n.dtype().is_uint8())
  {
    Node_to_ArrayView_single_uint8<type_selected(Types, conduit::DataType::UINT8_ID)>(n, func);
  }
  else if(n.dtype().is_uint16())
  {
    Node_to_ArrayView_single_uint16<type_selected(Types, conduit::DataType::UINT16_ID)>(n, func);
  }
  else if(n.dtype().is_uint32())
  {
    Node_to_ArrayView_single_uint32<type_selected(Types, conduit::DataType::UINT32_ID)>(n, func);
  }
  else if(n.dtype().is_uint64())
  {
    Node_to_ArrayView_single_uint64<type_selected(Types, conduit::DataType::UINT64_ID)>(n, func);
  }
  else if(n.dtype().is_float32())
  {
    Node_to_ArrayView_single_float32<type_selected(Types, conduit::DataType::FLOAT32_ID)>(n, func);
  }
  else if(n.dtype().is_float64())
  {
    Node_to_ArrayView_single_float64<type_selected(Types, conduit::DataType::FLOAT64_ID)>(n, func);
  }
}

template <int Types, typename FuncType, typename ... View>
void Node_to_ArrayView_internal(FuncType &&func, Delimiter, View&... views)
{
  func(views...);
}

template <int Types = select_all_types(), typename ... Args>
void Node_to_ArrayView_internal(const conduit::Node &first, Args&&... args)
{
  Node_to_ArrayView_single<Types>(first, [&](auto view)
  {
    Node_to_ArrayView_internal<Types>(args..., view);
  });
}

template <int Types = select_all_types(), typename ... Args>
void Node_to_ArrayView_internal(conduit::Node &first, Args&&... args)
{
  Node_to_ArrayView_single<Types>(first, [&](auto view)
  {
    Node_to_ArrayView_internal<Types>(args..., view);
  });
}

//------------------------------------------------------------------------------
template <bool Enabled, typename FuncType, typename... Args> std::enable_if_t<Enabled, void>
Node_to_ArrayView_same_internal_int8(FuncType &&func, Args&&... args)
{
  func(axom::ArrayView<conduit::int8>(args.as_int8_ptr(), args.dtype().number_of_elements())...);
}

template <bool Enabled, typename FuncType, typename... Args> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_same_internal_int8(FuncType &&func, Args&&... args)
{
}

template <bool Enabled, typename FuncType, typename... Args> std::enable_if_t<Enabled, void>
Node_to_ArrayView_same_internal_int16(FuncType &&func, Args&&... args)
{
  func(axom::ArrayView<conduit::int16>(args.as_int16_ptr(), args.dtype().number_of_elements())...);
}

template <bool Enabled, typename FuncType, typename... Args> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_same_internal_int16(FuncType &&func, Args&&... args)
{
}

template <bool Enabled, typename FuncType, typename... Args> std::enable_if_t<Enabled, void>
Node_to_ArrayView_same_internal_int32(FuncType &&func, Args&&... args)
{
  func(axom::ArrayView<conduit::int32>(args.as_int32_ptr(), args.dtype().number_of_elements())...);
}

template <bool Enabled, typename FuncType, typename... Args> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_same_internal_int32(FuncType &&func, Args&&... args)
{
}

template <bool Enabled, typename FuncType, typename... Args> std::enable_if_t<Enabled, void>
Node_to_ArrayView_same_internal_int64(FuncType &&func, Args&&... args)
{
  func(axom::ArrayView<conduit::int64>(args.as_int64_ptr(), args.dtype().number_of_elements())...);
}

template <bool Enabled, typename FuncType, typename... Args> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_same_internal_int64(FuncType &&func, Args&&... args)
{
}

template <bool Enabled, typename FuncType, typename... Args> std::enable_if_t<Enabled, void>
Node_to_ArrayView_same_internal_uint8(FuncType &&func, Args&&... args)
{
  func(axom::ArrayView<conduit::uint8>(args.as_uint8_ptr(), args.dtype().number_of_elements())...);
}

template <bool Enabled, typename FuncType, typename... Args> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_same_internal_uint8(FuncType &&func, Args&&... args)
{
}

template <bool Enabled, typename FuncType, typename... Args> std::enable_if_t<Enabled, void>
Node_to_ArrayView_same_internal_uint16(FuncType &&func, Args&&... args)
{
  func(axom::ArrayView<conduit::uint16>(args.as_uint16_ptr(), args.dtype().number_of_elements())...);
}

template <bool Enabled, typename FuncType, typename... Args> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_same_internal_uint16(FuncType &&func, Args&&... args)
{
}

template <bool Enabled, typename FuncType, typename... Args> std::enable_if_t<Enabled, void>
Node_to_ArrayView_same_internal_uint32(FuncType &&func, Args&&... args)
{
  func(axom::ArrayView<conduit::uint32>(args.as_uint32_ptr(), args.dtype().number_of_elements())...);
}

template <bool Enabled, typename FuncType, typename... Args> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_same_internal_uint32(FuncType &&func, Args&&... args)
{
}

template <bool Enabled, typename FuncType, typename... Args> std::enable_if_t<Enabled, void>
Node_to_ArrayView_same_internal_uint64(FuncType &&func, Args&&... args)
{
  func(axom::ArrayView<conduit::uint64>(args.as_uint64_ptr(), args.dtype().number_of_elements())...);
}

template <bool Enabled, typename FuncType, typename... Args> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_same_internal_uint64(FuncType &&func, Args&&... args)
{
}

template <bool Enabled, typename FuncType, typename... Args> std::enable_if_t<Enabled, void>
Node_to_ArrayView_same_internal_float32(FuncType &&func, Args&&... args)
{
  func(axom::ArrayView<conduit::float32>(args.as_float32_ptr(), args.dtype().number_of_elements())...);
}

template <bool Enabled, typename FuncType, typename... Args> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_same_internal_float32(FuncType &&func, Args&&... args)
{
}

template <bool Enabled, typename FuncType, typename... Args> std::enable_if_t<Enabled, void>
Node_to_ArrayView_same_internal_float64(FuncType &&func, Args&&... args)
{
  func(axom::ArrayView<conduit::float64>(args.as_float64_ptr(), args.dtype().number_of_elements())...);
}

template <bool Enabled, typename FuncType, typename... Args> std::enable_if_t<!Enabled, void>
Node_to_ArrayView_same_internal_float64(FuncType &&func, Args&&... args)
{
}

template <int Types = select_all_types(), typename FuncType, typename... Args>
void Node_to_ArrayView_same_internal(FuncType &&func, Delimiter, const conduit::Node &first, Args&&... args)
{
  if(first.dtype().is_int8())
  {
    Node_to_ArrayView_same_internal_int8<type_selected(Types, conduit::DataType::INT8_ID)>(func, first, args...);
  }
  else if(first.dtype().is_int16())
  {
    Node_to_ArrayView_same_internal_int16<type_selected(Types, conduit::DataType::INT16_ID)>(func, first, args...);
  }
  else if(first.dtype().is_int32())
  {
    Node_to_ArrayView_same_internal_int32<type_selected(Types, conduit::DataType::INT32_ID)>(func, first, args...);
  }
  else if(first.dtype().is_int64())
  {
    Node_to_ArrayView_same_internal_int64<type_selected(Types, conduit::DataType::INT64_ID)>(func, first, args...);
  }
  else if(first.dtype().is_uint8())
  {
    Node_to_ArrayView_same_internal_uint8<type_selected(Types, conduit::DataType::UINT8_ID)>(func, first, args...);
  }
  else if(first.dtype().is_uint16())
  {
    Node_to_ArrayView_same_internal_uint16<type_selected(Types, conduit::DataType::UINT16_ID)>(func, first, args...);
  }
  else if(first.dtype().is_uint32())
  {
    Node_to_ArrayView_same_internal_uint32<type_selected(Types, conduit::DataType::UINT32_ID)>(func, first, args...);
  }
  else if(first.dtype().is_uint64())
  {
    Node_to_ArrayView_same_internal_uint64<type_selected(Types, conduit::DataType::UINT64_ID)>(func, first, args...);
  }
  else if(first.dtype().is_float32())
  {
    Node_to_ArrayView_same_internal_float32<type_selected(Types, conduit::DataType::FLOAT32_ID)>(func, first, args...);
  }
  else if(first.dtype().is_float64())
  {
    Node_to_ArrayView_same_internal_float64<type_selected(Types, conduit::DataType::FLOAT64_ID)>(func, first, args...);
  }
}

template <int Types = select_all_types(), typename FuncType, typename... Args>
void Node_to_ArrayView_same_internal(FuncType &&func, Delimiter, conduit::Node &first, Args&&... args)
{
  if(first.dtype().is_int8())
  {
    Node_to_ArrayView_same_internal_int8<type_selected(Types, conduit::DataType::INT8_ID)>(func, first, args...);
  }
  else if(first.dtype().is_int16())
  {
    Node_to_ArrayView_same_internal_int16<type_selected(Types, conduit::DataType::INT16_ID)>(func, first, args...);
  }
  else if(first.dtype().is_int32())
  {
    Node_to_ArrayView_same_internal_int32<type_selected(Types, conduit::DataType::INT32_ID)>(func, first, args...);
  }
  else if(first.dtype().is_int64())
  {
    Node_to_ArrayView_same_internal_int64<type_selected(Types, conduit::DataType::INT64_ID)>(func, first, args...);
  }
  else if(first.dtype().is_uint8())
  {
    Node_to_ArrayView_same_internal_uint8<type_selected(Types, conduit::DataType::UINT8_ID)>(func, first, args...);
  }
  else if(first.dtype().is_uint16())
  {
    Node_to_ArrayView_same_internal_uint16<type_selected(Types, conduit::DataType::UINT16_ID)>(func, first, args...);
  }
  else if(first.dtype().is_uint32())
  {
    Node_to_ArrayView_same_internal_uint32<type_selected(Types, conduit::DataType::UINT32_ID)>(func, first, args...);
  }
  else if(first.dtype().is_uint64())
  {
    Node_to_ArrayView_same_internal_uint64<type_selected(Types, conduit::DataType::UINT64_ID)>(func, first, args...);
  }
  else if(first.dtype().is_float32())
  {
    Node_to_ArrayView_same_internal_float32<type_selected(Types, conduit::DataType::FLOAT32_ID)>(func, first, args...);
  }
  else if(first.dtype().is_float64())
  {
    Node_to_ArrayView_same_internal_float64<type_selected(Types, conduit::DataType::FLOAT64_ID)>(func, first, args...);
  }
}

/// Reorder args
template <int Types = select_all_types(), typename ... Args>
void Node_to_ArrayView_same_internal(const conduit::Node &first, Args&&... args)
{
  Node_to_ArrayView_same_internal<Types>(args..., first);
}

template <int Types = select_all_types(), typename ... Args>
void Node_to_ArrayView_same_internal(conduit::Node &first, Args&&... args)
{
  Node_to_ArrayView_same_internal<Types>(args..., first);
}

} // namespace detail

//------------------------------------------------------------------------------
// Node to ArrayView. Handle all types.
//------------------------------------------------------------------------------

/*!
 * \brief Convert a series of Conduit nodes to axom::ArrayView and pass the concrete
 *        views to a lambda function passed as the last argument.
 *
 * \note  This method handles all Conduit array types and will instantiate views
 *        of any type. In other words, mixed node types can be used.
 *
 * \param first A Conduit node to be convered to a view.
 * \param args  A sequence of Conduit nodes, followed by a lambda that can accept
 *              the same number of views of any type.
 *
 * Node_to_ArrayView(node1, node2, [](auto &view1, auto &view2) { });
 * 
 */
template <typename ... Args>
void Node_to_ArrayView(const conduit::Node &first, Args&&... args)
{
  detail::Node_to_ArrayView_internal(first, args..., detail::ArgumentDelimiter);
}

template <typename ... Args>
void Node_to_ArrayView(conduit::Node &first, Args&&... args)
{
  detail::Node_to_ArrayView_internal(first, args..., detail::ArgumentDelimiter);
}

/*!
 * \brief Convert a series of Conduit nodes to axom::ArrayView and pass the concrete
 *        views to a lambda function passed as the last argument.
 *
 * \note  This method handles all Conduit array types. All nodes will be treated
 *        as the same type as the first Conduit node. Use this when all nodes
 *        are assumed to contain the same type since it results in less code.
 *
 * \param first A Conduit node to be convered to a view.
 * \param args  A sequence of Conduit nodes, followed by a lambda that can accept
 *              the same number of views of any type.
 *
 * Node_to_ArrayView_same(node1, node2, [](auto &view1, auto &view2) { });
 * 
 */
template <typename... Args>
void Node_to_ArrayView_same(const conduit::Node &first, Args&&... args)
{
  detail::Node_to_ArrayView_same_internal(first, args..., detail::ArgumentDelimiter);
}

template <typename... Args>
void Node_to_ArrayView_same(conduit::Node &first, Args&&... args)
{
  detail::Node_to_ArrayView_same_internal(first, args..., detail::ArgumentDelimiter);
}

//------------------------------------------------------------------------------
// Index Node to ArrayView. Handle types used for indexing.
//------------------------------------------------------------------------------

template <typename ... Args>
void IndexNode_to_ArrayView(const conduit::Node &first, Args&&... args)
{
  detail::Node_to_ArrayView_internal<detail::select_index_types()>(first, args..., detail::ArgumentDelimiter);
}

template <typename ... Args>
void IndexNode_to_ArrayView(conduit::Node &first, Args&&... args)
{
  detail::Node_to_ArrayView_internal<detail::select_index_types()>(first, args..., detail::ArgumentDelimiter);
}

template <typename... Args>
void IndexNode_to_ArrayView_same(const conduit::Node &first, Args&&... args)
{
  detail::Node_to_ArrayView_same_internal<detail::select_index_types()>(first, args..., detail::ArgumentDelimiter);
}

template <typename... Args>
void IndexNode_to_ArrayView_same(conduit::Node &first, Args&&... args)
{
  detail::Node_to_ArrayView_same_internal<detail::select_index_types()>(first, args..., detail::ArgumentDelimiter);
}

//------------------------------------------------------------------------------
// Float Node to ArrayView. Handle float types.
//------------------------------------------------------------------------------
template <typename ... Args>
void FloatNode_to_ArrayView(const conduit::Node &first, Args&&... args)
{
  detail::Node_to_ArrayView_internal<detail::select_float_types()>(first, args..., detail::ArgumentDelimiter);
}

template <typename ... Args>
void FloatNode_to_ArrayView(conduit::Node &first, Args&&... args)
{
  detail::Node_to_ArrayView_internal<detail::select_float_types()>(first, args..., detail::ArgumentDelimiter);
}

template <typename... Args>
void FloatNode_to_ArrayView_same(const conduit::Node &first, Args&&... args)
{
  detail::Node_to_ArrayView_same_internal<detail::select_float_types()>(first, args..., detail::ArgumentDelimiter);
}

template <typename... Args>
void FloatNode_to_ArrayView_same(conduit::Node &first, Args&&... args)
{
  detail::Node_to_ArrayView_same_internal<detail::select_float_types()>(first, args..., detail::ArgumentDelimiter);
}

} // namespace views
} // namespace mir
} // namespace axom

#endif
