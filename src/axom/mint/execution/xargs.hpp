/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifndef MINT_EXECUTION_ARGS_HPP_
#define MINT_EXECUTION_ARGS_HPP_

/*!
 * \file
 *
 * \brief The xargs::* types are used as template arguments to the execution
 *  model loop/mesh traversal functions. They are used to indicate additional
 *  arguments to the lambda expression that is passed to the loop/mesh traversal
 *  function.
 *
 *  \see interface.hpp
 */

namespace axom
{
namespace mint
{

namespace xargs
{
  /*!
   * \brief Indicates that the the lambda expression takes the loop index.
   *
   * \note This is the default argument passed to the lambda expression and
   *  does not need to be specified explicitly.
   */
  struct index{ };

  /*!
   * \brief Indicates that the lambda expression takes the corresponding
   *  IJ lattice grid coordinates in addition to the associated linear index.
   *
   * \note This option can be used for node and/or cell loop traversals on
   *  Structured meshes.
   */
  struct ij{ };

  /*!
   * \brief Indicates that the lambda expression takes the corresponding IJK
   *  lattice grid coordinates in addition to the associated linear index.
   *
   * \note This option can be used for node and/or cell loop traversals on
   *  Structured meshes.
   */
  struct ijk{ };

  /*!
   * \brief Indicates that the lambda expression takes the corresponding X
   *  coordinate of a mesh node in addition to the associated node index.
   *
   * \note This option can be be use for node mesh traversals with any 1D mesh.
   */
  struct x{};

  /*!
   * \brief Indicates that the lambda expression takes the corresponding X,Y,
   *  coordinates of a mesh node in addition to the associated node index.
   *
   * \note This option can be be used for node mesh traversals with any 2D mesh.
   */
  struct xy{ };

  /*!
   * \brief Indicates that the lambda expression takes the corresponding X,Y,Z
   * coordinates of a mesh node in addition to the associated node index.
   *
   * \note This option can be used for node mesh traversals with any 3D mesh.
   */
  struct xyz{ };

  /*!
   * \brief Indicates that the lambda expression also takes the node
   *  connectivity information, i.e., the node IDs of the mesh entity being
   *  traversed, e.g., cells or faces and the number of nodes.
   *
   * \note This option can be use for cell/face mesh traversals with any mesh.
   */
  struct nodeids{ };

} /* namespace xargs */

/*!
 * \brief Traits class used for compile-time checking of xargs types.
 */
template < typename ArgType >
struct xargs_traits
{
  static constexpr bool valid() { return false; };
  static constexpr char* name() { return (char*)""; };
};

template < >
struct xargs_traits< xargs::index >
{
  static constexpr bool valid() { return true; };
  static constexpr char* name() { return (char*)("xargs::index"); };
};

template < >
struct xargs_traits< xargs::ij >
{
  static constexpr bool valid() { return true; };
  static constexpr char* name() { return (char*)("xargs::ij"); };
};

template < >
struct xargs_traits< xargs::ijk >
{
  static constexpr bool valid() { return true; };
  static constexpr char* name() { return (char*)("xargs::ijk"); };
};

template < >
struct xargs_traits< xargs::x >
{
  static constexpr bool valid() { return true; };
  static constexpr char* name() { return (char*)("xargs::x"); };
};

template < >
struct xargs_traits< xargs::xy >
{
  static constexpr bool valid() { return true; };
  static constexpr char* name() { return (char*)("xargs::xy"); };
};

template < >
struct xargs_traits< xargs::xyz >
{
  static constexpr bool valid() { return true; };
  static constexpr char* name() { return (char*)("xargs::xyz"); };
};

template < >
struct xargs_traits< xargs::nodeids >
{
  static constexpr bool valid() { return true; };
  static constexpr char* name() { return (char*)("xargs::nodeids"); };
};


} /* namespace mint */
} /* namespace axom */

#endif /* MINT_EXEC_ARGS_HPP_ */
