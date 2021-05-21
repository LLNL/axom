// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_FOR_ALL_NODES_HPP_
#define MINT_FOR_ALL_NODES_HPP_

// Axom core includes
#include "axom/config.hpp"                          // compile time definitions
#include "axom/core/execution/execution_space.hpp"  // for execution_space traits
#include "axom/core/execution/for_all.hpp"          // for axom::for_all

// mint includes
#include "axom/mint/execution/xargs.hpp"       // for xargs
#include "axom/mint/config.hpp"                // for compile-time definitions
#include "axom/mint/mesh/Mesh.hpp"             // for mint::Mesh
#include "axom/mint/mesh/RectilinearMesh.hpp"  // for mint::RectilinearMesh
#include "axom/mint/mesh/StructuredMesh.hpp"   // for mint::StructuredMesh
#include "axom/mint/mesh/UniformMesh.hpp"      // for mint::UniformMesh
#include "axom/mint/execution/internal/structured_exec.hpp"

#include "axom/core/StackArray.hpp"  // for axom::StackArray

#ifdef AXOM_USE_RAJA
  #include "RAJA/RAJA.hpp"
#endif

namespace axom
{
namespace mint
{
namespace internal
{
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes_impl(xargs::index, const Mesh& m, KernelType&& kernel)
{
  const IndexType numNodes = m.getNumberOfNodes();
  axom::for_all<ExecPolicy>(numNodes, std::forward<KernelType>(kernel));
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes(xargs::index, const Mesh& m, KernelType&& kernel)
{
  for_all_nodes_impl<ExecPolicy>(xargs::index(),
                                 m,
                                 std::forward<KernelType>(kernel));
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes_impl(xargs::ij,
                               const StructuredMesh& m,
                               KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 2,
                "xargs::ij is only valid for 2D structured meshes!");

  const IndexType jp = m.nodeJp();
  const IndexType Ni = m.getNodeResolution(I_DIRECTION);
  const IndexType Nj = m.getNodeResolution(J_DIRECTION);

#ifdef AXOM_USE_RAJA

  RAJA::RangeSegment i_range(0, Ni);
  RAJA::RangeSegment j_range(0, Nj);
  using exec_pol = typename structured_exec<ExecPolicy>::loop2d_policy;

  RAJA::kernel<exec_pol>(
    RAJA::make_tuple(i_range, j_range),
    AXOM_LAMBDA(IndexType i, IndexType j) {
      const IndexType nodeIdx = i + j * jp;
      kernel(nodeIdx, i, j);
    });

#else

  constexpr bool is_serial = std::is_same<ExecPolicy, axom::SEQ_EXEC>::value;
  AXOM_STATIC_ASSERT(is_serial);

  for(IndexType j = 0; j < Nj; ++j)
  {
    const IndexType j_offset = j * jp;
    for(IndexType i = 0; i < Ni; ++i)
    {
      const IndexType nodeIdx = i + j_offset;
      kernel(nodeIdx, i, j);
    }  // END for all i
  }    // END for all j
#endif
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes(xargs::ij, const Mesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(!m.isStructured(),
                "xargs::ijk is only valid for 2D structured meshes!");
  SLIC_ERROR_IF(m.getDimension() != 2,
                "xargs::ijk is only valid for 2D structured meshes!");

  const StructuredMesh& sm = static_cast<const StructuredMesh&>(m);
  for_all_nodes_impl<ExecPolicy>(xargs::ij(),
                                 sm,
                                 std::forward<KernelType>(kernel));
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes_impl(xargs::ijk,
                               const StructuredMesh& m,
                               KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 3,
                "xargs::ijk is only valid for 3D structured meshes!");

  const IndexType jp = m.nodeJp();
  const IndexType kp = m.nodeKp();
  const IndexType Ni = m.getNodeResolution(I_DIRECTION);
  const IndexType Nj = m.getNodeResolution(J_DIRECTION);
  const IndexType Nk = m.getNodeResolution(K_DIRECTION);

#ifdef AXOM_USE_RAJA

  RAJA::RangeSegment i_range(0, Ni);
  RAJA::RangeSegment j_range(0, Nj);
  RAJA::RangeSegment k_range(0, Nk);
  using exec_pol = typename structured_exec<ExecPolicy>::loop3d_policy;

  RAJA::kernel<exec_pol>(
    RAJA::make_tuple(i_range, j_range, k_range),
    AXOM_LAMBDA(IndexType i, IndexType j, IndexType k) {
      const IndexType nodeIdx = i + j * jp + k * kp;
      kernel(nodeIdx, i, j, k);
    });

#else

  constexpr bool is_serial = std::is_same<ExecPolicy, axom::SEQ_EXEC>::value;
  AXOM_STATIC_ASSERT(is_serial);

  for(IndexType k = 0; k < Nk; ++k)
  {
    const IndexType k_offset = k * kp;
    for(IndexType j = 0; j < Nj; ++j)
    {
      const IndexType j_offset = j * jp;
      for(IndexType i = 0; i < Ni; ++i)
      {
        const IndexType nodeIdx = i + j_offset + k_offset;
        kernel(nodeIdx, i, j, k);
      }  // END for all i
    }    // END for all j
  }      // END for all k

#endif
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes(xargs::ijk, const Mesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(!m.isStructured(),
                "xargs::ijk is only valid for 3D structured meshes!");
  SLIC_ERROR_IF(m.getDimension() != 3,
                "xargs::ijk is only valid for 3D structured meshes!");

  const StructuredMesh& sm = static_cast<const StructuredMesh&>(m);
  for_all_nodes_impl<ExecPolicy>(xargs::ijk(),
                                 sm,
                                 std::forward<KernelType>(kernel));
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes_impl(xargs::x, const UniformMesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 1, "xargs::x is only valid for 1D meshes");

  const double x0 = m.getOrigin()[0];
  const double dx = m.getSpacing()[0];

  for_all_nodes_impl<ExecPolicy>(
    xargs::index(),
    m,
    AXOM_LAMBDA(IndexType nodeID) {
      const double x = x0 + nodeID * dx;
      kernel(nodeID, x);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes_impl(xargs::x, const Mesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 1, "xargs::x is only valid for 1D meshes");
  SLIC_ERROR_IF(m.getMeshType() == STRUCTURED_UNIFORM_MESH,
                "Not valid for UniformMesh.");

  const double* x = m.getCoordinateArray(X_COORDINATE);
  SLIC_ASSERT(x != nullptr);

  for_all_nodes_impl<ExecPolicy>(
    xargs::index(),
    m,
    AXOM_LAMBDA(IndexType nodeID) { kernel(nodeID, x[nodeID]); });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes(xargs::x, const Mesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 1, "xargs::x is only valid for 1D meshes");

  const int mesh_type = m.getMeshType();
  if(mesh_type == STRUCTURED_UNIFORM_MESH)
  {
    const UniformMesh& um = static_cast<const UniformMesh&>(m);
    for_all_nodes_impl<ExecPolicy>(xargs::x(),
                                   um,
                                   std::forward<KernelType>(kernel));
  }
  else
  {
    for_all_nodes_impl<ExecPolicy>(xargs::x(),
                                   m,
                                   std::forward<KernelType>(kernel));
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes_impl(xargs::xy, const UniformMesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 2, "xargs::xy is only valid for 2D meshes");

  const StackArray<double, 3>& origin = m.getOrigin();
  const StackArray<double, 3>& spacing = m.getSpacing();

  for_all_nodes_impl<ExecPolicy>(
    xargs::ij(),
    m,
    AXOM_LAMBDA(IndexType nodeID, IndexType i, IndexType j) {
      const double x = origin[0] + i * spacing[0];
      const double y = origin[1] + j * spacing[1];
      kernel(nodeID, x, y);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes_impl(xargs::xy,
                               const RectilinearMesh& m,
                               KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 2, "xargs::xy is only valid for 2D meshes");

  const double* x = m.getCoordinateArray(X_COORDINATE);
  const double* y = m.getCoordinateArray(Y_COORDINATE);
  SLIC_ASSERT(x != nullptr);
  SLIC_ASSERT(y != nullptr);

  for_all_nodes_impl<ExecPolicy>(
    xargs::ij(),
    m,
    AXOM_LAMBDA(IndexType nodeID, IndexType i, IndexType j) {
      kernel(nodeID, x[i], y[j]);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes_impl(xargs::xy, const Mesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 2, "xargs::xy is only valid for 2D meshes");
  SLIC_ERROR_IF(m.getMeshType() == STRUCTURED_UNIFORM_MESH,
                "Not valid for UniformMesh.");
  SLIC_ERROR_IF(m.getMeshType() == STRUCTURED_RECTILINEAR_MESH,
                "Not valid for RectilinearMesh.");

  const double* x = m.getCoordinateArray(X_COORDINATE);
  const double* y = m.getCoordinateArray(Y_COORDINATE);
  SLIC_ASSERT(x != nullptr);
  SLIC_ASSERT(y != nullptr);

  for_all_nodes_impl<ExecPolicy>(
    xargs::index(),
    m,
    AXOM_LAMBDA(IndexType nodeID) { kernel(nodeID, x[nodeID], y[nodeID]); });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes(xargs::xy, const Mesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 2, "xargs::xy is only valid for 3D meshes");

  const int mesh_type = m.getMeshType();
  if(mesh_type == STRUCTURED_RECTILINEAR_MESH)
  {
    const RectilinearMesh& rm = static_cast<const RectilinearMesh&>(m);
    for_all_nodes_impl<ExecPolicy>(xargs::xy(),
                                   rm,
                                   std::forward<KernelType>(kernel));
  }
  else if(mesh_type == STRUCTURED_UNIFORM_MESH)
  {
    const UniformMesh& um = static_cast<const UniformMesh&>(m);
    for_all_nodes_impl<ExecPolicy>(xargs::xy(),
                                   um,
                                   std::forward<KernelType>(kernel));
  }
  else
  {
    for_all_nodes_impl<ExecPolicy>(xargs::xy(),
                                   m,
                                   std::forward<KernelType>(kernel));
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes_impl(xargs::xyz, const UniformMesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 3, "xargs::xyz is only valid for 3D meshes");

  const StackArray<double, 3> origin = m.getOrigin();
  const StackArray<double, 3> spacing = m.getSpacing();

  for_all_nodes_impl<ExecPolicy>(
    xargs::ijk(),
    m,
    AXOM_LAMBDA(IndexType nodeID, IndexType i, IndexType j, IndexType k) {
      const double x = origin[0] + i * spacing[0];
      const double y = origin[1] + j * spacing[1];
      const double z = origin[2] + k * spacing[2];
      kernel(nodeID, x, y, z);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes_impl(xargs::xyz,
                               const RectilinearMesh& m,
                               KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 3, "xargs::xyz is only valid for 3D meshes");
  const double* x = m.getCoordinateArray(X_COORDINATE);
  const double* y = m.getCoordinateArray(Y_COORDINATE);
  const double* z = m.getCoordinateArray(Z_COORDINATE);

  for_all_nodes_impl<ExecPolicy>(
    xargs::ijk(),
    m,
    AXOM_LAMBDA(IndexType nodeID, IndexType i, IndexType j, IndexType k) {
      kernel(nodeID, x[i], y[j], z[k]);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes_impl(xargs::xyz, const Mesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 3, "xargs::xyz is only valid for 3D meshes");
  SLIC_ERROR_IF(m.getMeshType() == STRUCTURED_UNIFORM_MESH,
                "Not valid for UniformMesh.");
  SLIC_ERROR_IF(m.getMeshType() == STRUCTURED_RECTILINEAR_MESH,
                "Not valid for RectilinearMesh.");

  const double* x = m.getCoordinateArray(X_COORDINATE);
  const double* y = m.getCoordinateArray(Y_COORDINATE);
  const double* z = m.getCoordinateArray(Z_COORDINATE);
  SLIC_ASSERT(x != nullptr);
  SLIC_ASSERT(y != nullptr);
  SLIC_ASSERT(z != nullptr);

  for_all_nodes_impl<ExecPolicy>(
    xargs::index(),
    m,
    AXOM_LAMBDA(IndexType nodeID) {
      kernel(nodeID, x[nodeID], y[nodeID], z[nodeID]);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes(xargs::xyz, const mint::Mesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 3, "xargs::xyz is only valid for 3D meshes");

  const int mesh_type = m.getMeshType();
  if(mesh_type == STRUCTURED_RECTILINEAR_MESH)
  {
    const RectilinearMesh& rm = static_cast<const RectilinearMesh&>(m);
    for_all_nodes_impl<ExecPolicy>(xargs::xyz(),
                                   rm,
                                   std::forward<KernelType>(kernel));
  }
  else if(mesh_type == STRUCTURED_UNIFORM_MESH)
  {
    const UniformMesh& um = static_cast<const UniformMesh&>(m);
    for_all_nodes_impl<ExecPolicy>(xargs::xyz(),
                                   um,
                                   std::forward<KernelType>(kernel));
  }
  else
  {
    for_all_nodes_impl<ExecPolicy>(xargs::xyz(),
                                   m,
                                   std::forward<KernelType>(kernel));
  }
}

} /* namespace internal */
} /* namespace mint     */
} /* namespace axom     */

#endif /* MINTFOR_ALL_NODES_STRUCTURED_HPP_ */
