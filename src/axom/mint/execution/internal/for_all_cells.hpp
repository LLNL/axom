// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_FOR_ALL_CELLS_HPP_
#define MINT_FOR_ALL_CELLS_HPP_

// Axom core includes
#include "axom/config.hpp"                          // compile time definitions
#include "axom/core/execution/execution_space.hpp"  // for execution_space traits
#include "axom/core/execution/for_all.hpp"          // for axom::for_all

// mint includes
#include "axom/mint/execution/xargs.hpp"  // for xargs

#include "axom/mint/config.hpp"                 // for compile-time definitions
#include "axom/mint/mesh/Mesh.hpp"              // for Mesh
#include "axom/mint/mesh/StructuredMesh.hpp"    // for StructuredMesh
#include "axom/mint/mesh/UniformMesh.hpp"       // for UniformMesh
#include "axom/mint/mesh/RectilinearMesh.hpp"   // for RectilinearMesh
#include "axom/mint/mesh/CurvilinearMesh.hpp"   // for CurvilinearMesh
#include "axom/mint/mesh/UnstructuredMesh.hpp"  // for UnstructuredMesh
#include "axom/mint/execution/internal/helpers.hpp"
#include "axom/mint/execution/internal/structured_exec.hpp"

#include "axom/core/StackArray.hpp"       // for axom::StackArray
#include "axom/core/numerics/Matrix.hpp"  // for Matrix

#ifdef AXOM_USE_RAJA
  #include "RAJA/RAJA.hpp"
#endif

namespace axom
{
namespace mint
{
namespace internal
{
//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells_impl(xargs::index, const Mesh& m, KernelType&& kernel)
{
  const IndexType numCells = m.getNumberOfCells();
  axom::for_all<ExecPolicy>(numCells, std::forward<KernelType>(kernel));
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells(xargs::index, const Mesh& m, KernelType&& kernel)
{
  for_all_cells_impl<ExecPolicy>(xargs::index(),
                                 m,
                                 std::forward<KernelType>(kernel));
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells_impl(xargs::ij,
                               const StructuredMesh& m,
                               KernelType&& kernel)
{
  // run-time checks
  SLIC_ERROR_IF(m.getDimension() != 2, "xargs::ij is only valid for 2D meshes!");

  const IndexType jp = m.cellJp();
  const IndexType Ni = m.getCellResolution(I_DIRECTION);
  const IndexType Nj = m.getCellResolution(J_DIRECTION);

#ifdef AXOM_USE_RAJA

  RAJA::RangeSegment i_range(0, Ni);
  RAJA::RangeSegment j_range(0, Nj);
  using exec_pol = typename structured_exec<ExecPolicy>::loop2d_policy;

  RAJA::kernel<exec_pol>(
    RAJA::make_tuple(i_range, j_range),
    AXOM_LAMBDA(IndexType i, IndexType j) {
      const IndexType cellID = i + j * jp;
      kernel(cellID, i, j);
    });

#else

  constexpr bool is_serial = std::is_same<ExecPolicy, axom::SEQ_EXEC>::value;
  AXOM_STATIC_ASSERT(is_serial);

  for(IndexType j = 0; j < Nj; ++j)
  {
    const IndexType j_offset = j * jp;
    for(IndexType i = 0; i < Ni; ++i)
    {
      const IndexType cellID = i + j_offset;
      kernel(cellID, i, j);
    }  // END for all i
  }    // END for all j

#endif
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells(xargs::ij, const Mesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(!m.isStructured(),
                "xargs::ij is only valid on structured meshes!");
  SLIC_ERROR_IF(m.getDimension() != 2, "xargs::ij is only valid on 2D meshes!");

  const StructuredMesh& sm = static_cast<const StructuredMesh&>(m);
  for_all_cells_impl<ExecPolicy>(xargs::ij(),
                                 sm,
                                 std::forward<KernelType>(kernel));
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells_impl(xargs::ijk,
                               const StructuredMesh& m,
                               KernelType&& kernel)
{
  const IndexType Ni = m.getCellResolution(I_DIRECTION);
  const IndexType Nj = m.getCellResolution(J_DIRECTION);
  const IndexType Nk = m.getCellResolution(K_DIRECTION);

  const IndexType jp = m.cellJp();
  const IndexType kp = m.cellKp();

#ifdef AXOM_USE_RAJA

  RAJA::RangeSegment i_range(0, Ni);
  RAJA::RangeSegment j_range(0, Nj);
  RAJA::RangeSegment k_range(0, Nk);
  using exec_pol = typename structured_exec<ExecPolicy>::loop3d_policy;

  RAJA::kernel<exec_pol>(
    RAJA::make_tuple(i_range, j_range, k_range),
    AXOM_LAMBDA(IndexType i, IndexType j, IndexType k) {
      const IndexType cellID = i + j * jp + k * kp;
      kernel(cellID, i, j, k);
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
        const IndexType cellID = i + j_offset + k_offset;
        kernel(cellID, i, j, k);
      }  // END for all i
    }    // END for all j
  }      // END for all k

#endif
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells(xargs::ijk, const Mesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(!m.isStructured(),
                "xargs::ijk is only valid on structured meshes!");
  SLIC_ERROR_IF(m.getDimension() != 3,
                "xargs::ijk is only valid for 3D structured meshes!");

  const StructuredMesh& sm = static_cast<const StructuredMesh&>(m);
  for_all_cells_impl<ExecPolicy>(xargs::ijk(),
                                 sm,
                                 std::forward<KernelType>(kernel));
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells_impl(xargs::nodeids,
                               const StructuredMesh& m,
                               KernelType&& kernel)
{
  const int dimension = m.getDimension();
  const IndexType nodeJp = m.nodeJp();
  const IndexType nodeKp = m.nodeKp();
  const StackArray<IndexType, 8>& offsets = m.getCellNodeOffsetsArray();

  if(dimension == 1)
  {
    for_all_cells_impl<ExecPolicy>(
      xargs::index(),
      m,
      AXOM_LAMBDA(IndexType cellID) {
        IndexType cell_connectivity[2] = {cellID, cellID + 1};
        kernel(cellID, cell_connectivity, 2);
      });
  }
  else if(dimension == 2)
  {
    for_all_cells_impl<ExecPolicy>(
      xargs::ij(),
      m,
      AXOM_LAMBDA(IndexType cellID, IndexType i, IndexType j) {
        const IndexType n0 = i + j * nodeJp;
        IndexType cell_connectivity[4];

        for(int ii = 0; ii < 4; ++ii)
        {
          cell_connectivity[ii] = n0 + offsets[ii];
        }

        kernel(cellID, cell_connectivity, 4);
      });
  }
  else
  {
    SLIC_ASSERT(dimension == 3);

    for_all_cells_impl<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType cellID, IndexType i, IndexType j, IndexType k) {
        const IndexType n0 = i + j * nodeJp + k * nodeKp;
        IndexType cell_connectivity[8];

        for(int ii = 0; ii < 8; ++ii)
        {
          cell_connectivity[ii] = n0 + offsets[ii];
        }

        kernel(cellID, cell_connectivity, 8);
      });
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells_impl(xargs::nodeids,
                               const UnstructuredMesh<MIXED_SHAPE>& m,
                               KernelType&& kernel)
{
  const IndexType* cell_connectivity = m.getCellNodesArray();
  const IndexType* cell_offsets = m.getCellNodesOffsetsArray();

  for_all_cells_impl<ExecPolicy>(
    xargs::index(),
    m,
    AXOM_LAMBDA(IndexType cellID) {
      const IndexType N = cell_offsets[cellID + 1] - cell_offsets[cellID];
      kernel(cellID, &cell_connectivity[cell_offsets[cellID]], N);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells_impl(xargs::nodeids,
                               const UnstructuredMesh<SINGLE_SHAPE>& m,
                               KernelType&& kernel)
{
  const IndexType* cell_connectivity = m.getCellNodesArray();
  const IndexType stride = m.getNumberOfCellNodes();

  for_all_cells_impl<ExecPolicy>(
    xargs::index(),
    m,
    AXOM_LAMBDA(IndexType cellID) {
      kernel(cellID, &cell_connectivity[cellID * stride], stride);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells(xargs::nodeids, const Mesh& m, KernelType&& kernel)
{
  if(m.isStructured())
  {
    const StructuredMesh& sm = static_cast<const StructuredMesh&>(m);
    for_all_cells_impl<ExecPolicy>(xargs::nodeids(),
                                   sm,
                                   std::forward<KernelType>(kernel));
  }
  else if(m.hasMixedCellTypes())
  {
    const UnstructuredMesh<MIXED_SHAPE>& um =
      static_cast<const UnstructuredMesh<MIXED_SHAPE>&>(m);
    for_all_cells_impl<ExecPolicy>(xargs::nodeids(),
                                   um,
                                   std::forward<KernelType>(kernel));
  }
  else
  {
    const UnstructuredMesh<SINGLE_SHAPE>& um =
      static_cast<const UnstructuredMesh<SINGLE_SHAPE>&>(m);
    for_all_cells_impl<ExecPolicy>(xargs::nodeids(),
                                   um,
                                   std::forward<KernelType>(kernel));
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells_impl(xargs::faceids,
                               const StructuredMesh& m,
                               KernelType&& kernel)
{
  const int dimension = m.getDimension();
  const IndexType ICellResolution = m.getCellResolution(I_DIRECTION);
  const IndexType numIFaces = m.getTotalNumFaces(I_DIRECTION);

  if(dimension == 2)
  {
    for_all_cells_impl<ExecPolicy>(
      xargs::ij(),
      m,
      AXOM_LAMBDA(IndexType cellID, IndexType AXOM_NOT_USED(i), IndexType j) {
        IndexType faces[4];

        /* The I_DIRECTION faces */
        faces[0] = cellID + j;
        faces[1] = faces[0] + 1;

        /* The J_DIRECTION faces */
        faces[2] = cellID + numIFaces;
        faces[3] = faces[2] + ICellResolution;

        kernel(cellID, faces, 4);
      });
  }
  else
  {
    SLIC_ASSERT(dimension == 3);

    const IndexType JCellResolution = m.getCellResolution(J_DIRECTION);
    const IndexType totalIJfaces = numIFaces + m.getTotalNumFaces(J_DIRECTION);
    const IndexType cellKp = m.cellKp();

    for_all_cells_impl<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType cellID,
                  IndexType AXOM_NOT_USED(i),
                  IndexType j,
                  IndexType k) {
        IndexType faces[6];

        /* The I_DIRECTION faces */
        faces[0] = cellID + j + JCellResolution * k;
        faces[1] = faces[0] + 1;

        /* The J_DIRECTION faces */
        faces[2] = cellID + numIFaces + ICellResolution * k;
        faces[3] = faces[2] + ICellResolution;

        /* The K_DIRECTION faces */
        faces[4] = cellID + totalIJfaces;
        faces[5] = faces[4] + cellKp;

        kernel(cellID, faces, 6);
      });
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells_impl(xargs::faceids,
                               const UnstructuredMesh<SINGLE_SHAPE>& m,
                               KernelType&& kernel)
{
  const IndexType* cells_to_faces = m.getCellFacesArray();
  const IndexType num_faces = m.getNumberOfCellFaces();

  for_all_cells_impl<ExecPolicy>(
    xargs::index(),
    m,
    AXOM_LAMBDA(IndexType cellID) {
      kernel(cellID, cells_to_faces + cellID * num_faces, num_faces);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells_impl(xargs::faceids,
                               const UnstructuredMesh<MIXED_SHAPE>& m,
                               KernelType&& kernel)
{
  const IndexType* cells_to_faces = m.getCellFacesArray();
  const IndexType* offsets = m.getCellFacesOffsetsArray();

  for_all_cells_impl<ExecPolicy>(
    xargs::index(),
    m,
    AXOM_LAMBDA(IndexType cellID) {
      const IndexType num_faces = offsets[cellID + 1] - offsets[cellID];
      kernel(cellID, cells_to_faces + offsets[cellID], num_faces);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells(xargs::faceids, const Mesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(
    m.getDimension() == 1,
    "For all cells with face IDs only supported for 2D and 3D meshes");

  if(m.isStructured())
  {
    const StructuredMesh& sm = static_cast<const StructuredMesh&>(m);
    for_all_cells_impl<ExecPolicy>(xargs::faceids(),
                                   sm,
                                   std::forward<KernelType>(kernel));
  }
  else if(m.hasMixedCellTypes())
  {
    const UnstructuredMesh<MIXED_SHAPE>& um =
      static_cast<const UnstructuredMesh<MIXED_SHAPE>&>(m);
    for_all_cells_impl<ExecPolicy>(xargs::faceids(),
                                   um,
                                   std::forward<KernelType>(kernel));
  }
  else
  {
    const UnstructuredMesh<SINGLE_SHAPE>& um =
      static_cast<const UnstructuredMesh<SINGLE_SHAPE>&>(m);
    for_all_cells_impl<ExecPolicy>(xargs::faceids(),
                                   um,
                                   std::forward<KernelType>(kernel));
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells_impl(xargs::coords,
                               const UniformMesh& m,
                               KernelType&& kernel)
{
  constexpr bool NO_COPY = true;

  const int dimension = m.getDimension();
  const double* origin = m.getOrigin();
  const double* spacing = m.getSpacing();
  const IndexType nodeJp = m.nodeJp();
  const IndexType nodeKp = m.nodeKp();

  const double x0 = origin[0];
  const double dx = spacing[0];

  const double y0 = origin[1];
  const double dy = spacing[1];

  const double z0 = origin[2];
  const double dz = spacing[2];

  if(dimension == 1)
  {
    for_all_cells_impl<ExecPolicy>(
      xargs::index(),
      m,
      AXOM_LAMBDA(IndexType cellID) {
        const IndexType nodeIDs[2] = {cellID, cellID + 1};
        double coords[2] = {x0 + nodeIDs[0] * dx, x0 + nodeIDs[1] * dx};

        numerics::Matrix<double> coordsMatrix(dimension, 2, coords, NO_COPY);
        kernel(cellID, coordsMatrix, nodeIDs);
      });
  }
  else if(dimension == 2)
  {
    for_all_cells_impl<ExecPolicy>(
      xargs::ij(),
      m,
      AXOM_LAMBDA(IndexType cellID, IndexType i, IndexType j) {
        const IndexType n0 = i + j * nodeJp;
        const IndexType nodeIDs[4] = {n0, n0 + 1, n0 + 1 + nodeJp, n0 + nodeJp};

        double coords[8] = {x0 + i * dx,
                            y0 + j * dy,
                            x0 + (i + 1) * dx,
                            y0 + j * dy,
                            x0 + (i + 1) * dx,
                            y0 + (j + 1) * dy,
                            x0 + i * dx,
                            y0 + (j + 1) * dy};

        numerics::Matrix<double> coordsMatrix(dimension, 4, coords, NO_COPY);
        kernel(cellID, coordsMatrix, nodeIDs);
      });
  }
  else
  {
    SLIC_ASSERT(dimension == 3);
    for_all_cells_impl<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType cellID, IndexType i, IndexType j, IndexType k) {
        const IndexType n0 = i + j * nodeJp + k * nodeKp;
        const IndexType nodeIDs[8] = {n0,
                                      n0 + 1,
                                      n0 + 1 + nodeJp,
                                      n0 + nodeJp,
                                      n0 + nodeKp,
                                      n0 + 1 + nodeKp,
                                      n0 + 1 + nodeJp + nodeKp,
                                      n0 + nodeJp + nodeKp};

        double coords[24] = {
          x0 + i * dx,       y0 + j * dy,       z0 + k * dz,
          x0 + (i + 1) * dx, y0 + j * dy,       z0 + k * dz,
          x0 + (i + 1) * dx, y0 + (j + 1) * dy, z0 + k * dz,
          x0 + i * dx,       y0 + (j + 1) * dy, z0 + k * dz,
          x0 + i * dx,       y0 + j * dy,       z0 + (k + 1) * dz,
          x0 + (i + 1) * dx, y0 + j * dy,       z0 + (k + 1) * dz,
          x0 + (i + 1) * dx, y0 + (j + 1) * dy, z0 + (k + 1) * dz,
          x0 + i * dx,       y0 + (j + 1) * dy, z0 + (k + 1) * dz};

        numerics::Matrix<double> coordsMatrix(dimension, 8, coords, NO_COPY);
        kernel(cellID, coordsMatrix, nodeIDs);
      });
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells_impl(xargs::coords,
                               const RectilinearMesh& m,
                               KernelType&& kernel)
{
  constexpr bool NO_COPY = true;

  const int dimension = m.getDimension();
  const IndexType nodeJp = m.nodeJp();
  const IndexType nodeKp = m.nodeKp();
  const double* x = m.getCoordinateArray(X_COORDINATE);

  if(dimension == 1)
  {
    for_all_cells_impl<ExecPolicy>(
      xargs::index(),
      m,
      AXOM_LAMBDA(IndexType cellID) {
        const IndexType nodeIDs[2] = {cellID, cellID + 1};
        double coords[2] = {x[nodeIDs[0]], x[nodeIDs[1]]};

        numerics::Matrix<double> coordsMatrix(dimension, 2, coords, NO_COPY);
        kernel(cellID, coordsMatrix, nodeIDs);
      });
  }
  else if(dimension == 2)
  {
    const double* y = m.getCoordinateArray(Y_COORDINATE);
    for_all_cells_impl<ExecPolicy>(
      xargs::ij(),
      m,
      AXOM_LAMBDA(IndexType cellID, IndexType i, IndexType j) {
        const IndexType n0 = i + j * nodeJp;
        const IndexType nodeIDs[4] = {n0, n0 + 1, n0 + 1 + nodeJp, n0 + nodeJp};

        double coords[8] =
          {x[i], y[j], x[i + 1], y[j], x[i + 1], y[j + 1], x[i], y[j + 1]};

        numerics::Matrix<double> coordsMatrix(dimension, 4, coords, NO_COPY);
        kernel(cellID, coordsMatrix, nodeIDs);
      });
  }
  else
  {
    SLIC_ASSERT(dimension == 3);
    const double* y = m.getCoordinateArray(Y_COORDINATE);
    const double* z = m.getCoordinateArray(Z_COORDINATE);
    for_all_cells_impl<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType cellID, IndexType i, IndexType j, IndexType k) {
        const IndexType n0 = i + j * nodeJp + k * nodeKp;
        const IndexType nodeIDs[8] = {n0,
                                      n0 + 1,
                                      n0 + 1 + nodeJp,
                                      n0 + nodeJp,
                                      n0 + nodeKp,
                                      n0 + 1 + nodeKp,
                                      n0 + 1 + nodeJp + nodeKp,
                                      n0 + nodeJp + nodeKp};

        double coords[24] = {x[i],     y[j],     z[k],     x[i + 1], y[j],
                             z[k],     x[i + 1], y[j + 1], z[k],     x[i],
                             y[j + 1], z[k],     x[i],     y[j],     z[k + 1],
                             x[i + 1], y[j],     z[k + 1], x[i + 1], y[j + 1],
                             z[k + 1], x[i],     y[j + 1], z[k + 1]};

        numerics::Matrix<double> coordsMatrix(dimension, 8, coords, NO_COPY);
        kernel(cellID, coordsMatrix, nodeIDs);
      });
  }
}

//------------------------------------------------------------------------------
struct for_all_cell_nodes_functor
{
  template <typename ExecPolicy, typename MeshType, typename KernelType>
  inline void operator()(ExecPolicy AXOM_NOT_USED(policy),
                         const MeshType& m,
                         KernelType&& kernel) const
  {
    constexpr bool valid_mesh_type = std::is_base_of<Mesh, MeshType>::value;
    AXOM_STATIC_ASSERT(valid_mesh_type);

    for_all_cells_impl<ExecPolicy>(xargs::nodeids(),
                                   m,
                                   std::forward<KernelType>(kernel));
  }
};

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells_impl(xargs::coords,
                               const CurvilinearMesh& m,
                               KernelType&& kernel)
{
  const int dimension = m.getDimension();
  if(dimension == 1)
  {
    for_all_coords<ExecPolicy, 1, 2>(for_all_cell_nodes_functor(),
                                     m,
                                     std::forward<KernelType>(kernel));
  }
  else if(dimension == 2)
  {
    for_all_coords<ExecPolicy, 2, 4>(for_all_cell_nodes_functor(),
                                     m,
                                     std::forward<KernelType>(kernel));
  }
  else
  {
    for_all_coords<ExecPolicy, 3, 8>(for_all_cell_nodes_functor(),
                                     m,
                                     std::forward<KernelType>(kernel));
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType, Topology TOPO>
inline void for_all_cells_impl(xargs::coords,
                               const UnstructuredMesh<TOPO>& m,
                               KernelType&& kernel)
{
  constexpr bool NO_COPY = true;

  const int dimension = m.getDimension();
  const double* x = m.getCoordinateArray(X_COORDINATE);

  if(dimension == 1)
  {
    for_all_cells_impl<ExecPolicy>(
      xargs::nodeids(),
      m,
      AXOM_LAMBDA(IndexType cellID,
                  const IndexType* nodeIDs,
                  IndexType AXOM_NOT_USED(numNodes)) {
        double coords[2] = {x[nodeIDs[0]], x[nodeIDs[1]]};

        numerics::Matrix<double> coordsMatrix(dimension, 2, coords, NO_COPY);
        kernel(cellID, coordsMatrix, nodeIDs);
      });
  }
  else if(dimension == 2)
  {
    const double* y = m.getCoordinateArray(Y_COORDINATE);
    for_all_cells_impl<ExecPolicy>(
      xargs::nodeids(),
      m,
      AXOM_LAMBDA(IndexType cellID, const IndexType* nodeIDs, IndexType numNodes) {
        double coords[2 * MAX_CELL_NODES];
        for(int i = 0; i < numNodes; ++i)
        {
          const IndexType nodeID = nodeIDs[i];
          coords[2 * i] = x[nodeID];
          coords[2 * i + 1] = y[nodeID];
        }

        numerics::Matrix<double> coordsMatrix(dimension, numNodes, coords, NO_COPY);
        kernel(cellID, coordsMatrix, nodeIDs);
      });
  }
  else
  {
    SLIC_ASSERT(dimension == 3);
    const double* y = m.getCoordinateArray(Y_COORDINATE);
    const double* z = m.getCoordinateArray(Z_COORDINATE);
    for_all_cells_impl<ExecPolicy>(
      xargs::nodeids(),
      m,
      AXOM_LAMBDA(IndexType cellID, const IndexType* nodeIDs, IndexType numNodes) {
        double coords[3 * MAX_CELL_NODES];
        for(int i = 0; i < numNodes; ++i)
        {
          const IndexType nodeID = nodeIDs[i];
          coords[3 * i] = x[nodeID];
          coords[3 * i + 1] = y[nodeID];
          coords[3 * i + 2] = z[nodeID];
        }

        numerics::Matrix<double> coordsMatrix(dimension, numNodes, coords, NO_COPY);
        kernel(cellID, coordsMatrix, nodeIDs);
      });
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells(xargs::coords, const Mesh& m, KernelType&& kernel)
{
  if(m.getMeshType() == STRUCTURED_UNIFORM_MESH)
  {
    const UniformMesh& um = static_cast<const UniformMesh&>(m);
    for_all_cells_impl<ExecPolicy>(xargs::coords(),
                                   um,
                                   std::forward<KernelType>(kernel));
  }
  else if(m.getMeshType() == STRUCTURED_RECTILINEAR_MESH)
  {
    const RectilinearMesh& rm = static_cast<const RectilinearMesh&>(m);
    for_all_cells_impl<ExecPolicy>(xargs::coords(),
                                   rm,
                                   std::forward<KernelType>(kernel));
  }
  else if(m.getMeshType() == STRUCTURED_CURVILINEAR_MESH)
  {
    const CurvilinearMesh& cm = static_cast<const CurvilinearMesh&>(m);
    for_all_cells_impl<ExecPolicy>(xargs::coords(),
                                   cm,
                                   std::forward<KernelType>(kernel));
  }
  else if(m.getMeshType() == UNSTRUCTURED_MESH)
  {
    if(m.hasMixedCellTypes())
    {
      const UnstructuredMesh<MIXED_SHAPE>& um =
        static_cast<const UnstructuredMesh<MIXED_SHAPE>&>(m);

      for_all_cells_impl<ExecPolicy>(xargs::coords(),
                                     um,
                                     std::forward<KernelType>(kernel));
    }
    else
    {
      const UnstructuredMesh<SINGLE_SHAPE>& um =
        static_cast<const UnstructuredMesh<SINGLE_SHAPE>&>(m);

      for_all_cells_impl<ExecPolicy>(xargs::coords(),
                                     um,
                                     std::forward<KernelType>(kernel));
    }
  }
  else
  {
    SLIC_ERROR("Unsupported mesh type: " << m.getMeshType());
  }
}

} /* namespace internal */
} /* namespace mint     */
} /* namespace axom     */

#endif /* MINT_FOR_ALL_CELLS_HPP_ */
