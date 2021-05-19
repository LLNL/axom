// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_FOR_ALL_FACES_HPP_
#define MINT_FOR_ALL_FACES_HPP_

// Axom core includes
#include "axom/config.hpp"                          // compile time definitions
#include "axom/core/execution/execution_space.hpp"  // for execution_space traits
#include "axom/core/execution/for_all.hpp"          // for axom::for_all

// mint includes
#include "axom/mint/execution/xargs.hpp"       // for xargs
#include "axom/mint/config.hpp"                // for compile-time definitions
#include "axom/mint/mesh/Mesh.hpp"             // for Mesh
#include "axom/mint/mesh/StructuredMesh.hpp"   // for StructuredMesh
#include "axom/mint/mesh/UniformMesh.hpp"      // for UniformMesh
#include "axom/mint/mesh/RectilinearMesh.hpp"  // for RectilinearMesh
#include "axom/mint/mesh/CurvilinearMesh.hpp"  // for CurvilinearMesh
#include "axom/mint/execution/internal/helpers.hpp"  // for for_all_coords
#include "axom/mint/execution/internal/structured_exec.hpp"

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
namespace helpers
{
//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_I_faces(xargs::ij, const StructuredMesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 2, "Mesh must be 2D.");

  const IndexType INodeResolution = m.getNodeResolution(I_DIRECTION);
  const IndexType Ni = INodeResolution;
  const IndexType Nj = m.getCellResolution(J_DIRECTION);

#ifdef AXOM_USE_RAJA

  RAJA::RangeSegment i_range(0, Ni);
  RAJA::RangeSegment j_range(0, Nj);

  using exec_pol = typename structured_exec<ExecPolicy>::loop2d_policy;
  RAJA::kernel<exec_pol>(
    RAJA::make_tuple(i_range, j_range),
    AXOM_LAMBDA(IndexType i, IndexType j) {
      const IndexType faceID = i + j * INodeResolution;
      kernel(faceID, i, j);
    });

#else

  constexpr bool is_serial = std::is_same<ExecPolicy, axom::SEQ_EXEC>::value;
  AXOM_STATIC_ASSERT(is_serial);

  for(IndexType j = 0; j < Nj; ++j)
  {
    const IndexType offset = j * INodeResolution;
    for(IndexType i = 0; i < Ni; ++i)
    {
      const IndexType faceID = i + offset;
      kernel(faceID, i, j);
    }
  }

#endif
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_I_faces(xargs::ijk, const StructuredMesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 3, "Mesh must be a 3D.");

  const IndexType INodeResolution = m.getNodeResolution(I_DIRECTION);
  const IndexType numIFacesInKSlice =
    INodeResolution * m.getCellResolution(J_DIRECTION);
  const IndexType Ni = INodeResolution;
  const IndexType Nj = m.getCellResolution(J_DIRECTION);
  const IndexType Nk = m.getCellResolution(K_DIRECTION);

#ifdef AXOM_USE_RAJA

  RAJA::RangeSegment i_range(0, Ni);
  RAJA::RangeSegment j_range(0, Nj);
  RAJA::RangeSegment k_range(0, Nk);

  using exec_pol = typename structured_exec<ExecPolicy>::loop3d_policy;
  RAJA::kernel<exec_pol>(
    RAJA::make_tuple(i_range, j_range, k_range),
    AXOM_LAMBDA(IndexType i, IndexType j, IndexType k) {
      const IndexType faceID = i + j * INodeResolution + k * numIFacesInKSlice;
      kernel(faceID, i, j, k);
    });

#else

  constexpr bool is_serial = std::is_same<ExecPolicy, axom::SEQ_EXEC>::value;
  AXOM_STATIC_ASSERT(is_serial);

  for(IndexType k = 0; k < Nk; ++k)
  {
    const IndexType k_offset = k * numIFacesInKSlice;
    for(IndexType j = 0; j < Nj; ++j)
    {
      const IndexType offset = j * INodeResolution + k_offset;
      for(IndexType i = 0; i < Ni; ++i)
      {
        const IndexType faceID = i + offset;
        kernel(faceID, i, j, k);
      }
    }
  }

#endif
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_J_faces(xargs::ij, const StructuredMesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 2, "Mesh must be 2D.");

  const IndexType ICellResolution = m.getCellResolution(I_DIRECTION);
  const IndexType numIFaces = m.getTotalNumFaces(I_DIRECTION);
  const IndexType Ni = ICellResolution;
  const IndexType Nj = m.getNodeResolution(J_DIRECTION);

#ifdef AXOM_USE_RAJA

  RAJA::RangeSegment i_range(0, Ni);
  RAJA::RangeSegment j_range(0, Nj);

  using exec_pol = typename structured_exec<ExecPolicy>::loop2d_policy;
  RAJA::kernel<exec_pol>(
    RAJA::make_tuple(i_range, j_range),
    AXOM_LAMBDA(IndexType i, IndexType j) {
      const IndexType faceID = numIFaces + i + j * ICellResolution;
      kernel(faceID, i, j);
    });

#else

  constexpr bool is_serial = std::is_same<ExecPolicy, axom::SEQ_EXEC>::value;
  AXOM_STATIC_ASSERT(is_serial);

  for(IndexType j = 0; j < Nj; ++j)
  {
    const IndexType offset = numIFaces + j * ICellResolution;
    for(IndexType i = 0; i < Ni; ++i)
    {
      const IndexType faceID = i + offset;
      kernel(faceID, i, j);
    }
  }

#endif
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_J_faces(xargs::ijk, const StructuredMesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 3, "Mesh must be 3D.");

  const IndexType numIFaces = m.getTotalNumFaces(I_DIRECTION);
  const IndexType ICellResolution = m.getCellResolution(I_DIRECTION);
  const IndexType numJFacesInKSlice =
    ICellResolution * m.getNodeResolution(J_DIRECTION);
  const IndexType Ni = ICellResolution;
  const IndexType Nj = m.getNodeResolution(J_DIRECTION);
  const IndexType Nk = m.getCellResolution(K_DIRECTION);

#ifdef AXOM_USE_RAJA

  RAJA::RangeSegment i_range(0, Ni);
  RAJA::RangeSegment j_range(0, Nj);
  RAJA::RangeSegment k_range(0, Nk);

  using exec_pol = typename structured_exec<ExecPolicy>::loop3d_policy;
  RAJA::kernel<exec_pol>(
    RAJA::make_tuple(i_range, j_range, k_range),
    AXOM_LAMBDA(IndexType i, IndexType j, IndexType k) {
      const IndexType jp = j * ICellResolution;
      const IndexType kp = k * numJFacesInKSlice;
      const IndexType faceID = numIFaces + i + jp + kp;
      kernel(faceID, i, j, k);
    });

#else

  constexpr bool is_serial = std::is_same<ExecPolicy, axom::SEQ_EXEC>::value;
  AXOM_STATIC_ASSERT(is_serial);

  for(IndexType k = 0; k < Nk; ++k)
  {
    const IndexType k_offset = k * numJFacesInKSlice + numIFaces;
    for(IndexType j = 0; j < Nj; ++j)
    {
      const IndexType offset = j * ICellResolution + k_offset;
      for(IndexType i = 0; i < Ni; ++i)
      {
        const IndexType faceID = i + offset;
        kernel(faceID, i, j, k);
      }
    }
  }

#endif
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_K_faces(xargs::ijk, const StructuredMesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 3, "Mesh must be 3D.");

  const IndexType numIJFaces =
    m.getTotalNumFaces(I_DIRECTION) + m.getTotalNumFaces(J_DIRECTION);
  const IndexType ICellResolution = m.getCellResolution(I_DIRECTION);
  const IndexType cellKp = m.cellKp();
  const IndexType Ni = ICellResolution;
  const IndexType Nj = m.getCellResolution(J_DIRECTION);
  const IndexType Nk = m.getNodeResolution(K_DIRECTION);

#ifdef AXOM_USE_RAJA

  RAJA::RangeSegment i_range(0, Ni);
  RAJA::RangeSegment j_range(0, Nj);
  RAJA::RangeSegment k_range(0, Nk);

  using exec_pol = typename structured_exec<ExecPolicy>::loop3d_policy;
  RAJA::kernel<exec_pol>(
    RAJA::make_tuple(i_range, j_range, k_range),
    AXOM_LAMBDA(IndexType i, IndexType j, IndexType k) {
      const IndexType jp = j * ICellResolution;
      const IndexType kp = k * cellKp;
      const IndexType faceID = numIJFaces + i + jp + kp;
      kernel(faceID, i, j, k);
    });

#else

  constexpr bool is_serial = std::is_same<ExecPolicy, axom::SEQ_EXEC>::value;
  AXOM_STATIC_ASSERT(is_serial);

  for(IndexType k = 0; k < Nk; ++k)
  {
    const IndexType k_offset = k * cellKp + numIJFaces;
    for(IndexType j = 0; j < Nj; ++j)
    {
      const IndexType offset = j * ICellResolution + k_offset;
      for(IndexType i = 0; i < Ni; ++i)
      {
        const IndexType faceID = i + offset;
        kernel(faceID, i, j, k);
      }
    }
  }

#endif
}

} /* namespace helpers */

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_faces_impl(xargs::index, const Mesh& m, KernelType&& kernel)
{
  const IndexType numFaces = m.getNumberOfFaces();
  axom::for_all<ExecPolicy>(numFaces, std::forward<KernelType>(kernel));
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_faces(xargs::index, const Mesh& m, KernelType&& kernel)
{
  return for_all_faces_impl<ExecPolicy>(xargs::index(),
                                        m,
                                        std::forward<KernelType>(kernel));
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_faces_impl(xargs::nodeids,
                               const StructuredMesh& m,
                               KernelType&& kernel)
{
  const IndexType dimension = m.getDimension();
  const IndexType* offsets = m.getCellNodeOffsetsArray();
  const IndexType cellNodeOffset3 = offsets[3];

  if(dimension == 2)
  {
    const IndexType numIFaces = m.getTotalNumFaces(I_DIRECTION);

    helpers::for_all_I_faces<ExecPolicy>(
      xargs::ij(),
      m,
      AXOM_LAMBDA(IndexType faceID,
                  IndexType AXOM_NOT_USED(i),
                  IndexType AXOM_NOT_USED(j)) {
        IndexType nodes[2];
        nodes[0] = faceID;
        nodes[1] = nodes[0] + cellNodeOffset3;
        kernel(faceID, nodes, 2);
      });

    helpers::for_all_J_faces<ExecPolicy>(
      xargs::ij(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType AXOM_NOT_USED(i), IndexType j) {
        const IndexType shiftedID = faceID - numIFaces;
        IndexType nodes[2];
        nodes[0] = shiftedID + j;
        nodes[1] = nodes[0] + 1;
        kernel(faceID, nodes, 2);
      });
  }
  else
  {
    SLIC_ERROR_IF(dimension != 3,
                  "for_all_faces is only valid for 2 or 3D meshes.");

    const IndexType numIFaces = m.getTotalNumFaces(I_DIRECTION);
    const IndexType numIJFaces = numIFaces + m.getTotalNumFaces(J_DIRECTION);
    const IndexType INodeResolution = m.getNodeResolution(I_DIRECTION);
    const IndexType JNodeResolution = m.getNodeResolution(J_DIRECTION);
    const IndexType KFaceNodeStride =
      m.getCellResolution(I_DIRECTION) + m.getCellResolution(J_DIRECTION) + 1;

    const IndexType cellNodeOffset2 = offsets[2];
    const IndexType cellNodeOffset4 = offsets[4];
    const IndexType cellNodeOffset5 = offsets[5];
    const IndexType cellNodeOffset7 = offsets[7];

    helpers::for_all_I_faces<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType faceID,
                  IndexType AXOM_NOT_USED(i),
                  IndexType AXOM_NOT_USED(j),
                  IndexType k) {
        IndexType nodes[4];
        nodes[0] = faceID + k * INodeResolution;
        nodes[1] = nodes[0] + cellNodeOffset4;
        nodes[2] = nodes[0] + cellNodeOffset7;
        nodes[3] = nodes[0] + cellNodeOffset3;
        kernel(faceID, nodes, 4);
      });

    helpers::for_all_J_faces<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType faceID,
                  IndexType AXOM_NOT_USED(i),
                  IndexType j,
                  IndexType k) {
        const IndexType shiftedID = faceID - numIFaces;
        IndexType nodes[4];
        nodes[0] = shiftedID + j + k * JNodeResolution;
        nodes[1] = nodes[0] + 1;
        nodes[2] = nodes[0] + cellNodeOffset5;
        nodes[3] = nodes[0] + cellNodeOffset4;
        kernel(faceID, nodes, 4);
      });

    helpers::for_all_K_faces<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType faceID,
                  IndexType AXOM_NOT_USED(i),
                  IndexType j,
                  IndexType k) {
        const IndexType shiftedID = faceID - numIJFaces;
        IndexType nodes[4];
        nodes[0] = shiftedID + j + k * KFaceNodeStride;
        nodes[1] = nodes[0] + 1;
        nodes[2] = nodes[0] + cellNodeOffset2;
        nodes[3] = nodes[0] + cellNodeOffset3;
        kernel(faceID, nodes, 4);
      });
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_faces_impl(xargs::nodeids,
                               const UnstructuredMesh<SINGLE_SHAPE>& m,
                               KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getNumberOfFaces() <= 0,
                "No faces in the mesh, perhaps you meant to call "
                  << "UnstructuredMesh::initializeFaceConnectivity first.");

  const IndexType* faces_to_nodes = m.getFaceNodesArray();
  const IndexType num_nodes = m.getNumberOfFaceNodes();

  for_all_faces_impl<ExecPolicy>(
    xargs::index(),
    m,
    AXOM_LAMBDA(IndexType faceID) {
      kernel(faceID, faces_to_nodes + faceID * num_nodes, num_nodes);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_faces_impl(xargs::nodeids,
                               const UnstructuredMesh<MIXED_SHAPE>& m,
                               KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getNumberOfFaces() <= 0,
                "No faces in the mesh, perhaps you meant to call "
                  << "UnstructuredMesh::initializeFaceConnectivity first.");

  const IndexType* faces_to_nodes = m.getFaceNodesArray();
  const IndexType* offsets = m.getFaceNodesOffsetsArray();

  for_all_faces_impl<ExecPolicy>(
    xargs::index(),
    m,
    AXOM_LAMBDA(IndexType faceID) {
      const IndexType num_nodes = offsets[faceID + 1] - offsets[faceID];
      kernel(faceID, faces_to_nodes + offsets[faceID], num_nodes);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_faces(xargs::nodeids, const Mesh& m, KernelType&& kernel)
{
  SLIC_ASSERT(m.getDimension() > 1 && m.getDimension() <= 3);

  if(m.isStructured())
  {
    const StructuredMesh& sm = static_cast<const StructuredMesh&>(m);
    for_all_faces_impl<ExecPolicy>(xargs::nodeids(),
                                   sm,
                                   std::forward<KernelType>(kernel));
  }
  else if(m.hasMixedCellTypes())
  {
    const UnstructuredMesh<MIXED_SHAPE>& um =
      static_cast<const UnstructuredMesh<MIXED_SHAPE>&>(m);
    for_all_faces_impl<ExecPolicy>(xargs::nodeids(),
                                   um,
                                   std::forward<KernelType>(kernel));
  }
  else
  {
    const UnstructuredMesh<SINGLE_SHAPE>& um =
      static_cast<const UnstructuredMesh<SINGLE_SHAPE>&>(m);
    for_all_faces_impl<ExecPolicy>(xargs::nodeids(),
                                   um,
                                   std::forward<KernelType>(kernel));
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_faces_impl(xargs::cellids,
                               const StructuredMesh& m,
                               KernelType&& kernel)
{
  const IndexType ICellResolution = m.getCellResolution(I_DIRECTION);
  const IndexType JCellResolution = m.getCellResolution(J_DIRECTION);
  const IndexType cellJp = m.cellJp();
  const int dimension = m.getDimension();

  if(dimension == 2)
  {
    helpers::for_all_I_faces<ExecPolicy>(
      xargs::ij(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j) {
        IndexType cellIDTwo = i + j * cellJp;
        IndexType cellIDOne = cellIDTwo - 1;
        if(i == 0)
        {
          cellIDOne = cellIDTwo;
          cellIDTwo = -1;
        }
        else if(i == ICellResolution)
        {
          cellIDTwo = -1;
        }

        kernel(faceID, cellIDOne, cellIDTwo);
      });

    helpers::for_all_J_faces<ExecPolicy>(
      xargs::ij(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j) {
        IndexType cellIDTwo = i + j * cellJp;
        IndexType cellIDOne = cellIDTwo - cellJp;
        if(j == 0)
        {
          cellIDOne = cellIDTwo;
          cellIDTwo = -1;
        }
        else if(j == JCellResolution)
        {
          cellIDTwo = -1;
        }

        kernel(faceID, cellIDOne, cellIDTwo);
      });
  }
  else
  {
    SLIC_ERROR_IF(dimension != 3,
                  "for_all_faces only valid for 2 or 3D meshes.");

    const IndexType KCellResolution = m.getCellResolution(K_DIRECTION);
    const IndexType cellKp = m.cellKp();

    helpers::for_all_I_faces<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j, IndexType k) {
        IndexType cellIDTwo = i + j * cellJp + k * cellKp;
        IndexType cellIDOne = cellIDTwo - 1;
        if(i == 0)
        {
          cellIDOne = cellIDTwo;
          cellIDTwo = -1;
        }
        else if(i == ICellResolution)
        {
          cellIDTwo = -1;
        }

        kernel(faceID, cellIDOne, cellIDTwo);
      });

    helpers::for_all_J_faces<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j, IndexType k) {
        IndexType cellIDTwo = i + j * cellJp + k * cellKp;
        IndexType cellIDOne = cellIDTwo - cellJp;
        if(j == 0)
        {
          cellIDOne = cellIDTwo;
          cellIDTwo = -1;
        }
        else if(j == JCellResolution)
        {
          cellIDTwo = -1;
        }

        kernel(faceID, cellIDOne, cellIDTwo);
      });

    helpers::for_all_K_faces<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j, IndexType k) {
        IndexType cellIDTwo = i + j * cellJp + k * cellKp;
        IndexType cellIDOne = cellIDTwo - cellKp;
        if(k == 0)
        {
          cellIDOne = cellIDTwo;
          cellIDTwo = -1;
        }
        else if(k == KCellResolution)
        {
          cellIDTwo = -1;
        }

        kernel(faceID, cellIDOne, cellIDTwo);
      });
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, Topology TOPO, typename KernelType>
inline void for_all_faces_impl(xargs::cellids,
                               const UnstructuredMesh<TOPO>& m,
                               KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getNumberOfFaces() <= 0,
                "No faces in the mesh, perhaps you meant to call "
                  << "UnstructuredMesh::initializeFaceConnectivity first.");

  const IndexType* faces_to_cells = m.getFaceCellsArray();

  for_all_faces_impl<ExecPolicy>(
    xargs::index(),
    m,
    AXOM_LAMBDA(IndexType faceID) {
      const IndexType offset = 2 * faceID;
      kernel(faceID, faces_to_cells[offset], faces_to_cells[offset + 1]);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_faces(xargs::cellids, const Mesh& m, KernelType&& kernel)
{
  SLIC_ASSERT(m.getDimension() > 1 && m.getDimension() <= 3);

  if(m.isStructured())
  {
    const StructuredMesh& sm = static_cast<const StructuredMesh&>(m);
    for_all_faces_impl<ExecPolicy>(xargs::cellids(),
                                   sm,
                                   std::forward<KernelType>(kernel));
  }
  else if(m.hasMixedCellTypes())
  {
    const UnstructuredMesh<MIXED_SHAPE>& um =
      static_cast<const UnstructuredMesh<MIXED_SHAPE>&>(m);
    for_all_faces_impl<ExecPolicy>(xargs::cellids(),
                                   um,
                                   std::forward<KernelType>(kernel));
  }
  else
  {
    const UnstructuredMesh<SINGLE_SHAPE>& um =
      static_cast<const UnstructuredMesh<SINGLE_SHAPE>&>(m);
    for_all_faces_impl<ExecPolicy>(xargs::cellids(),
                                   um,
                                   std::forward<KernelType>(kernel));
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_faces_impl(xargs::coords,
                               const UniformMesh& m,
                               KernelType&& kernel)
{
  constexpr bool NO_COPY = true;

  SLIC_ASSERT(m.getDimension() > 1 && m.getDimension() <= 3);

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

  if(dimension == 2)
  {
    helpers::for_all_I_faces<ExecPolicy>(
      xargs::ij(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j) {
        const IndexType n0 = i + j * nodeJp;
        const IndexType nodeIDs[2] = {n0, n0 + nodeJp};

        double coords[4] = {x0 + i * dx,
                            y0 + j * dy,
                            x0 + i * dx,
                            y0 + (j + 1) * dy};

        numerics::Matrix<double> coordsMatrix(dimension, 2, coords, NO_COPY);
        kernel(faceID, coordsMatrix, nodeIDs);
      });

    helpers::for_all_J_faces<ExecPolicy>(
      xargs::ij(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j) {
        const IndexType n0 = i + j * nodeJp;
        const IndexType nodeIDs[2] = {n0, n0 + 1};

        double coords[4] = {x0 + i * dx,
                            y0 + j * dy,
                            x0 + (i + 1) * dx,
                            y0 + j * dy};

        numerics::Matrix<double> coordsMatrix(dimension, 2, coords, NO_COPY);
        kernel(faceID, coordsMatrix, nodeIDs);
      });
  }
  else
  {
    helpers::for_all_I_faces<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j, IndexType k) {
        const IndexType n0 = i + j * nodeJp + k * nodeKp;
        const IndexType nodeIDs[4] = {n0,
                                      n0 + nodeKp,
                                      n0 + nodeJp + nodeKp,
                                      n0 + nodeJp};

        double coords[12] = {x0 + i * dx,
                             y0 + j * dy,
                             z0 + k * dz,
                             x0 + i * dx,
                             y0 + j * dy,
                             z0 + (k + 1) * dz,
                             x0 + i * dx,
                             y0 + (j + 1) * dy,
                             z0 + (k + 1) * dz,
                             x0 + i * dx,
                             y0 + (j + 1) * dy,
                             z0 + k * dz};

        numerics::Matrix<double> coordsMatrix(dimension, 4, coords, NO_COPY);
        kernel(faceID, coordsMatrix, nodeIDs);
      });

    helpers::for_all_J_faces<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j, IndexType k) {
        const IndexType n0 = i + j * nodeJp + k * nodeKp;
        const IndexType nodeIDs[4] = {n0, n0 + 1, n0 + 1 + nodeKp, n0 + nodeKp};

        double coords[12] = {x0 + i * dx,
                             y0 + j * dy,
                             z0 + k * dz,
                             x0 + (i + 1) * dx,
                             y0 + j * dy,
                             z0 + k * dz,
                             x0 + (i + 1) * dx,
                             y0 + j * dy,
                             z0 + (k + 1) * dz,
                             x0 + i * dx,
                             y0 + j * dy,
                             z0 + (k + 1) * dz};

        numerics::Matrix<double> coordsMatrix(dimension, 4, coords, NO_COPY);
        kernel(faceID, coordsMatrix, nodeIDs);
      });

    helpers::for_all_K_faces<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j, IndexType k) {
        const IndexType n0 = i + j * nodeJp + k * nodeKp;
        const IndexType nodeIDs[4] = {n0, n0 + 1, n0 + 1 + nodeJp, n0 + nodeJp};

        double coords[12] = {x0 + i * dx,
                             y0 + j * dy,
                             z0 + k * dz,
                             x0 + (i + 1) * dx,
                             y0 + j * dy,
                             z0 + k * dz,
                             x0 + (i + 1) * dx,
                             y0 + (j + 1) * dy,
                             z0 + k * dz,
                             x0 + i * dx,
                             y0 + (j + 1) * dy,
                             z0 + k * dz};

        numerics::Matrix<double> coordsMatrix(dimension, 4, coords, NO_COPY);
        kernel(faceID, coordsMatrix, nodeIDs);
      });
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_faces_impl(xargs::coords,
                               const RectilinearMesh& m,
                               KernelType&& kernel)
{
  constexpr bool NO_COPY = true;

  SLIC_ASSERT(m.getDimension() > 1 && m.getDimension() <= 3);

  const int dimension = m.getDimension();
  const IndexType nodeJp = m.nodeJp();
  const IndexType nodeKp = m.nodeKp();
  const double* x = m.getCoordinateArray(X_COORDINATE);
  const double* y = m.getCoordinateArray(Y_COORDINATE);

  if(dimension == 2)
  {
    helpers::for_all_I_faces<ExecPolicy>(
      xargs::ij(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j) {
        const IndexType n0 = i + j * nodeJp;
        const IndexType nodeIDs[2] = {n0, n0 + nodeJp};

        double coords[4] = {x[i], y[j], x[i], y[j + 1]};

        numerics::Matrix<double> coordsMatrix(dimension, 2, coords, NO_COPY);
        kernel(faceID, coordsMatrix, nodeIDs);
      });

    helpers::for_all_J_faces<ExecPolicy>(
      xargs::ij(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j) {
        const IndexType n0 = i + j * nodeJp;
        const IndexType nodeIDs[2] = {n0, n0 + 1};

        double coords[4] = {x[i], y[j], x[i + 1], y[j]};

        numerics::Matrix<double> coordsMatrix(dimension, 2, coords, NO_COPY);
        kernel(faceID, coordsMatrix, nodeIDs);
      });
  }
  else
  {
    const double* z = m.getCoordinateArray(Z_COORDINATE);
    helpers::for_all_I_faces<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j, IndexType k) {
        const IndexType n0 = i + j * nodeJp + k * nodeKp;
        const IndexType nodeIDs[4] = {n0,
                                      n0 + nodeKp,
                                      n0 + nodeJp + nodeKp,
                                      n0 + nodeJp};

        double coords[12] = {x[i],
                             y[j],
                             z[k],
                             x[i],
                             y[j],
                             z[k + 1],
                             x[i],
                             y[j + 1],
                             z[k + 1],
                             x[i],
                             y[j + 1],
                             z[k]};

        numerics::Matrix<double> coordsMatrix(dimension, 4, coords, NO_COPY);
        kernel(faceID, coordsMatrix, nodeIDs);
      });

    helpers::for_all_J_faces<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j, IndexType k) {
        const IndexType n0 = i + j * nodeJp + k * nodeKp;
        const IndexType nodeIDs[4] = {n0, n0 + 1, n0 + 1 + nodeKp, n0 + nodeKp};

        double coords[12] = {x[i],
                             y[j],
                             z[k],
                             x[i + 1],
                             y[j],
                             z[k],
                             x[i + 1],
                             y[j],
                             z[k + 1],
                             x[i],
                             y[j],
                             z[k + 1]};

        numerics::Matrix<double> coordsMatrix(dimension, 4, coords, NO_COPY);
        kernel(faceID, coordsMatrix, nodeIDs);
      });

    helpers::for_all_K_faces<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j, IndexType k) {
        const IndexType n0 = i + j * nodeJp + k * nodeKp;
        const IndexType nodeIDs[4] = {n0, n0 + 1, n0 + 1 + nodeJp, n0 + nodeJp};

        double coords[12] = {x[i],
                             y[j],
                             z[k],
                             x[i + 1],
                             y[j],
                             z[k],
                             x[i + 1],
                             y[j + 1],
                             z[k],
                             x[i],
                             y[j + 1],
                             z[k]};

        numerics::Matrix<double> coordsMatrix(dimension, 4, coords, NO_COPY);
        kernel(faceID, coordsMatrix, nodeIDs);
      });
  }
}

//------------------------------------------------------------------------------
struct for_all_face_nodes_functor
{
  template <typename ExecPolicy, typename MeshType, typename KernelType>
  inline void operator()(ExecPolicy AXOM_NOT_USED(policy),
                         const MeshType& m,
                         KernelType&& kernel) const
  {
    constexpr bool valid_mesh_type = std::is_base_of<Mesh, MeshType>::value;
    AXOM_STATIC_ASSERT(valid_mesh_type);

    for_all_faces_impl<ExecPolicy>(xargs::nodeids(),
                                   m,
                                   std::forward<KernelType>(kernel));
  }
};

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_faces_impl(xargs::coords,
                               const CurvilinearMesh& m,
                               KernelType&& kernel)
{
  SLIC_ASSERT(m.getDimension() > 1 && m.getDimension() <= 3);

  const int dimension = m.getDimension();
  if(dimension == 2)
  {
    for_all_coords<ExecPolicy, 2, 2>(for_all_face_nodes_functor(),
                                     m,
                                     std::forward<KernelType>(kernel));
  }
  else
  {
    for_all_coords<ExecPolicy, 3, 4>(for_all_face_nodes_functor(),
                                     m,
                                     std::forward<KernelType>(kernel));
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType, Topology TOPO>
inline void for_all_faces_impl(xargs::coords,
                               const UnstructuredMesh<TOPO>& m,
                               KernelType&& kernel)
{
  constexpr bool NO_COPY = true;

  SLIC_ASSERT(m.getDimension() > 1 && m.getDimension() <= 3);

  const int dimension = m.getDimension();
  const double* x = m.getCoordinateArray(X_COORDINATE);
  const double* y = m.getCoordinateArray(Y_COORDINATE);

  if(dimension == 2)
  {
    for_all_faces_impl<ExecPolicy>(
      xargs::nodeids(),
      m,
      AXOM_LAMBDA(IndexType faceID, const IndexType* nodeIDs, IndexType numNodes) {
        double coords[2 * MAX_FACE_NODES];
        for(int i = 0; i < numNodes; ++i)
        {
          const IndexType nodeID = nodeIDs[i];
          coords[2 * i] = x[nodeID];
          coords[2 * i + 1] = y[nodeID];
        }

        numerics::Matrix<double> coordsMatrix(dimension, 2, coords, NO_COPY);
        kernel(faceID, coordsMatrix, nodeIDs);
      });
  }
  else
  {
    const double* z = m.getCoordinateArray(Z_COORDINATE);
    for_all_faces_impl<ExecPolicy>(
      xargs::nodeids(),
      m,
      AXOM_LAMBDA(IndexType faceID, const IndexType* nodeIDs, IndexType numNodes) {
        double coords[3 * MAX_FACE_NODES];
        for(int i = 0; i < numNodes; ++i)
        {
          const IndexType nodeID = nodeIDs[i];
          coords[3 * i] = x[nodeID];
          coords[3 * i + 1] = y[nodeID];
          coords[3 * i + 2] = z[nodeID];
        }

        numerics::Matrix<double> coordsMatrix(dimension, 4, coords, NO_COPY);
        kernel(faceID, coordsMatrix, nodeIDs);
      });
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_faces(xargs::coords, const Mesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() <= 1 || m.getDimension() > 3,
                "Invalid dimension");

  if(m.getMeshType() == STRUCTURED_UNIFORM_MESH)
  {
    const UniformMesh& um = static_cast<const UniformMesh&>(m);
    for_all_faces_impl<ExecPolicy>(xargs::coords(),
                                   um,
                                   std::forward<KernelType>(kernel));
  }
  else if(m.getMeshType() == STRUCTURED_RECTILINEAR_MESH)
  {
    const RectilinearMesh& rm = static_cast<const RectilinearMesh&>(m);
    for_all_faces_impl<ExecPolicy>(xargs::coords(),
                                   rm,
                                   std::forward<KernelType>(kernel));
  }
  else if(m.getMeshType() == STRUCTURED_CURVILINEAR_MESH)
  {
    const CurvilinearMesh& cm = static_cast<const CurvilinearMesh&>(m);
    for_all_faces_impl<ExecPolicy>(xargs::coords(),
                                   cm,
                                   std::forward<KernelType>(kernel));
  }
  else if(m.getMeshType() == UNSTRUCTURED_MESH)
  {
    if(m.hasMixedCellTypes())
    {
      const UnstructuredMesh<MIXED_SHAPE>& um =
        static_cast<const UnstructuredMesh<MIXED_SHAPE>&>(m);

      for_all_faces_impl<ExecPolicy>(xargs::coords(),
                                     um,
                                     std::forward<KernelType>(kernel));
    }
    else
    {
      const UnstructuredMesh<SINGLE_SHAPE>& um =
        static_cast<const UnstructuredMesh<SINGLE_SHAPE>&>(m);

      for_all_faces_impl<ExecPolicy>(xargs::coords(),
                                     um,
                                     std::forward<KernelType>(kernel));
    }
  }
  else
  {
    SLIC_ERROR("Unknown mesh type.");
  }
}

} /* namespace internal */
} /* namespace mint     */
} /* namespace axom     */

#endif /* MINT_FOR_ALL_FACES_HPP_ */
