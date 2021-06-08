// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// axom/core includes
#include "axom/config.hpp"

#include "axom/core/Types.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/core/numerics/Matrix.hpp"

// axom/primal includes
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"

// axom/spin includes
#include "axom/spin/BVH.hpp"
#include "axom/spin/UniformGrid.hpp"

#include "axom/spin/internal/linear_bvh/QueryAccessor.hpp"

// axom/mint includes
#include "axom/mint/mesh/Mesh.hpp"
#include "axom/mint/mesh/UniformMesh.hpp"
#include "axom/mint/execution/interface.hpp"
#include "axom/mint/utils/vtk_utils.hpp"

// gtest includes
#include "gtest/gtest.h"

using namespace axom;
namespace xargs = mint::xargs;

// Uncomment the following for debugging
//#define VTK_DEBUG

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
//------------------------------------------------------------------------------
template <typename FloatType>
void dump_ray(const std::string& file,
              FloatType t,
              FloatType x0,
              FloatType nx,
              FloatType y0,
              FloatType ny,
              FloatType z0 = 0.0,
              FloatType nz = 0.0)
{
  std::ofstream ofs(file.c_str());

  ofs << "# vtk DataFile Version 3.0\n";
  ofs << "Ray Data\n";
  ofs << "ASCII\n";
  ofs << "DATASET UNSTRUCTURED_GRID\n";

  ofs << "POINTS 2 double\n";
  ofs << x0 << " " << y0 << " " << z0 << std::endl;

  ofs << (x0 + t * nx) << " ";
  ofs << (y0 + t * ny) << " ";
  ofs << (z0 + t * nz) << std::endl;

  ofs << "CELLS 1 3\n";
  ofs << "2 0 1\n";

  ofs << "CELL_TYPES 1\n";
  ofs << "3\n";

  ofs << "CELL_DATA 1\n";
  ofs << "VECTORS normal double\n";
  ofs << nx << " " << ny << " " << nz << std::endl;

  ofs.close();
}

//------------------------------------------------------------------------------

/*!
 * \brief Give a 2D mesh object, this method generates an array of axis-aligned
 *  bounding boxes corresponding to each constituent cell of the mesh and
 *  a corresponding centroid.
 *
 * \param [in]  mesh pointer to the mesh object.
 * \param [out] aabbs flat array of bounding boxes [xmin,ymin,xmax,ymax....]
 *
 * \note The intent of this method is to generate synthetic input test data
 *  to test the functionality of the BVH.
 *
 * \warning This method allocates aabbs internally. The caller is responsible
 *  for properly deallocating aabbs.
 *
 * \pre aabbs == nullptr
 * \post aabbs != nullptr
 */
template <typename FloatType>
void generate_aabbs2d(const mint::Mesh* mesh, FloatType*& aabbs)
{
  // sanity checks
  EXPECT_TRUE(aabbs == nullptr);

  // calculate some constants
  constexpr int ndims = 2;
  constexpr int stride = 2 * ndims;

  EXPECT_EQ(ndims, mesh->getDimension());

  // allocate output arrays
  const IndexType ncells = mesh->getNumberOfCells();
  aabbs = axom::allocate<FloatType>(ncells * stride);

  using exec_policy = axom::SEQ_EXEC;
  mint::for_all_cells<exec_policy, xargs::coords>(
    mesh,
    AXOM_LAMBDA(IndexType cellIdx,
                numerics::Matrix<double> & coords,
                const IndexType* AXOM_NOT_USED(nodeIds)) {
      primal::BoundingBox<double, 2> range;

      for(IndexType inode = 0; inode < 4; ++inode)
      {
        primal::Point<double, 2> node {coords.getColumn(inode)};

        range.addPoint(node);
      }  // END for all cells nodes

      const IndexType offset = cellIdx * stride;
      aabbs[offset] = range.getMin()[mint::X_COORDINATE];
      aabbs[offset + 1] = range.getMin()[mint::Y_COORDINATE];
      aabbs[offset + 2] = range.getMax()[mint::X_COORDINATE];
      aabbs[offset + 3] = range.getMax()[mint::Y_COORDINATE];
    });

  // post-condition sanity checks
  EXPECT_TRUE(aabbs != nullptr);
}

//------------------------------------------------------------------------------

/*!
 * \brief Give a 3D mesh object, this method generates an array of axis-aligned
 *  bounding boxes (AABBS) corresponding to each constituent cell of the mesh
 *  and a corresponding centroid.
 *
 * \param [in]  mesh pointer to the mesh object.
 * \param [out] aabbs flat array of AABBS [xmin,ymin,zmin,xmax,ymax,zmax....]
 *
 * \note The intent of this method is to generate synthetic input test data
 *  to test the functionality of the BVH.
 *
 * \warning This method allocates aabbs internally. The caller is responsible
 *  for properly deallocating aabbs.
 *
 * \pre aabbs == nullptr
 * \pre xc != nullptr
 * \pre yc != nullptr
 * \pre zc != nullptr
 *
 * \post aabbs != nullptr
 */
template <typename FloatType>
void generate_aabbs3d(const mint::Mesh* mesh, FloatType*& aabbs)
{
  // sanity checks
  EXPECT_TRUE(aabbs == nullptr);

  // calculate some constants
  constexpr int ndims = 3;
  constexpr int stride = 2 * ndims;

  EXPECT_EQ(ndims, mesh->getDimension());

  // allocate output arrays
  const IndexType ncells = mesh->getNumberOfCells();
  aabbs = axom::allocate<FloatType>(ncells * stride);

  using exec_policy = axom::SEQ_EXEC;
  mint::for_all_cells<exec_policy, xargs::coords>(
    mesh,
    AXOM_LAMBDA(IndexType cellIdx,
                numerics::Matrix<double> & coords,
                const IndexType* AXOM_NOT_USED(nodeIds)) {
      primal::BoundingBox<double, 3> range;

      for(IndexType inode = 0; inode < 8; ++inode)
      {
        primal::Point<double, 3> node {coords.getColumn(inode)};

        range.addPoint(node);
      }  // END for all cells nodes

      const IndexType offset = cellIdx * stride;
      aabbs[offset] = range.getMin()[mint::X_COORDINATE];
      aabbs[offset + 1] = range.getMin()[mint::Y_COORDINATE];
      aabbs[offset + 2] = range.getMin()[mint::Z_COORDINATE];
      aabbs[offset + 3] = range.getMax()[mint::X_COORDINATE];
      aabbs[offset + 4] = range.getMax()[mint::Y_COORDINATE];
      aabbs[offset + 5] = range.getMax()[mint::Z_COORDINATE];
    });

  // post-condition sanity checks
  EXPECT_TRUE(aabbs != nullptr);
}

//------------------------------------------------------------------------------

/*!
 * \brief Give a 2D mesh object, this method generates an array of axis-aligned
 *  bounding boxes corresponding to each constituent cell of the mesh and
 *  a corresponding centroid.
 *
 * \param [in]  mesh pointer to the mesh object.
 * \param [out] aabbs flat array of bounding boxes [xmin,ymin,xmax,ymax....]
 * \param [out] xc buffer to store the x-component of the cell centroids
 * \param [out] yc buffer to store the y-component of the cell centroids
 *
 * \note The intent of this method is to generate synthetic input test data
 *  to test the functionality of the BVH.
 *
 * \warning This method allocates aabbs internally. The caller is responsible
 *  for properly deallocating aabbs.
 *
 * \pre aabbs == nullptr
 * \pre xc != nullptr
 * \pre yc != nullptr
 *
 * \post aabbs != nullptr
 */
template <typename FloatType>
void generate_aabbs_and_centroids2d(const mint::Mesh* mesh,
                                    FloatType*& aabbs,
                                    FloatType* xc,
                                    FloatType* yc)
{
  // sanity checks
  EXPECT_TRUE(aabbs == nullptr);
  EXPECT_TRUE(xc != nullptr);
  EXPECT_TRUE(yc != nullptr);

  // calculate some constants
  constexpr double ONE_OVER_4 = 1.f / 4.f;
  constexpr int NDIMS = 2;
  constexpr int STRIDE = 2 * NDIMS;

  EXPECT_EQ(NDIMS, mesh->getDimension());

  // allocate output arrays
  const IndexType ncells = mesh->getNumberOfCells();
  aabbs = axom::allocate<FloatType>(ncells * STRIDE);

#ifdef __ibmxl__

  /* workaround for IBM XL compiler */
  constexpr int NUM_NODES_PER_CELL = 4;
  for(IndexType icell = 0; icell < ncells; ++icell)
  {
    primal::BoundingBox<FloatType, 2> range;

    IndexType nodeIDs[NUM_NODES_PER_CELL];
    const IndexType numNodesPerCell = mesh->getCellNodeIDs(icell, nodeIDs);
    EXPECT_EQ(numNodesPerCell, NUM_NODES_PER_CELL);

    FloatType xsum = 0.0;
    FloatType ysum = 0.0;
    for(IndexType inode = 0; inode < NUM_NODES_PER_CELL; ++inode)
    {
      const IndexType& nodeIdx = nodeIDs[inode];

      double coords[NDIMS];
      mesh->getNode(nodeIdx, coords);

      xsum += coords[mint::X_COORDINATE];
      ysum += coords[mint::Y_COORDINATE];

      range.addPoint(primal::Point<FloatType, 2> {
        static_cast<FloatType>(coords[mint::X_COORDINATE]),
        static_cast<FloatType>(coords[mint::Y_COORDINATE])});
    }  // END for all cell nodes

    xc[icell] = xsum * ONE_OVER_4;
    yc[icell] = ysum * ONE_OVER_4;

    const IndexType offset = icell * STRIDE;
    range.getMin().to_array(aabbs + offset);
    range.getMax().to_array(aabbs + offset + 2);
  }  // END for all cells

#else

  using exec_policy = axom::SEQ_EXEC;
  mint::for_all_cells<exec_policy, xargs::coords>(
    mesh,
    AXOM_LAMBDA(IndexType cellIdx,
                numerics::Matrix<double> & coords,
                const IndexType* AXOM_NOT_USED(nodeIds)) {
      primal::BoundingBox<FloatType, 2> range;

      FloatType xsum = 0.0;
      FloatType ysum = 0.0;

      for(IndexType inode = 0; inode < 4; ++inode)
      {
        const double* node = coords.getColumn(inode);
        xsum += node[mint::X_COORDINATE];
        ysum += node[mint::Y_COORDINATE];

        range.addPoint(primal::Point<FloatType, 2> {
          static_cast<FloatType>(node[mint::X_COORDINATE]),
          static_cast<FloatType>(node[mint::Y_COORDINATE])});
      }  // END for all cells nodes

      xc[cellIdx] = xsum * ONE_OVER_4;
      yc[cellIdx] = ysum * ONE_OVER_4;

      const IndexType offset = cellIdx * STRIDE;
      aabbs[offset] = range.getMin()[mint::X_COORDINATE];
      aabbs[offset + 1] = range.getMin()[mint::Y_COORDINATE];
      aabbs[offset + 2] = range.getMax()[mint::X_COORDINATE];
      aabbs[offset + 3] = range.getMax()[mint::Y_COORDINATE];
    });

#endif  // __ibmxl__

  // post-condition sanity checks
  EXPECT_TRUE(aabbs != nullptr);
}

//------------------------------------------------------------------------------

/*!
 * \brief Give a 3D mesh object, this method generates an array of axis-aligned
 *  bounding boxes (AABBS) corresponding to each constituent cell of the mesh
 *  and a corresponding centroid.
 *
 * \param [in]  mesh pointer to the mesh object.
 * \param [out] aabbs flat array of AABBS [xmin,ymin,zmin,xmax,ymax,zmax....]
 * \param [out] xc buffer to store the x-component of the cell centroids
 * \param [out] yc buffer to store the y-component of the cell centroids
 * \param [out] zc buffer to store the z-component of the cell centroids
 *
 * \note The intent of this method is to generate synthetic input test data
 *  to test the functionality of the BVH.
 *
 * \warning This method allocates aabbs internally. The caller is responsible
 *  for properly deallocating aabbs.
 *
 * \pre aabbs == nullptr
 * \pre xc != nullptr
 * \pre yc != nullptr
 * \pre zc != nullptr
 *
 * \post aabbs != nullptr
 */
template <typename FloatType>
void generate_aabbs_and_centroids3d(const mint::Mesh* mesh,
                                    FloatType*& aabbs,
                                    FloatType* xc,
                                    FloatType* yc,
                                    FloatType* zc)
{
  // sanity checks
  EXPECT_TRUE(aabbs == nullptr);
  EXPECT_TRUE(xc != nullptr);
  EXPECT_TRUE(yc != nullptr);
  EXPECT_TRUE(zc != nullptr);

  // calculate some constants
  constexpr double ONE_OVER_8 = 1.f / 8.f;
  constexpr int NDIMS = 3;
  constexpr int STRIDE = 2 * NDIMS;

  EXPECT_EQ(NDIMS, mesh->getDimension());

  // allocate output arrays
  const IndexType ncells = mesh->getNumberOfCells();
  aabbs = axom::allocate<FloatType>(ncells * STRIDE);

#ifdef __ibmxl__

  /* workaround for IBM XL compiler */
  constexpr int NUM_NODES_PER_CELL = 8;
  for(IndexType icell = 0; icell < ncells; ++icell)
  {
    primal::BoundingBox<FloatType, 3> range;

    IndexType nodeIDs[NUM_NODES_PER_CELL];
    const IndexType numNodesPerCell = mesh->getCellNodeIDs(icell, nodeIDs);
    EXPECT_EQ(numNodesPerCell, NUM_NODES_PER_CELL);

    FloatType xsum = 0.0;
    FloatType ysum = 0.0;
    FloatType zsum = 0.0;
    for(IndexType inode = 0; inode < NUM_NODES_PER_CELL; ++inode)
    {
      const IndexType& nodeIdx = nodeIDs[inode];

      double coords[NDIMS];
      mesh->getNode(nodeIdx, coords);

      xsum += coords[mint::X_COORDINATE];
      ysum += coords[mint::Y_COORDINATE];
      zsum += coords[mint::Z_COORDINATE];

      range.addPoint(primal::Point<FloatType, 3> {
        static_cast<FloatType>(coords[mint::X_COORDINATE]),
        static_cast<FloatType>(coords[mint::Y_COORDINATE]),
        static_cast<FloatType>(coords[mint::Z_COORDINATE])});
    }  // END for all cell nodes

    xc[icell] = xsum * ONE_OVER_8;
    yc[icell] = ysum * ONE_OVER_8;
    zc[icell] = zsum * ONE_OVER_8;

    const IndexType offset = icell * STRIDE;
    range.getMin().to_array(aabbs + offset);
    range.getMax().to_array(aabbs + offset + 3);
  }  // END for all cells

#else

  using exec_policy = axom::SEQ_EXEC;
  mint::for_all_cells<exec_policy, xargs::coords>(
    mesh,
    AXOM_LAMBDA(IndexType cellIdx,
                numerics::Matrix<double> & coords,
                const IndexType* AXOM_NOT_USED(nodeIds)) {
      primal::BoundingBox<FloatType, 3> range;

      FloatType xsum = 0.0;
      FloatType ysum = 0.0;
      FloatType zsum = 0.0;

      for(IndexType inode = 0; inode < 8; ++inode)
      {
        const double* node = coords.getColumn(inode);
        xsum += node[mint::X_COORDINATE];
        ysum += node[mint::Y_COORDINATE];
        zsum += node[mint::Z_COORDINATE];

        range.addPoint(primal::Point<FloatType, 3> {
          static_cast<FloatType>(node[mint::X_COORDINATE]),
          static_cast<FloatType>(node[mint::Y_COORDINATE]),
          static_cast<FloatType>(node[mint::Z_COORDINATE])});
      }  // END for all cells nodes

      xc[cellIdx] = xsum * ONE_OVER_8;
      yc[cellIdx] = ysum * ONE_OVER_8;
      zc[cellIdx] = zsum * ONE_OVER_8;

      const IndexType offset = cellIdx * STRIDE;
      aabbs[offset] = range.getMin()[mint::X_COORDINATE];
      aabbs[offset + 1] = range.getMin()[mint::Y_COORDINATE];
      aabbs[offset + 2] = range.getMin()[mint::Z_COORDINATE];
      aabbs[offset + 3] = range.getMax()[mint::X_COORDINATE];
      aabbs[offset + 4] = range.getMax()[mint::Y_COORDINATE];
      aabbs[offset + 5] = range.getMax()[mint::Z_COORDINATE];
    });

#endif  // __ibmxl__

  // post-condition sanity checks
  EXPECT_TRUE(aabbs != nullptr);
}

//------------------------------------------------------------------------------

/*!
 * \brief Tests the construction of the BVH in 2D by inserting two bounding
 *  boxes in the BVH and ensuring that the bounds of the BVH are as expected.
 */
template <typename ExecSpace, typename FloatType>
void check_build_bvh2d()
{
  constexpr int NUM_BOXES = 2;
  constexpr int NDIMS = 2;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  FloatType* boxes = axom::allocate<FloatType>(8);
  boxes[0] = boxes[1] = 0.;
  boxes[2] = boxes[3] = 1.;
  boxes[4] = boxes[5] = 1.;
  boxes[6] = boxes[7] = 2.;

  spin::BVH<NDIMS, ExecSpace, FloatType> bvh(boxes, NUM_BOXES);

  int allocatorID = bvh.getAllocatorID();
  EXPECT_EQ(allocatorID, axom::execution_space<ExecSpace>::allocatorID());

  bvh.setScaleFactor(1.0);  // i.e., no scaling
  bvh.build();

  FloatType lo[NDIMS];
  FloatType hi[NDIMS];
  bvh.getBounds(lo, hi);

  for(int idim = 0; idim < NDIMS; ++idim)
  {
    EXPECT_DOUBLE_EQ(lo[idim], 0.0);
    EXPECT_DOUBLE_EQ(hi[idim], 2.0);
  }

  axom::deallocate(boxes);
  axom::setDefaultAllocator(current_allocator);
}

//------------------------------------------------------------------------------

/*!
 * \brief Tests the construction of the BVH in 3D by inserting two bounding
 *  boxes in the BVH and ensuring that the bounds of the BVH are as expected.
 */
template <typename ExecSpace, typename FloatType>
void check_build_bvh3d()
{
  constexpr int NUM_BOXES = 2;
  constexpr int NDIMS = 3;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  FloatType* boxes = axom::allocate<FloatType>(12);
  boxes[0] = boxes[1] = boxes[2] = 0.;
  boxes[3] = boxes[4] = boxes[5] = 1.;
  boxes[6] = boxes[7] = boxes[8] = 1.;
  boxes[9] = boxes[10] = boxes[11] = 2.;

  spin::BVH<NDIMS, ExecSpace, FloatType> bvh(boxes, NUM_BOXES);

  int allocatorID = bvh.getAllocatorID();
  EXPECT_EQ(allocatorID, axom::execution_space<ExecSpace>::allocatorID());

  bvh.setScaleFactor(1.0);  // i.e., no scaling
  bvh.build();

  FloatType lo[NDIMS];
  FloatType hi[NDIMS];
  bvh.getBounds(lo, hi);

  for(int idim = 0; idim < NDIMS; ++idim)
  {
    EXPECT_DOUBLE_EQ(lo[idim], 0.0);
    EXPECT_DOUBLE_EQ(hi[idim], 2.0);
  }

  axom::deallocate(boxes);
  axom::setDefaultAllocator(current_allocator);
}

//------------------------------------------------------------------------------
template <typename ExecSpace, typename FloatType>
void check_find_bounding_boxes3d()
{
  constexpr int NDIMS = 3;
  constexpr IndexType N = 2;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  // setup query bounding boxes: both boxes have lower left min at
  // (-1.0,-1.0,-1.0) but different upper right max.
  // The first bounding box is setup such that it intersects
  // 18 bounding boxes of the mesh.
  FloatType* xmin = axom::allocate<FloatType>(N);
  FloatType* ymin = axom::allocate<FloatType>(N);
  FloatType* zmin = axom::allocate<FloatType>(N);
  xmin[0] = xmin[1] = ymin[0] = ymin[1] = zmin[0] = zmin[1] = -1.0;

  FloatType* xmax = axom::allocate<FloatType>(N);
  FloatType* ymax = axom::allocate<FloatType>(N);
  FloatType* zmax = axom::allocate<FloatType>(N);
  xmax[0] = 2.5;
  ymax[0] = 2.5;
  zmax[0] = 1.5;
  xmax[1] = ymax[1] = zmax[1] = -0.5;

  // setup a test mesh (3 x 3 x 3)
  double lo[NDIMS] = {0.0, 0.0, 0.0};
  double hi[NDIMS] = {3.0, 3.0, 3.0};
  mint::UniformMesh mesh(lo, hi, 4, 4, 4);
  const IndexType ncells = mesh.getNumberOfCells();
  FloatType* aabbs = nullptr;
  generate_aabbs3d(&mesh, aabbs);
  EXPECT_TRUE(aabbs != nullptr);

  // construct the BVH
  spin::BVH<NDIMS, ExecSpace, FloatType> bvh(aabbs, ncells);
  bvh.setScaleFactor(1.0);  // i.e., no scaling
  bvh.build();

  // check BVH bounding box
  FloatType min[NDIMS];
  FloatType max[NDIMS];
  bvh.getBounds(min, max);
  for(int i = 0; i < NDIMS; ++i)
  {
    EXPECT_DOUBLE_EQ(min[i], lo[i]);
    EXPECT_DOUBLE_EQ(max[i], hi[i]);
  }

  // traverse the BVH to find the candidates for all the bounding boxes
  IndexType* offsets = axom::allocate<IndexType>(N);
  IndexType* counts = axom::allocate<IndexType>(N);
  IndexType* candidates = nullptr;
  bvh.findBoundingBoxes(offsets,
                        counts,
                        candidates,
                        N,
                        xmin,
                        xmax,
                        ymin,
                        ymax,
                        zmin,
                        zmax);
  EXPECT_TRUE(candidates != nullptr);

  // flag cells that are found by the bounding box ID
  int* iblank = mesh.createField<int>("iblank", mint::CELL_CENTERED);
  mint::for_all_cells<ExecSpace>(
    &mesh,
    AXOM_LAMBDA(IndexType cellIdx) { iblank[cellIdx] = -1; });

  for(int i = 0; i < N; ++i)
  {
    IndexType ncounts = counts[i];
    IndexType offset = offsets[i];
    for(int j = 0; j < ncounts; ++j)
    {
      IndexType idx = candidates[offset + j];
      iblank[idx] = i;
    }  // END for all cells the bounding box intersects
  }    // END for all bounding boxes

  // check answer with results verified manually by inspection
  constexpr int INTERSECTS_BB = 0;
  constexpr int DOES_NOT_INTERSECT_BB = -1;
  constexpr int EXPECTED_BB1_INTERSECTIONS = 18;
  EXPECT_EQ(counts[0], EXPECTED_BB1_INTERSECTIONS);
  EXPECT_EQ(counts[1], 0);

  for(IndexType i = 0; i < ncells; ++i)
  {
    switch(i)
    {
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
    case 13:
    case 14:
    case 15:
    case 16:
    case 17:
      EXPECT_EQ(iblank[i], INTERSECTS_BB);
      break;
    default:
      EXPECT_EQ(iblank[i], DOES_NOT_INTERSECT_BB);
    }
  }  // END for all mesh cells

  // deallocate
  axom::deallocate(offsets);
  axom::deallocate(candidates);
  axom::deallocate(counts);
  axom::deallocate(aabbs);

  axom::deallocate(xmin);
  axom::deallocate(xmax);
  axom::deallocate(ymin);
  axom::deallocate(ymax);
  axom::deallocate(zmin);
  axom::deallocate(zmax);

  axom::setDefaultAllocator(current_allocator);
}

//------------------------------------------------------------------------------

template <typename ExecSpace, typename FloatType>
void check_find_bounding_boxes2d()
{
  constexpr int NDIMS = 2;
  constexpr IndexType N = 2;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  // setup query bounding boxes: both boxes are source at (-1.0,-1.0) but have
  // different max (upper right). The first box is setup to intersect six
  // bounding boxes of the mesh.
  FloatType* xmin = axom::allocate<FloatType>(N);
  FloatType* ymin = axom::allocate<FloatType>(N);
  xmin[0] = xmin[1] = ymin[0] = ymin[1] = -1.0;

  FloatType* xmax = axom::allocate<FloatType>(N);
  FloatType* ymax = axom::allocate<FloatType>(N);
  xmax[0] = 2.5;
  ymax[0] = 1.5;
  xmax[1] = ymax[1] = -0.1;

  // setup a test mesh (3 x 3)
  double lo[NDIMS] = {0.0, 0.0};
  double hi[NDIMS] = {3.0, 3.0};
  mint::UniformMesh mesh(lo, hi, 4, 4);
  const IndexType ncells = mesh.getNumberOfCells();
  FloatType* aabbs = nullptr;
  generate_aabbs2d(&mesh, aabbs);
  EXPECT_TRUE(aabbs != nullptr);

  // construct the BVH
  spin::BVH<NDIMS, ExecSpace, FloatType> bvh(aabbs, ncells);
  bvh.setScaleFactor(1.0);  // i.e., no scaling
  bvh.build();

  // check BVH bounding box
  FloatType min[NDIMS];
  FloatType max[NDIMS];
  bvh.getBounds(min, max);
  for(int i = 0; i < NDIMS; ++i)
  {
    EXPECT_DOUBLE_EQ(min[i], lo[i]);
    EXPECT_DOUBLE_EQ(max[i], hi[i]);
  }

  // traverse the BVH to find the candidates for all the bounding boxes
  IndexType* offsets = axom::allocate<IndexType>(N);
  IndexType* counts = axom::allocate<IndexType>(N);
  IndexType* candidates = nullptr;
  bvh.findBoundingBoxes(offsets, counts, candidates, N, xmin, xmax, ymin, ymax);
  EXPECT_TRUE(candidates != nullptr);

  // flag cells that are found by the bounding box ID
  int* iblank = mesh.createField<int>("iblank", mint::CELL_CENTERED);
  mint::for_all_cells<ExecSpace>(
    &mesh,
    AXOM_LAMBDA(IndexType cellIdx) { iblank[cellIdx] = -1; });

  for(int i = 0; i < N; ++i)
  {
    IndexType ncounts = counts[i];
    IndexType offset = offsets[i];
    for(int j = 0; j < ncounts; ++j)
    {
      IndexType idx = candidates[offset + j];
      iblank[idx] = i;
    }  // END for all cells the bounding boxes intersects with
  }    // END for all bounding boxes

  // check answer with results verified manually by inspection
  constexpr int INTERSECTS_BB = 0;
  constexpr int DOES_NOT_INTERSECT_BB = -1;
  constexpr int EXPECTED_BB1_INTERSECTIONS = 6;
  EXPECT_EQ(counts[0], EXPECTED_BB1_INTERSECTIONS);
  EXPECT_EQ(counts[1], 0);

  for(IndexType i = 0; i < ncells; ++i)
  {
    if(i == 6 || i == 7 || i == 8)
    {
      EXPECT_EQ(iblank[i], DOES_NOT_INTERSECT_BB);
    }
    else
    {
      EXPECT_EQ(iblank[i], INTERSECTS_BB);
    }

  }  // END for all mesh cells

  // deallocate
  axom::deallocate(offsets);
  axom::deallocate(candidates);
  axom::deallocate(counts);
  axom::deallocate(aabbs);

  axom::deallocate(xmin);
  axom::deallocate(xmax);
  axom::deallocate(ymin);
  axom::deallocate(ymax);

  axom::setDefaultAllocator(current_allocator);
}

//------------------------------------------------------------------------------
template <typename ExecSpace, typename FloatType>
void check_find_rays3d()
{
  constexpr int NDIMS = 3;
  constexpr IndexType N = 2;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  // setup query rays: both rays are source at (-1.0,-1.0) but point in
  // opposite directions. The first ray is setup such that it intersects
  // three bounding boxes of the mesh.
  FloatType* x0 = axom::allocate<FloatType>(N);
  FloatType* y0 = axom::allocate<FloatType>(N);
  FloatType* z0 = axom::allocate<FloatType>(N);
  x0[0] = x0[1] = y0[0] = y0[1] = z0[0] = z0[1] = -1.0;

  FloatType* nx = axom::allocate<FloatType>(N);
  FloatType* ny = axom::allocate<FloatType>(N);
  FloatType* nz = axom::allocate<FloatType>(N);
  nx[0] = ny[0] = nz[0] = 1.0;
  nx[1] = ny[1] = nz[1] = -1.0;

#ifdef VTK_DEBUG
  FloatType t = 5.0;
  dump_ray("ray0.vtk", t, x0[0], nx[0], y0[0], ny[0], z0[0], nz[0]);
  dump_ray("ray1.vtk", t, x0[1], nx[1], y0[1], ny[1], z0[1], nz[1]);
#endif

  // setup a test mesh
  double lo[NDIMS] = {0.0, 0.0, 0.0};
  double hi[NDIMS] = {3.0, 3.0, 3.0};
  mint::UniformMesh mesh(lo, hi, 4, 4, 4);
  const IndexType ncells = mesh.getNumberOfCells();
  FloatType* aabbs = nullptr;
  generate_aabbs3d(&mesh, aabbs);
  EXPECT_TRUE(aabbs != nullptr);

  // construct the BVH
  spin::BVH<NDIMS, ExecSpace, FloatType> bvh(aabbs, ncells);
  bvh.setScaleFactor(1.0);  // i.e., no scaling
  bvh.build();

  // check BVH bounding box
  FloatType min[NDIMS];
  FloatType max[NDIMS];
  bvh.getBounds(min, max);
  for(int i = 0; i < NDIMS; ++i)
  {
    EXPECT_DOUBLE_EQ(min[i], lo[i]);
    EXPECT_DOUBLE_EQ(max[i], hi[i]);
  }

  // traverse the BVH to find the candidates for all the centroids
  IndexType* offsets = axom::allocate<IndexType>(N);
  IndexType* counts = axom::allocate<IndexType>(N);
  IndexType* candidates = nullptr;
  bvh.findRays(offsets, counts, candidates, N, x0, nx, y0, ny, z0, nz);
  EXPECT_TRUE(candidates != nullptr);

  // flag cells that are found by the ray ID
  int* iblank = mesh.createField<int>("iblank", mint::CELL_CENTERED);
  mint::for_all_cells<ExecSpace>(
    &mesh,
    AXOM_LAMBDA(IndexType cellIdx) { iblank[cellIdx] = -1; });

  for(int i = 0; i < N; ++i)
  {
    IndexType ncounts = counts[i];
    IndexType offset = offsets[i];
    for(int j = 0; j < ncounts; ++j)
    {
      IndexType idx = candidates[offset + j];
      iblank[idx] = i;
    }  // END for all cells the ray hits
  }    // END for all rays

#ifdef VTK_DEBUG
  mint::write_vtk(&mesh, "uniform_mesh.vtk");
#endif

  // check answer with results verified manually by inspection
  constexpr int INTERSECTS_RAY = 0;
  constexpr int DOES_NOT_INTERSECT_RAY = -1;
  constexpr int EXPECTED_RAY1_INTERSECTIONS = 15;
  EXPECT_EQ(counts[0], EXPECTED_RAY1_INTERSECTIONS);
  EXPECT_EQ(counts[1], 0);

  for(IndexType i = 0; i < ncells; ++i)
  {
    switch(i)
    {
    case 0:
    case 1:
    case 3:
    case 4:
    case 9:
    case 10:
    case 12:
    case 13:
    case 14:
    case 16:
    case 17:
    case 22:
    case 23:
    case 25:
    case 26:
      EXPECT_EQ(iblank[i], INTERSECTS_RAY);
      break;
    default:
      EXPECT_EQ(iblank[i], DOES_NOT_INTERSECT_RAY);
    }
  }  // END for all mesh cells

  // deallocate
  axom::deallocate(offsets);
  axom::deallocate(candidates);
  axom::deallocate(counts);
  axom::deallocate(aabbs);

  axom::deallocate(x0);
  axom::deallocate(nx);
  axom::deallocate(y0);
  axom::deallocate(ny);
  axom::deallocate(z0);
  axom::deallocate(nz);

  axom::setDefaultAllocator(current_allocator);
}

//------------------------------------------------------------------------------

template <typename ExecSpace, typename FloatType>
void check_find_rays2d()
{
  constexpr int NDIMS = 2;
  constexpr IndexType N = 2;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  // setup query rays: both rays are source at (-1.0,-1.0) but point in
  // opposite directions. The first ray is setup such that it intersects
  // seven bounding boxes of the mesh.
  FloatType* x0 = axom::allocate<FloatType>(N);
  FloatType* y0 = axom::allocate<FloatType>(N);
  x0[0] = x0[1] = y0[0] = y0[1] = -1.0;

  FloatType* nx = axom::allocate<FloatType>(N);
  FloatType* ny = axom::allocate<FloatType>(N);
  nx[0] = ny[0] = 1.0;
  nx[1] = ny[1] = -1.0;

#ifdef VTK_DEBUG
  FloatType t = 5.0;
  dump_ray("ray0.vtk", t, x0[0], nx[0], y0[0], ny[0]);
  dump_ray("ray1.vtk", t, x0[1], nx[1], y0[1], ny[1]);
#endif

  // setup a test mesh
  double lo[NDIMS] = {0.0, 0.0};
  double hi[NDIMS] = {3.0, 3.0};
  mint::UniformMesh mesh(lo, hi, 4, 4);
  const IndexType ncells = mesh.getNumberOfCells();
  FloatType* aabbs = nullptr;
  generate_aabbs2d(&mesh, aabbs);
  EXPECT_TRUE(aabbs != nullptr);

  // construct the BVH
  spin::BVH<NDIMS, ExecSpace, FloatType> bvh(aabbs, ncells);
  bvh.setScaleFactor(1.0);  // i.e., no scaling
  bvh.build();

  // check BVH bounding box
  FloatType min[NDIMS];
  FloatType max[NDIMS];
  bvh.getBounds(min, max);
  for(int i = 0; i < NDIMS; ++i)
  {
    EXPECT_DOUBLE_EQ(min[i], lo[i]);
    EXPECT_DOUBLE_EQ(max[i], hi[i]);
  }

  // traverse the BVH to find the candidates for all the centroids
  IndexType* offsets = axom::allocate<IndexType>(N);
  IndexType* counts = axom::allocate<IndexType>(N);
  IndexType* candidates = nullptr;
  bvh.findRays(offsets, counts, candidates, N, x0, nx, y0, ny);
  EXPECT_TRUE(candidates != nullptr);

  // flag cells that are found by the ray ID
  int* iblank = mesh.createField<int>("iblank", mint::CELL_CENTERED);
  mint::for_all_cells<ExecSpace>(
    &mesh,
    AXOM_LAMBDA(IndexType cellIdx) { iblank[cellIdx] = -1; });

  for(int i = 0; i < N; ++i)
  {
    IndexType ncounts = counts[i];
    IndexType offset = offsets[i];
    for(int j = 0; j < ncounts; ++j)
    {
      IndexType idx = candidates[offset + j];
      iblank[idx] = i;
    }  // END for all cells the ray hits
  }    // END for all rays

#ifdef VTK_DEBUG
  mint::write_vtk(&mesh, "uniform_mesh.vtk");
#endif

  // check answer with results verified manually by inspection
  constexpr int INTERSECTS_RAY = 0;
  constexpr int DOES_NOT_INTERSECT_RAY = -1;
  constexpr int EXPECTED_RAY1_INTERSECTIONS = 7;
  EXPECT_EQ(counts[0], EXPECTED_RAY1_INTERSECTIONS);
  EXPECT_EQ(counts[1], 0);

  for(IndexType i = 0; i < ncells; ++i)
  {
    if(i == 2 || i == 6)
    {
      EXPECT_EQ(iblank[i], DOES_NOT_INTERSECT_RAY);
    }
    else
    {
      EXPECT_EQ(iblank[i], INTERSECTS_RAY);
    }

  }  // END for all mesh cells

  // deallocate
  axom::deallocate(offsets);
  axom::deallocate(candidates);
  axom::deallocate(counts);
  axom::deallocate(aabbs);

  axom::deallocate(x0);
  axom::deallocate(nx);
  axom::deallocate(y0);
  axom::deallocate(ny);

  axom::setDefaultAllocator(current_allocator);
}

//------------------------------------------------------------------------------

/*!
 * \brief Tests the find algorithm of the BVH in 3D.
 *
 *  A uniform mesh is used for this test where the query points are generated
 *  by taking the centroids of the constituent cells of the mesh. Since the
 *  mesh is a uniform, cartesian mesh the find algorithm should return exactly
 *  one candidate for each query point corresponding to the cell on the mesh
 *  that generated the centroid. This property is checked by using the
 *  spin::UniformGrid class.
 *
 *  In addition, the test shifts the points by an offset to ensure that points
 *  outside the mesh return no candidate.
 *
 */
template <typename ExecSpace, typename FloatType>
void check_find_points3d()
{
  constexpr int NDIMS = 3;
  constexpr IndexType N = 4;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  using PointType = primal::Point<double, NDIMS>;

  double lo[NDIMS] = {0.0, 0.0, 0.0};
  double hi[NDIMS] = {3.0, 3.0, 3.0};
  int res[NDIMS] = {N - 1, N - 1, N - 1};

  mint::UniformMesh mesh(lo, hi, N, N, N);
  FloatType* xc = mesh.createField<FloatType>("xc", mint::CELL_CENTERED);
  FloatType* yc = mesh.createField<FloatType>("yc", mint::CELL_CENTERED);
  FloatType* zc = mesh.createField<FloatType>("zc", mint::CELL_CENTERED);
  const IndexType ncells = mesh.getNumberOfCells();

  FloatType* aabbs = nullptr;
  generate_aabbs_and_centroids3d(&mesh, aabbs, xc, yc, zc);

  // construct the BVH
  spin::BVH<NDIMS, ExecSpace, FloatType> bvh(aabbs, ncells);
  bvh.setScaleFactor(1.0);  // i.e., no scaling
  bvh.build();

  FloatType min[NDIMS];
  FloatType max[NDIMS];
  bvh.getBounds(min, max);
  for(int i = 0; i < NDIMS; ++i)
  {
    EXPECT_DOUBLE_EQ(min[i], lo[i]);
    EXPECT_DOUBLE_EQ(max[i], hi[i]);
  }

  // traverse the BVH to find the candidates for all the centroids
  IndexType* offsets = axom::allocate<IndexType>(ncells);
  IndexType* counts = axom::allocate<IndexType>(ncells);
  IndexType* candidates = nullptr;
  bvh.findPoints(offsets, counts, candidates, ncells, xc, yc, zc);

  EXPECT_TRUE(candidates != nullptr);

  spin::UniformGrid<IndexType, NDIMS> ug(lo, hi, res);

  for(IndexType i = 0; i < ncells; ++i)
  {
    PointType q = PointType::make_point(xc[i], yc[i], zc[i]);
    const int donorCellIdx = ug.getBinIndex(q);
    EXPECT_EQ(counts[i], 1);
    EXPECT_EQ(donorCellIdx, candidates[offsets[i]]);
  }  // END for all cell centroids

  axom::deallocate(candidates);

  // check points that are outside by shifting the query points
  constexpr double OFFSET = 10.0;
  for(IndexType i = 0; i < ncells; ++i)
  {
    xc[i] += OFFSET;
    yc[i] += OFFSET;
    zc[i] += OFFSET;
  }

  bvh.findPoints(offsets, counts, candidates, ncells, xc, yc, zc);

  for(IndexType i = 0; i < ncells; ++i)
  {
    EXPECT_EQ(counts[i], 0);
  }

  axom::deallocate(offsets);
  axom::deallocate(candidates);
  axom::deallocate(counts);
  axom::deallocate(aabbs);

  axom::setDefaultAllocator(current_allocator);
}

//------------------------------------------------------------------------------

/*!
 * \brief Tests the find algorithm of the BVH in 2D.
 *
 *  A uniform mesh is used for this test where the query points are generated
 *  by taking the centroids of the constituent cells of the mesh. Since the
 *  mesh is a uniform, cartesian mesh the find algorithm should return exactly
 *  one candidate for each query point corresponding to the cell on the mesh
 *  that generated the centroid. This property is checked by using the
 *  spin::UniformGrid class.
 *
 *  In addition, the test shifts the points by an offset to ensure that points
 *  outside the mesh return no candidate.
 *
 */
template <typename ExecSpace, typename FloatType>
void check_find_points2d()
{
  constexpr int NDIMS = 2;
  constexpr IndexType N = 4;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  using PointType = primal::Point<double, NDIMS>;

  double lo[NDIMS] = {0.0, 0.0};
  double hi[NDIMS] = {3.0, 3.0};
  int res[NDIMS] = {N - 1, N - 1};

  mint::UniformMesh mesh(lo, hi, N, N);
  FloatType* xc = mesh.createField<FloatType>("xc", mint::CELL_CENTERED);
  FloatType* yc = mesh.createField<FloatType>("yc", mint::CELL_CENTERED);
  const IndexType ncells = mesh.getNumberOfCells();

  FloatType* aabbs = nullptr;
  generate_aabbs_and_centroids2d(&mesh, aabbs, xc, yc);

  // construct the BVH
  spin::BVH<NDIMS, ExecSpace, FloatType> bvh(aabbs, ncells);
  bvh.setScaleFactor(1.0);  // i.e., no scaling
  bvh.build();

  FloatType min[NDIMS];
  FloatType max[NDIMS];
  bvh.getBounds(min, max);
  for(int i = 0; i < NDIMS; ++i)
  {
    EXPECT_DOUBLE_EQ(min[i], lo[i]);
    EXPECT_DOUBLE_EQ(max[i], hi[i]);
  }

  // traverse the BVH to find the candidates for all the centroids
  IndexType* offsets = axom::allocate<IndexType>(ncells);
  IndexType* counts = axom::allocate<IndexType>(ncells);
  IndexType* candidates = nullptr;
  bvh.findPoints(offsets, counts, candidates, ncells, xc, yc);

  EXPECT_TRUE(candidates != nullptr);

  spin::UniformGrid<IndexType, NDIMS> ug(lo, hi, res);

  for(IndexType i = 0; i < ncells; ++i)
  {
    PointType q = PointType::make_point(xc[i], yc[i]);
    const int donorCellIdx = ug.getBinIndex(q);
    EXPECT_EQ(counts[i], 1);
    EXPECT_EQ(donorCellIdx, candidates[offsets[i]]);
  }  // END for all cell centroids

  axom::deallocate(candidates);

  // check points that are outside by shifting the query points
  constexpr double OFFSET = 10.0;
  for(IndexType i = 0; i < ncells; ++i)
  {
    xc[i] += OFFSET;
    yc[i] += OFFSET;
  }

  bvh.findPoints(offsets, counts, candidates, ncells, xc, yc);

  for(IndexType i = 0; i < ncells; ++i)
  {
    EXPECT_EQ(counts[i], 0);
  }

  axom::deallocate(offsets);
  axom::deallocate(candidates);
  axom::deallocate(counts);
  axom::deallocate(aabbs);

  axom::setDefaultAllocator(current_allocator);
}

//------------------------------------------------------------------------------

/*!
 * \brief Checks that the BVH behaves properly when user supplies a single box.
 *
 *  The Test inserts a single bounding box that spans [0,1] x [0,1] and tests
 *  that the centroid xc=(0.5,0.5) can be found. Moreover, it shifts the point
 *  by an offset and ensures that a point that is outside is not found.
 */
template <typename ExecSpace, typename FloatType>
void check_single_box2d()
{
  constexpr int NUM_BOXES = 1;
  constexpr int NDIMS = 2;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  // single bounding box in [0,1] x [0,1]
  FloatType* boxes = axom::allocate<FloatType>(4);
  boxes[0] = boxes[1] = 0.;
  boxes[2] = boxes[3] = 1.;

  // construct a BVH with a single box
  spin::BVH<NDIMS, ExecSpace, FloatType> bvh(boxes, NUM_BOXES);
  bvh.setScaleFactor(1.0);  // i.e., no scaling
  bvh.build();

  // check the bounds -- should match the bounds of the input bounding box
  FloatType lo[NDIMS];
  FloatType hi[NDIMS];
  bvh.getBounds(lo, hi);

  for(int idim = 0; idim < NDIMS; ++idim)
  {
    EXPECT_DOUBLE_EQ(lo[idim], 0.0);
    EXPECT_DOUBLE_EQ(hi[idim], 1.0);
  }

  // run the find algorithm w/ the centroid of the bounding box as input.
  // Should return one and only one candidate that corresponds to the
  // single bounding box.
  FloatType* xc = axom::allocate<FloatType>(NUM_BOXES);
  FloatType* yc = axom::allocate<FloatType>(NUM_BOXES);
  xc[0] = yc[0] = 0.5;

  IndexType* offsets = axom::allocate<IndexType>(NUM_BOXES);
  IndexType* counts = axom::allocate<IndexType>(NUM_BOXES);
  IndexType* candidates = nullptr;
  bvh.findPoints(offsets, counts, candidates, NUM_BOXES, xc, yc);
  EXPECT_TRUE(candidates != nullptr);
  EXPECT_EQ(counts[0], 1);
  EXPECT_EQ(0, candidates[offsets[0]]);
  axom::deallocate(candidates);

  // shift centroid outside of the BVH, should return no candidates.
  xc[0] += 10.0;
  yc[0] += 10.0;
  bvh.findPoints(offsets, counts, candidates, NUM_BOXES, xc, yc);
  EXPECT_EQ(counts[0], 0);

  axom::deallocate(xc);
  axom::deallocate(yc);
  axom::deallocate(boxes);
  axom::deallocate(offsets);
  axom::deallocate(counts);
  axom::deallocate(candidates);
  axom::setDefaultAllocator(current_allocator);
}

//------------------------------------------------------------------------------
/*!
 * \brief Checks that the BVH behaves properly when user supplies a single box.
 *
 *  The Test inserts a single bounding box that spans [0,1] x [0,1] x [0,1] and
 *  tests that the centroid xc=(0.5,0.5) can be found. Moreover, it shifts the
 *  point by an offset and ensures that a point that is outside is not found.
 */
template <typename ExecSpace, typename FloatType>
void check_single_box3d()
{
  constexpr int NUM_BOXES = 1;
  constexpr int NDIMS = 3;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  // single bounding box in [0,1] x [0,1] x [0,1]
  FloatType* boxes = axom::allocate<FloatType>(6);
  boxes[0] = boxes[1] = boxes[2] = 0.;
  boxes[3] = boxes[4] = boxes[5] = 1.;

  // construct a BVH with a single box
  spin::BVH<NDIMS, ExecSpace, FloatType> bvh(boxes, NUM_BOXES);
  bvh.setScaleFactor(1.0);  // i.e., no scaling
  bvh.build();

  // check the bounds -- should match the bounds of the input bounding box
  FloatType lo[NDIMS];
  FloatType hi[NDIMS];
  bvh.getBounds(lo, hi);

  for(int idim = 0; idim < NDIMS; ++idim)
  {
    EXPECT_DOUBLE_EQ(lo[idim], 0.0);
    EXPECT_DOUBLE_EQ(hi[idim], 1.0);
  }

  // run the find algorithm w/ the centroid of the bounding box as input.
  // Should return one and only one candidate that corresponds to the
  // single bounding box.
  FloatType* xc = axom::allocate<FloatType>(NUM_BOXES);
  FloatType* yc = axom::allocate<FloatType>(NUM_BOXES);
  FloatType* zc = axom::allocate<FloatType>(NUM_BOXES);
  xc[0] = yc[0] = zc[0] = 0.5;

  IndexType* offsets = axom::allocate<IndexType>(NUM_BOXES);
  IndexType* counts = axom::allocate<IndexType>(NUM_BOXES);
  IndexType* candidates = nullptr;
  bvh.findPoints(offsets, counts, candidates, NUM_BOXES, xc, yc, zc);
  EXPECT_TRUE(candidates != nullptr);
  EXPECT_EQ(counts[0], 1);
  EXPECT_EQ(0, candidates[offsets[0]]);
  axom::deallocate(candidates);

  // shift centroid outside of the BVH, should return no candidates.
  xc[0] += 10.0;
  yc[0] += 10.0;
  bvh.findPoints(offsets, counts, candidates, NUM_BOXES, xc, yc, zc);
  EXPECT_EQ(counts[0], 0);

  axom::deallocate(xc);
  axom::deallocate(yc);
  axom::deallocate(zc);
  axom::deallocate(boxes);
  axom::deallocate(offsets);
  axom::deallocate(counts);
  axom::deallocate(candidates);
  axom::setDefaultAllocator(current_allocator);
}

} /* end unnamed namespace */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(spin_bvh, query_point_accessor)
{
  constexpr double VAL = 42.0;
  constexpr IndexType ID = 0;

  double x[] = {VAL};
  double y[] = {VAL + 1.5};
  double z[] = {VAL + 2.5};

  namespace bvh = axom::spin::internal::linear_bvh;
  using QueryAccessor2D = bvh::QueryAccessor<2, double>;
  using QueryAccessor3D = bvh::QueryAccessor<3, double>;
  using Point2D = primal::Vector<double, 2>;
  using Point3D = primal::Vector<double, 3>;

  Point2D test_point2d;
  QueryAccessor2D::getPoint(test_point2d, ID, x, y, nullptr);
  EXPECT_DOUBLE_EQ(test_point2d[0], x[ID]);
  EXPECT_DOUBLE_EQ(test_point2d[1], y[ID]);

  Point3D test_point3d;
  QueryAccessor3D::getPoint(test_point3d, 0, x, y, z);
  EXPECT_DOUBLE_EQ(test_point3d[0], x[ID]);
  EXPECT_DOUBLE_EQ(test_point3d[1], y[ID]);
  EXPECT_DOUBLE_EQ(test_point3d[2], z[ID]);
}

//------------------------------------------------------------------------------
TEST(spin_bvh, query_ray_accessor)
{
  constexpr double VAL = 42.0;
  constexpr IndexType ID = 0;

  double x[] = {VAL};
  double y[] = {VAL + 1.5};
  double z[] = {VAL + 2.5};

  double nx[] = {VAL};
  double ny[] = {VAL + 1.5};
  double nz[] = {VAL + 2.5};

  namespace bvh = axom::spin::internal::linear_bvh;
  using QueryAccessor2D = bvh::QueryAccessor<2, double>;
  using QueryAccessor3D = bvh::QueryAccessor<3, double>;
  using Ray2D = primal::Vector<double, 4>;
  using Ray3D = primal::Vector<double, 6>;

  Ray2D ray2d;
  QueryAccessor2D::getRay(ray2d, ID, x, nx, y, ny, nullptr, nullptr);
  EXPECT_DOUBLE_EQ(ray2d[0], x[ID]);
  EXPECT_DOUBLE_EQ(ray2d[1], y[ID]);

  EXPECT_DOUBLE_EQ(ray2d[2], nx[ID]);
  EXPECT_DOUBLE_EQ(ray2d[3], ny[ID]);

  Ray3D ray3d;
  QueryAccessor3D::getRay(ray3d, ID, x, nx, y, ny, z, nz);
  EXPECT_DOUBLE_EQ(ray3d[0], x[ID]);
  EXPECT_DOUBLE_EQ(ray3d[1], y[ID]);
  EXPECT_DOUBLE_EQ(ray3d[2], z[ID]);

  EXPECT_DOUBLE_EQ(ray3d[3], nx[ID]);
  EXPECT_DOUBLE_EQ(ray3d[4], ny[ID]);
  EXPECT_DOUBLE_EQ(ray3d[5], nz[ID]);
}

TEST(spin_bvh, query_bounding_box_accessor)
{
  constexpr double VAL = 42.0;
  constexpr IndexType ID = 0;

  double xmin[] = {VAL};
  double ymin[] = {VAL + 1.5};
  double zmin[] = {VAL + 2.5};

  double xmax[] = {VAL};
  double ymax[] = {VAL + 1.5};
  double zmax[] = {VAL + 2.5};

  namespace bvh = axom::spin::internal::linear_bvh;
  using QueryAccessor2D = bvh::QueryAccessor<2, double>;
  using QueryAccessor3D = bvh::QueryAccessor<3, double>;
  using BoundingBox2D = primal::BoundingBox<double, 2>;
  using BoundingBox3D = primal::BoundingBox<double, 3>;

  BoundingBox2D box2D;
  QueryAccessor2D::getBoundingBox(box2D, ID, xmin, xmax, ymin, ymax, nullptr, nullptr);

  EXPECT_DOUBLE_EQ(box2D.getMin()[0], xmin[ID]);
  EXPECT_DOUBLE_EQ(box2D.getMin()[1], ymin[ID]);

  EXPECT_DOUBLE_EQ(box2D.getMax()[0], xmax[ID]);
  EXPECT_DOUBLE_EQ(box2D.getMax()[1], ymax[ID]);

  BoundingBox3D box3D;
  QueryAccessor3D::getBoundingBox(box3D, ID, xmin, xmax, ymin, ymax, zmin, zmax);

  EXPECT_DOUBLE_EQ(box3D.getMin()[0], xmin[ID]);
  EXPECT_DOUBLE_EQ(box3D.getMin()[1], ymin[ID]);
  EXPECT_DOUBLE_EQ(box3D.getMin()[2], zmin[ID]);

  EXPECT_DOUBLE_EQ(box3D.getMax()[0], xmax[ID]);
  EXPECT_DOUBLE_EQ(box3D.getMax()[1], ymax[ID]);
  EXPECT_DOUBLE_EQ(box3D.getMax()[2], zmax[ID]);
}

//------------------------------------------------------------------------------
TEST(spin_bvh, construct2D_sequential)
{
  check_build_bvh2d<axom::SEQ_EXEC, double>();
  check_build_bvh2d<axom::SEQ_EXEC, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, construct3D_sequential)
{
  check_build_bvh3d<axom::SEQ_EXEC, double>();
  check_build_bvh3d<axom::SEQ_EXEC, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, find_bounding_boxes_3d_sequential)
{
  check_find_bounding_boxes3d<axom::SEQ_EXEC, double>();
  check_find_bounding_boxes3d<axom::SEQ_EXEC, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, find_bounding_boxes_2d_sequential)
{
  check_find_bounding_boxes2d<axom::SEQ_EXEC, double>();
  check_find_bounding_boxes2d<axom::SEQ_EXEC, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, find_rays_3d_sequential)
{
  check_find_rays3d<axom::SEQ_EXEC, double>();
  check_find_rays3d<axom::SEQ_EXEC, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, find_rays_2d_sequential)
{
  check_find_rays2d<axom::SEQ_EXEC, double>();
  check_find_rays2d<axom::SEQ_EXEC, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, find_points_3d_sequential)
{
  check_find_points3d<axom::SEQ_EXEC, double>();
  check_find_points3d<axom::SEQ_EXEC, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, find_points_2d_sequential)
{
  check_find_points2d<axom::SEQ_EXEC, double>();
  check_find_points2d<axom::SEQ_EXEC, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, single_box2d_sequential)
{
  check_single_box2d<axom::SEQ_EXEC, double>();
  check_single_box2d<axom::SEQ_EXEC, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, single_box3d_sequential)
{
  check_single_box3d<axom::SEQ_EXEC, double>();
  check_single_box3d<axom::SEQ_EXEC, float>();
}

//------------------------------------------------------------------------------
#ifdef AXOM_USE_OPENMP

TEST(spin_bvh, construct2D_omp)
{
  check_build_bvh2d<axom::OMP_EXEC, double>();
  check_build_bvh2d<axom::OMP_EXEC, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, construct3D_omp)
{
  check_build_bvh3d<axom::OMP_EXEC, double>();
  check_build_bvh3d<axom::OMP_EXEC, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, find_bounding_boxes_3d_omp)
{
  check_find_bounding_boxes3d<axom::OMP_EXEC, double>();
  check_find_bounding_boxes3d<axom::OMP_EXEC, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, find_bounding_boxes_2d_omp)
{
  check_find_bounding_boxes2d<axom::OMP_EXEC, double>();
  check_find_bounding_boxes2d<axom::OMP_EXEC, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, find_rays_3d_omp)
{
  check_find_rays3d<axom::OMP_EXEC, double>();
  check_find_rays3d<axom::OMP_EXEC, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, find_rays_2d_omp)
{
  check_find_rays2d<axom::OMP_EXEC, double>();
  check_find_rays2d<axom::OMP_EXEC, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, find_points_3d_omp)
{
  check_find_points3d<axom::OMP_EXEC, double>();
  check_find_points3d<axom::OMP_EXEC, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, find_points_2d_omp)
{
  check_find_points2d<axom::OMP_EXEC, double>();
  check_find_points2d<axom::OMP_EXEC, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, single_box2d_omp)
{
  check_single_box2d<axom::OMP_EXEC, double>();
  check_single_box2d<axom::OMP_EXEC, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, single_box3d_omp)
{
  check_single_box3d<axom::OMP_EXEC, double>();
}

#endif

//------------------------------------------------------------------------------
#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)

AXOM_CUDA_TEST(spin_bvh, construct2D_cuda)
{
  constexpr int BLOCK_SIZE = 256;
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;

  check_build_bvh2d<exec, double>();
  check_build_bvh2d<exec, float>();
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(spin_bvh, construct3D_cuda)
{
  constexpr int BLOCK_SIZE = 256;
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;

  check_build_bvh3d<exec, double>();
  check_build_bvh3d<exec, float>();
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(spin_bvh, find_bounding_boxes_3d_cuda)
{
  constexpr int BLOCK_SIZE = 256;
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;

  check_find_bounding_boxes3d<exec, double>();
  check_find_bounding_boxes3d<exec, float>();
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(spin_bvh, find_bounding_boxes_2d_cuda)
{
  constexpr int BLOCK_SIZE = 256;
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;

  check_find_bounding_boxes2d<exec, double>();
  check_find_bounding_boxes2d<exec, float>();
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(spin_bvh, find_rays_3d_cuda)
{
  constexpr int BLOCK_SIZE = 256;
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;

  check_find_rays3d<exec, double>();
  check_find_rays3d<exec, float>();
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(spin_bvh, find_rays_2d_cuda)
{
  constexpr int BLOCK_SIZE = 256;
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;

  check_find_rays2d<exec, double>();
  check_find_rays2d<exec, float>();
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(spin_bvh, find_points_3d_cuda)
{
  constexpr int BLOCK_SIZE = 256;
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;

  check_find_points3d<exec, double>();
  check_find_points3d<exec, float>();
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(spin_bvh, find_points_2d_cuda)
{
  constexpr int BLOCK_SIZE = 256;
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;

  check_find_points2d<exec, double>();
  check_find_points2d<exec, float>();
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(spin_bvh, single_box2d_cuda)
{
  constexpr int BLOCK_SIZE = 256;
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;

  check_single_box2d<exec, double>();
  check_single_box2d<exec, float>();
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(spin_bvh, single_box3d_cuda)
{
  constexpr int BLOCK_SIZE = 256;
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;

  check_single_box3d<exec, double>();
  check_single_box3d<exec, float>();
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(spin_bvh, use_pool_allocator)
{
  // Create a pool allocator on the device
  constexpr size_t POOL_SIZE = (1024 * 1024 * 1024) + 1;
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  umpire::Allocator allocator = rm.getAllocator(umpire::resource::Unified);

  umpire::Allocator pool_allocator =
    rm.makeAllocator<umpire::strategy::DynamicPool>("DEVICE_POOL",
                                                    allocator,
                                                    POOL_SIZE);

  const int allocID = pool_allocator.getId();

  // aliases & constants
  constexpr int NDIMS = 3;

  constexpr int BLOCK_SIZE = 256;
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;

  constexpr int NUM_BOXES = 1;

  using FloatType = double;

  // single bounding box in [0,1] x [0,1] x [0,1]
  FloatType* boxes = axom::allocate<FloatType>(6, allocID);
  axom::for_all<exec>(
    0,
    6,
    AXOM_LAMBDA(axom::IndexType idx) { boxes[idx] = (idx < 3) ? 0.0 : 1.0; });

  // construct a BVH with a single box
  spin::BVH<NDIMS, exec, FloatType> bvh(boxes, NUM_BOXES, allocID);
  EXPECT_EQ(bvh.getAllocatorID(), allocID);

  bvh.setScaleFactor(1.0);  // i.e., no scaling
  bvh.build();

  // run the find algorithm w/ the centroid of the bounding box as input.
  // Should return one and only one candidate that corresponds to the
  // single bounding box.
  FloatType* xc = axom::allocate<FloatType>(NUM_BOXES, allocID);
  FloatType* yc = axom::allocate<FloatType>(NUM_BOXES, allocID);
  FloatType* zc = axom::allocate<FloatType>(NUM_BOXES, allocID);
  axom::for_all<exec>(
    0,
    1,
    AXOM_LAMBDA(axom::IndexType idx) { xc[idx] = yc[idx] = zc[idx] = 0.5; });

  IndexType* offsets = axom::allocate<IndexType>(NUM_BOXES, allocID);
  IndexType* counts = axom::allocate<IndexType>(NUM_BOXES, allocID);
  IndexType* candidates = nullptr;
  bvh.findPoints(offsets, counts, candidates, NUM_BOXES, xc, yc, zc);
  EXPECT_TRUE(candidates != nullptr);

  // Ensure the BVH uses interally the supplied pool allocator
  EXPECT_EQ(rm.getAllocator(candidates).getId(), allocID);

  axom::deallocate(xc);
  axom::deallocate(yc);
  axom::deallocate(zc);
  axom::deallocate(boxes);
  axom::deallocate(offsets);
  axom::deallocate(counts);
  axom::deallocate(candidates);
}

#endif /* AXOM_USE_CUDA && AXOM_USE_RAJA && AXOM_USE_UMPIRE */

//------------------------------------------------------------------------------
#include "axom/slic/core/SimpleLogger.hpp"
using axom::slic::SimpleLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  SimpleLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}
