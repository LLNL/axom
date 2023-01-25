// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// axom includes
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"

#include "axom/spin/BVH.hpp"
#include "axom/spin/UniformGrid.hpp"

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
template <typename FloatType, int NDIMS>
void dump_ray(const std::string& file,
              FloatType t,
              const primal::Ray<FloatType, NDIMS>& ray)
{
  std::ofstream ofs(file.c_str());

  // Expand 2D points into 3D
  primal::Point<FloatType, 3> orig_3d(0.0), at_3d(0.0);
  primal::Vector<FloatType, 3> norm_3d(0.0);

  primal::Point<FloatType, NDIMS> ray_at_t = ray.at(t);
  for(int dim = 0; dim < NDIMS; dim++)
  {
    orig_3d[dim] = ray.origin()[dim];
    norm_3d[dim] = ray.direction()[dim];
    at_3d[dim] = ray_at_t[dim];
  }

  ofs << "# vtk DataFile Version 3.0\n";
  ofs << "Ray Data\n";
  ofs << "ASCII\n";
  ofs << "DATASET UNSTRUCTURED_GRID\n";

  ofs << "POINTS 2 double\n";
  ofs << orig_3d[0] << " " << orig_3d[1] << " " << orig_3d[2] << std::endl;
  ofs << at_3d[0] << " " << at_3d[1] << " " << at_3d[2] << std::endl;

  ofs << "CELLS 1 3\n";
  ofs << "2 0 1\n";

  ofs << "CELL_TYPES 1\n";
  ofs << "3\n";

  ofs << "CELL_DATA 1\n";
  ofs << "VECTORS normal double\n";
  ofs << norm_3d[0] << " " << norm_3d[1] << " " << norm_3d[2] << std::endl;

  ofs.close();
}

//------------------------------------------------------------------------------

/*!
 * \brief Given a mesh object, this method generates an array of axis-aligned
 *  bounding boxes corresponding to each constituent cell of the mesh and
 *  a corresponding centroid.
 *
 * \param [in]  mesh pointer to the mesh object.
 * \param [out] aabbs array of bounding boxes
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
template <typename FloatType, int NDIMS>
void generate_aabbs(const mint::Mesh* mesh,
                    primal::BoundingBox<FloatType, NDIMS>*& aabbs)
{
  // sanity checks
  EXPECT_TRUE(aabbs == nullptr);
  EXPECT_EQ(NDIMS, mesh->getDimension());

  using BoxType = typename primal::BoundingBox<FloatType, NDIMS>;

  // calculate some constants
  constexpr int nodes_per_dim = 1 << NDIMS;

  // allocate output arrays
  const IndexType ncells = mesh->getNumberOfCells();
  aabbs = axom::allocate<BoxType>(ncells);

  using exec_policy = axom::SEQ_EXEC;
  mint::for_all_cells<exec_policy, xargs::coords>(
    mesh,
    AXOM_LAMBDA(IndexType cellIdx,
                numerics::Matrix<double> & coords,
                const IndexType* AXOM_UNUSED_PARAM(nodeIds)) {
      primal::BoundingBox<double, NDIMS> range;

      for(IndexType inode = 0; inode < nodes_per_dim; ++inode)
      {
        primal::Point<double, NDIMS> node {coords.getColumn(inode)};

        range.addPoint(node);
      }  // END for all cells nodes

      aabbs[cellIdx].clear();
      aabbs[cellIdx].addBox(range);
    });

  // post-condition sanity checks
  EXPECT_TRUE(aabbs != nullptr);
}

//------------------------------------------------------------------------------

/*!
 * \brief Given a mesh object, this method generates an array of axis-aligned
 *  bounding boxes corresponding to each constituent cell of the mesh and
 *  a corresponding centroid.
 *
 * \param [in]  mesh pointer to the mesh object.
 * \param [out] aabbs array of bounding boxes
 * \param [out] c buffer to store the cell centroids
 *
 * \note The intent of this method is to generate synthetic input test data
 *  to test the functionality of the BVH.
 *
 * \warning This method allocates aabbs internally. The caller is responsible
 *  for properly deallocating aabbs.
 *
 * \pre aabbs == nullptr
 * \pre c != nullptr
 *
 * \post aabbs != nullptr
 */
template <typename FloatType, int NDIMS>
void generate_aabbs_and_centroids(const mint::Mesh* mesh,
                                  primal::BoundingBox<FloatType, NDIMS>*& aabbs,
                                  primal::Point<FloatType, NDIMS>* c)
{
  // sanity checks
  EXPECT_TRUE(aabbs == nullptr);
  EXPECT_TRUE(c != nullptr);

  // calculate some constants
  constexpr int NUM_NODES_PER_CELL = 1 << NDIMS;

  using BoxType = typename primal::BoundingBox<FloatType, NDIMS>;
  using PointType = typename primal::Point<FloatType, NDIMS>;

  EXPECT_EQ(NDIMS, mesh->getDimension());

  // allocate output arrays
  const IndexType ncells = mesh->getNumberOfCells();
  aabbs = axom::allocate<BoxType>(ncells);

  using exec_policy = axom::SEQ_EXEC;
  mint::for_all_cells<exec_policy, xargs::coords>(
    mesh,
    AXOM_LAMBDA(IndexType cellIdx,
                numerics::Matrix<double> & coords,
                const IndexType* AXOM_UNUSED_PARAM(nodeIds)) {
      BoxType range;

      PointType sum(0.0);

      for(IndexType inode = 0; inode < NUM_NODES_PER_CELL; ++inode)
      {
        const double* node = coords.getColumn(inode);

        PointType coords;
        for(int dim = 0; dim < NDIMS; dim++)
        {
          coords[dim] = node[dim];
          sum[dim] += node[dim];
        }

        range.addPoint(coords);
      }  // END for all cells nodes
      for(int dim = 0; dim < NDIMS; dim++)
      {
        sum[dim] /= NUM_NODES_PER_CELL;
      }

      c[cellIdx] = range.getCentroid();
#ifndef AXOM_DEVICE_CODE
      EXPECT_EQ(c[cellIdx], sum);
#endif

      aabbs[cellIdx] = range;
    });

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

  using BoxType = typename primal::BoundingBox<FloatType, NDIMS>;
  using PointType = typename primal::Point<FloatType, NDIMS>;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  BoxType* boxes = axom::allocate<BoxType>(2);
  boxes[0] = BoxType {PointType(0.), PointType(1.)};
  boxes[1] = BoxType {PointType(1.), PointType(2.)};

  spin::BVH<NDIMS, ExecSpace, FloatType> bvh;

  bvh.setScaleFactor(1.0);  // i.e., no scaling
  bvh.initialize(boxes, NUM_BOXES);

  int allocatorID = bvh.getAllocatorID();
  EXPECT_EQ(allocatorID, axom::execution_space<ExecSpace>::allocatorID());

  BoxType bounds = bvh.getBounds();

  for(int idim = 0; idim < NDIMS; ++idim)
  {
    EXPECT_DOUBLE_EQ(bounds.getMin()[idim], 0.0);
    EXPECT_DOUBLE_EQ(bounds.getMax()[idim], 2.0);
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

  using BoxType = typename primal::BoundingBox<FloatType, NDIMS>;
  using PointType = typename primal::Point<FloatType, NDIMS>;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  BoxType* boxes = axom::allocate<BoxType>(2);
  boxes[0] = BoxType {PointType(0.), PointType(1.)};
  boxes[1] = BoxType {PointType(1.), PointType(2.)};

  spin::BVH<NDIMS, ExecSpace, FloatType> bvh;

  bvh.setScaleFactor(1.0);  // i.e., no scaling
  bvh.initialize(boxes, NUM_BOXES);

  int allocatorID = bvh.getAllocatorID();
  EXPECT_EQ(allocatorID, axom::execution_space<ExecSpace>::allocatorID());

  BoxType bounds = bvh.getBounds();

  for(int idim = 0; idim < NDIMS; ++idim)
  {
    EXPECT_DOUBLE_EQ(bounds.getMin()[idim], 0.0);
    EXPECT_DOUBLE_EQ(bounds.getMax()[idim], 2.0);
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

  using BoxType = typename primal::BoundingBox<FloatType, NDIMS>;
  using PointType = typename primal::Point<FloatType, NDIMS>;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  // setup query bounding boxes: both boxes have lower left min at
  // (-1.0,-1.0,-1.0) but different upper right max.
  // The first bounding box is setup such that it intersects
  // 18 bounding boxes of the mesh.
  BoxType* query_boxes = axom::allocate<BoxType>(N);
  query_boxes[0].clear();
  query_boxes[1].clear();

  query_boxes[0].addPoint(PointType {-1., -1., -1.});
  query_boxes[1].addPoint(PointType {-1., -1., -1.});

  query_boxes[0].addPoint(PointType {2.5, 2.5, 1.5});
  query_boxes[1].addPoint(PointType {-0.5, -0.5, -0.5});

  // setup a test mesh (3 x 3 x 3)
  double lo[NDIMS] = {0.0, 0.0, 0.0};
  double hi[NDIMS] = {3.0, 3.0, 3.0};
  mint::UniformMesh mesh(lo, hi, 4, 4, 4);
  const IndexType ncells = mesh.getNumberOfCells();
  BoxType* aabbs = nullptr;
  generate_aabbs(&mesh, aabbs);
  EXPECT_TRUE(aabbs != nullptr);

  // construct the BVH
  spin::BVH<NDIMS, ExecSpace, FloatType> bvh;
  bvh.setScaleFactor(1.0);  // i.e., no scaling
  bvh.initialize(aabbs, ncells);

  // check BVH bounding box
  BoxType bounds = bvh.getBounds();

  for(int idim = 0; idim < NDIMS; ++idim)
  {
    EXPECT_DOUBLE_EQ(bounds.getMin()[idim], lo[idim]);
    EXPECT_DOUBLE_EQ(bounds.getMax()[idim], hi[idim]);
  }

  // traverse the BVH to find the candidates for all the bounding boxes
  axom::Array<IndexType> offsets(N);
  axom::Array<IndexType> counts(N);
  axom::Array<IndexType> candidates;
  bvh.findBoundingBoxes(offsets, counts, candidates, N, query_boxes);

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
  axom::deallocate(aabbs);

  axom::deallocate(query_boxes);

  axom::setDefaultAllocator(current_allocator);
}

//------------------------------------------------------------------------------

template <typename ExecSpace, typename FloatType>
void check_find_bounding_boxes2d()
{
  constexpr int NDIMS = 2;
  constexpr IndexType N = 2;

  using BoxType = typename primal::BoundingBox<FloatType, NDIMS>;
  using PointType = typename primal::Point<FloatType, NDIMS>;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  // setup query bounding boxes: both boxes are source at (-1.0,-1.0) but have
  // different max (upper right). The first box is setup to intersect six
  // bounding boxes of the mesh.
  BoxType* query_boxes = axom::allocate<BoxType>(N);
  query_boxes[0].clear();
  query_boxes[1].clear();

  query_boxes[0].addPoint(PointType {-1., -1.});
  query_boxes[1].addPoint(PointType {-1., -1.});

  query_boxes[0].addPoint(PointType {2.5, 1.5});
  query_boxes[1].addPoint(PointType {-0.1, -0.1});

  // setup a test mesh (3 x 3)
  double lo[NDIMS] = {0.0, 0.0};
  double hi[NDIMS] = {3.0, 3.0};
  mint::UniformMesh mesh(lo, hi, 4, 4);
  const IndexType ncells = mesh.getNumberOfCells();
  BoxType* aabbs = nullptr;
  generate_aabbs(&mesh, aabbs);
  EXPECT_TRUE(aabbs != nullptr);

  // construct the BVH
  spin::BVH<NDIMS, ExecSpace, FloatType> bvh;
  bvh.setScaleFactor(1.0);  // i.e., no scaling
  bvh.initialize(aabbs, ncells);

  // check BVH bounding box
  BoxType bounds = bvh.getBounds();

  for(int idim = 0; idim < NDIMS; ++idim)
  {
    EXPECT_DOUBLE_EQ(bounds.getMin()[idim], lo[idim]);
    EXPECT_DOUBLE_EQ(bounds.getMax()[idim], hi[idim]);
  }

  // traverse the BVH to find the candidates for all the bounding boxes
  axom::Array<IndexType> offsets(N);
  axom::Array<IndexType> counts(N);
  axom::Array<IndexType> candidates;
  bvh.findBoundingBoxes(offsets, counts, candidates, N, query_boxes);

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
  axom::deallocate(aabbs);

  axom::deallocate(query_boxes);

  axom::setDefaultAllocator(current_allocator);
}

//------------------------------------------------------------------------------
template <typename ExecSpace, typename FloatType>
void check_find_rays3d()
{
  constexpr int NDIMS = 3;
  constexpr IndexType N = 2;

  using BoxType = typename primal::BoundingBox<FloatType, NDIMS>;
  using RayType = typename primal::Ray<FloatType, NDIMS>;
  using PointType = typename primal::Point<FloatType, NDIMS>;
  using VectorType = typename primal::Vector<FloatType, NDIMS>;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  // setup query rays: both rays are source at (-1.0,-1.0) but point in
  // opposite directions. The first ray is setup such that it intersects
  // three bounding boxes of the mesh.
  RayType* query_rays = axom::allocate<RayType>(N);
  query_rays[0] =
    RayType {PointType {-1.0, -1.0, -1.0}, VectorType {1.0, 1.0, 1.0}};
  query_rays[1] =
    RayType {PointType {-1.0, -1.0, -1.0}, VectorType {-1.0, -1.0, -1.0}};

#ifdef VTK_DEBUG
  FloatType t = 5.0;
  dump_ray("ray0.vtk", t, query_rays[0]);
  dump_ray("ray1.vtk", t, query_rays[1]);
#endif

  // setup a test mesh
  double lo[NDIMS] = {0.0, 0.0, 0.0};
  double hi[NDIMS] = {3.0, 3.0, 3.0};
  mint::UniformMesh mesh(lo, hi, 4, 4, 4);
  const IndexType ncells = mesh.getNumberOfCells();
  BoxType* aabbs = nullptr;
  generate_aabbs(&mesh, aabbs);
  EXPECT_TRUE(aabbs != nullptr);

  // construct the BVH
  spin::BVH<NDIMS, ExecSpace, FloatType> bvh;
  bvh.setScaleFactor(1.0);  // i.e., no scaling
  bvh.initialize(aabbs, ncells);

  // check BVH bounding box
  BoxType bounds = bvh.getBounds();

  for(int idim = 0; idim < NDIMS; ++idim)
  {
    EXPECT_DOUBLE_EQ(bounds.getMin()[idim], lo[idim]);
    EXPECT_DOUBLE_EQ(bounds.getMax()[idim], hi[idim]);
  }

  // traverse the BVH to find the candidates for all the centroids
  axom::Array<IndexType> offsets(N);
  axom::Array<IndexType> counts(N);
  axom::Array<IndexType> candidates;
  bvh.findRays(offsets, counts, candidates, N, query_rays);

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
  axom::deallocate(aabbs);

  axom::deallocate(query_rays);

  axom::setDefaultAllocator(current_allocator);
}

//------------------------------------------------------------------------------

template <typename ExecSpace, typename FloatType>
void check_find_rays2d()
{
  constexpr int NDIMS = 2;
  constexpr IndexType N = 2;

  using BoxType = typename primal::BoundingBox<FloatType, NDIMS>;
  using RayType = typename primal::Ray<FloatType, NDIMS>;
  using PointType = typename primal::Point<FloatType, NDIMS>;
  using VectorType = typename primal::Vector<FloatType, NDIMS>;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  // setup query rays: both rays are source at (-1.0,-1.0) but point in
  // opposite directions. The first ray is setup such that it intersects
  // seven bounding boxes of the mesh.
  RayType* query_rays = axom::allocate<RayType>(N);
  query_rays[0] = RayType {PointType {-1.0, -1.0}, VectorType {1.0, 1.0}};
  query_rays[1] = RayType {PointType {-1.0, -1.0}, VectorType {-1.0, -1.0}};

#ifdef VTK_DEBUG
  FloatType t = 5.0;
  dump_ray("ray0.vtk",
           t,
           query_rays[0].origin()[0],
           query_rays[0].direction()[0],
           query_rays[0].origin()[1],
           query_rays[0].direction()[1]);
  dump_ray("ray1.vtk",
           t,
           query_rays[1].origin()[0],
           query_rays[1].direction()[0],
           query_rays[1].origin()[1],
           query_rays[1].direction()[1]);
#endif

  // setup a test mesh
  double lo[NDIMS] = {0.0, 0.0};
  double hi[NDIMS] = {3.0, 3.0};
  mint::UniformMesh mesh(lo, hi, 4, 4);
  const IndexType ncells = mesh.getNumberOfCells();
  BoxType* aabbs = nullptr;
  generate_aabbs(&mesh, aabbs);
  EXPECT_TRUE(aabbs != nullptr);

  // construct the BVH
  spin::BVH<NDIMS, ExecSpace, FloatType> bvh;
  bvh.setScaleFactor(1.0);  // i.e., no scaling
  bvh.initialize(aabbs, ncells);

  // check BVH bounding box
  BoxType bounds = bvh.getBounds();

  for(int idim = 0; idim < NDIMS; ++idim)
  {
    EXPECT_DOUBLE_EQ(bounds.getMin()[idim], lo[idim]);
    EXPECT_DOUBLE_EQ(bounds.getMax()[idim], hi[idim]);
  }

  // traverse the BVH to find the candidates for all the centroids
  axom::Array<IndexType> offsets(N);
  axom::Array<IndexType> counts(N);
  axom::Array<IndexType> candidates;
  bvh.findRays(offsets, counts, candidates, N, query_rays);

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
  axom::deallocate(aabbs);

  axom::deallocate(query_rays);

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

  using BoxType = typename primal::BoundingBox<FloatType, NDIMS>;
  using PointType = primal::Point<FloatType, NDIMS>;
  using PointDbl = primal::Point<double, NDIMS>;

  double lo[NDIMS] = {0.0, 0.0, 0.0};
  double hi[NDIMS] = {3.0, 3.0, 3.0};
  int res[NDIMS] = {N - 1, N - 1, N - 1};

  mint::UniformMesh mesh(lo, hi, N, N, N);
  FloatType* centroids_raw =
    mesh.createField<FloatType>("centroid", mint::CELL_CENTERED, NDIMS);
  PointType* centroids = reinterpret_cast<PointType*>(centroids_raw);
  const IndexType ncells = mesh.getNumberOfCells();

  BoxType* aabbs = nullptr;
  generate_aabbs_and_centroids(&mesh, aabbs, centroids);

  // construct the BVH
  spin::BVH<NDIMS, ExecSpace, FloatType> bvh;
  bvh.setScaleFactor(1.0);  // i.e., no scaling
  bvh.initialize(aabbs, ncells);

  BoxType bounds = bvh.getBounds();

  for(int idim = 0; idim < NDIMS; ++idim)
  {
    EXPECT_DOUBLE_EQ(bounds.getMin()[idim], lo[idim]);
    EXPECT_DOUBLE_EQ(bounds.getMax()[idim], hi[idim]);
  }

  // traverse the BVH to find the candidates for all the centroids
  axom::Array<IndexType> offsets(ncells);
  axom::Array<IndexType> counts(ncells);
  axom::Array<IndexType> candidates;
  bvh.findPoints(offsets, counts, candidates, ncells, centroids);

  spin::UniformGrid<IndexType, NDIMS> ug(lo, hi, res);

  for(IndexType i = 0; i < ncells; ++i)
  {
    PointDbl q = PointDbl {centroids[i][0], centroids[i][1], centroids[i][2]};
    const int donorCellIdx = ug.getBinIndex(q);
    EXPECT_EQ(counts[i], 1);
    EXPECT_EQ(donorCellIdx, candidates[offsets[i]]);
  }  // END for all cell centroids

  // check points that are outside by shifting the query points
  constexpr double OFFSET = 10.0;
  for(IndexType i = 0; i < ncells; ++i)
  {
    centroids[i][0] += OFFSET;
    centroids[i][1] += OFFSET;
    centroids[i][2] += OFFSET;
  }

  bvh.findPoints(offsets, counts, candidates, ncells, centroids);

  for(IndexType i = 0; i < ncells; ++i)
  {
    EXPECT_EQ(counts[i], 0);
  }

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

  using BoxType = typename primal::BoundingBox<FloatType, NDIMS>;
  using PointType = primal::Point<FloatType, NDIMS>;
  using PointDbl = primal::Point<double, NDIMS>;

  double lo[NDIMS] = {0.0, 0.0};
  double hi[NDIMS] = {3.0, 3.0};
  int res[NDIMS] = {N - 1, N - 1};

  mint::UniformMesh mesh(lo, hi, N, N);
  FloatType* centroids_raw =
    mesh.createField<FloatType>("centroid", mint::CELL_CENTERED, NDIMS);
  PointType* centroids = reinterpret_cast<PointType*>(centroids_raw);
  const IndexType ncells = mesh.getNumberOfCells();

  BoxType* aabbs = nullptr;
  generate_aabbs_and_centroids(&mesh, aabbs, centroids);

  // construct the BVH
  spin::BVH<NDIMS, ExecSpace, FloatType> bvh;
  bvh.setScaleFactor(1.0);  // i.e., no scaling
  bvh.initialize(aabbs, ncells);

  BoxType bounds = bvh.getBounds();

  for(int idim = 0; idim < NDIMS; ++idim)
  {
    EXPECT_DOUBLE_EQ(bounds.getMin()[idim], lo[idim]);
    EXPECT_DOUBLE_EQ(bounds.getMax()[idim], hi[idim]);
  }

  // traverse the BVH to find the candidates for all the centroids
  axom::Array<IndexType> offsets(ncells);
  axom::Array<IndexType> counts(ncells);
  axom::Array<IndexType> candidates;
  bvh.findPoints(offsets, counts, candidates, ncells, centroids);

  spin::UniformGrid<IndexType, NDIMS> ug(lo, hi, res);

  for(IndexType i = 0; i < ncells; ++i)
  {
    PointDbl q = PointDbl {centroids[i][0], centroids[i][1]};
    const int donorCellIdx = ug.getBinIndex(q);
    EXPECT_EQ(counts[i], 1);
    EXPECT_EQ(donorCellIdx, candidates[offsets[i]]);
  }  // END for all cell centroids

  // check points that are outside by shifting the query points
  constexpr double OFFSET = 10.0;
  for(IndexType i = 0; i < ncells; ++i)
  {
    centroids[i][0] += OFFSET;
    centroids[i][1] += OFFSET;
  }

  bvh.findPoints(offsets, counts, candidates, ncells, centroids);

  for(IndexType i = 0; i < ncells; ++i)
  {
    EXPECT_EQ(counts[i], 0);
  }

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

  using BoxType = typename primal::BoundingBox<FloatType, NDIMS>;
  using PointType = primal::Point<FloatType, NDIMS>;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  // single bounding box in [0,1] x [0,1]
  BoxType* boxes = axom::allocate<BoxType>(2);
  boxes[0] = BoxType {PointType(0.), PointType(1.)};

  // construct a BVH with a single box
  spin::BVH<NDIMS, ExecSpace, FloatType> bvh;
  bvh.setScaleFactor(1.0);  // i.e., no scaling
  bvh.initialize(boxes, NUM_BOXES);

  // check the bounds -- should match the bounds of the input bounding box
  BoxType bounds = bvh.getBounds();

  for(int idim = 0; idim < NDIMS; ++idim)
  {
    EXPECT_DOUBLE_EQ(bounds.getMin()[idim], 0.0);
    EXPECT_DOUBLE_EQ(bounds.getMax()[idim], 1.0);
  }

  // run the find algorithm w/ the centroid of the bounding box as input.
  // Should return one and only one candidate that corresponds to the
  // single bounding box.
  PointType centroid {0.5, 0.5};
  PointType* centroid_device = axom::allocate<PointType>(1);
  axom::copy(centroid_device, &centroid, sizeof(PointType));

  axom::Array<IndexType> offsets(NUM_BOXES);
  axom::Array<IndexType> counts(NUM_BOXES);
  axom::Array<IndexType> candidates;
  bvh.findPoints(offsets, counts, candidates, NUM_BOXES, centroid_device);
  EXPECT_EQ(counts[0], 1);
  EXPECT_EQ(0, candidates[offsets[0]]);

  // shift centroid outside of the BVH, should return no candidates.
  centroid[0] += 10.0;
  centroid[1] += 10.0;
  axom::copy(centroid_device, &centroid, sizeof(PointType));
  bvh.findPoints(offsets, counts, candidates, NUM_BOXES, centroid_device);

  EXPECT_EQ(counts[0], 0);

  axom::deallocate(centroid_device);
  axom::deallocate(boxes);

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

  using BoxType = typename primal::BoundingBox<FloatType, NDIMS>;
  using PointType = primal::Point<FloatType, NDIMS>;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  // single bounding box in [0,1] x [0,1] x [0,1]
  BoxType* boxes = axom::allocate<BoxType>(2);
  boxes[0] = BoxType {PointType(0.), PointType(1.)};

  // construct a BVH with a single box
  spin::BVH<NDIMS, ExecSpace, FloatType> bvh;
  bvh.setScaleFactor(1.0);  // i.e., no scaling
  bvh.initialize(boxes, NUM_BOXES);

  // check the bounds -- should match the bounds of the input bounding box
  BoxType bounds = bvh.getBounds();

  for(int idim = 0; idim < NDIMS; ++idim)
  {
    EXPECT_DOUBLE_EQ(bounds.getMin()[idim], 0.0);
    EXPECT_DOUBLE_EQ(bounds.getMax()[idim], 1.0);
  }

  // run the find algorithm w/ the centroid of the bounding box as input.
  // Should return one and only one candidate that corresponds to the
  // single bounding box.
  PointType centroid {0.5, 0.5, 0.5};
  PointType* centroid_device = axom::allocate<PointType>(1);
  axom::copy(centroid_device, &centroid, sizeof(PointType));

  axom::Array<IndexType> offsets(NUM_BOXES);
  axom::Array<IndexType> counts(NUM_BOXES);
  axom::Array<IndexType> candidates;
  bvh.findPoints(offsets, counts, candidates, NUM_BOXES, centroid_device);

  EXPECT_EQ(counts[0], 1);
  EXPECT_EQ(0, candidates[offsets[0]]);

  // shift centroid outside of the BVH, should return no candidates.
  centroid[0] += 10.0;
  centroid[1] += 10.0;
  axom::copy(centroid_device, &centroid, sizeof(PointType));
  bvh.findPoints(offsets, counts, candidates, NUM_BOXES, centroid_device);

  EXPECT_EQ(counts[0], 0);

  axom::deallocate(centroid_device);
  axom::deallocate(boxes);

  axom::setDefaultAllocator(current_allocator);
}

//------------------------------------------------------------------------------

/*!
 * \brief Tests the construction of the BVH in 3D by inserting two bounding
 *  boxes in the BVH and ensuring that the bounds of the BVH are as expected.
 */
template <typename ExecSpace, typename FloatType>
void check_build_bvh_zip3d()
{
  constexpr int NUM_BOXES = 2;
  constexpr int NDIMS = 3;

  using BoxType = typename primal::BoundingBox<FloatType, NDIMS>;
  using ZipIter = typename primal::ZipIndexable<BoxType>;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  FloatType* xmin = axom::allocate<FloatType>(NUM_BOXES);
  FloatType* ymin = axom::allocate<FloatType>(NUM_BOXES);
  FloatType* zmin = axom::allocate<FloatType>(NUM_BOXES);

  FloatType* xmax = axom::allocate<FloatType>(NUM_BOXES);
  FloatType* ymax = axom::allocate<FloatType>(NUM_BOXES);
  FloatType* zmax = axom::allocate<FloatType>(NUM_BOXES);

  xmin[0] = ymin[0] = zmin[0] = 0.;
  xmax[0] = ymax[0] = zmax[0] = 1.;
  xmin[1] = ymin[1] = zmin[1] = 1.;
  xmax[1] = ymax[1] = zmax[1] = 2.;

  ZipIter it {{xmin, ymin, zmin}, {xmax, ymax, zmax}};

  spin::BVH<NDIMS, ExecSpace, FloatType> bvh;

  bvh.setScaleFactor(1.0);  // i.e., no scaling
  bvh.initialize(it, NUM_BOXES);

  int allocatorID = bvh.getAllocatorID();
  EXPECT_EQ(allocatorID, axom::execution_space<ExecSpace>::allocatorID());

  BoxType bounds = bvh.getBounds();

  for(int idim = 0; idim < NDIMS; ++idim)
  {
    EXPECT_DOUBLE_EQ(bounds.getMin()[idim], 0.0);
    EXPECT_DOUBLE_EQ(bounds.getMax()[idim], 2.0);
  }

  axom::deallocate(xmin);
  axom::deallocate(ymin);
  axom::deallocate(zmin);

  axom::deallocate(xmax);
  axom::deallocate(ymax);
  axom::deallocate(zmax);

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
void check_find_points_zip3d()
{
  constexpr int NDIMS = 3;
  constexpr IndexType N = 4;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  using BoxType = typename primal::BoundingBox<FloatType, NDIMS>;
  using PointType = primal::Point<FloatType, NDIMS>;
  using PointDbl = primal::Point<double, NDIMS>;

  double lo[NDIMS] = {0.0, 0.0, 0.0};
  double hi[NDIMS] = {3.0, 3.0, 3.0};
  int res[NDIMS] = {N - 1, N - 1, N - 1};

  mint::UniformMesh mesh(lo, hi, N, N, N);
  FloatType* centroids_raw =
    mesh.createField<FloatType>("centroid", mint::CELL_CENTERED, NDIMS);
  PointType* centroids = reinterpret_cast<PointType*>(centroids_raw);
  const IndexType ncells = mesh.getNumberOfCells();

  BoxType* aabbs = nullptr;
  generate_aabbs_and_centroids(&mesh, aabbs, centroids);

  FloatType* xs = axom::allocate<FloatType>(ncells);
  FloatType* ys = axom::allocate<FloatType>(ncells);
  FloatType* zs = axom::allocate<FloatType>(ncells);

  for(int icell = 0; icell < ncells; icell++)
  {
    xs[icell] = centroids[icell][0];
    ys[icell] = centroids[icell][1];
    zs[icell] = centroids[icell][2];
  }

  primal::ZipIndexable<PointType> zip_test {{xs, ys, zs}};

  // construct the BVH
  spin::BVH<NDIMS, ExecSpace, FloatType> bvh;
  bvh.setScaleFactor(1.0);  // i.e., no scaling
  bvh.initialize(aabbs, ncells);

  BoxType bounds = bvh.getBounds();

  for(int idim = 0; idim < NDIMS; ++idim)
  {
    EXPECT_DOUBLE_EQ(bounds.getMin()[idim], lo[idim]);
    EXPECT_DOUBLE_EQ(bounds.getMax()[idim], hi[idim]);
  }

  // traverse the BVH to find the candidates for all the centroids
  axom::Array<IndexType> offsets(ncells);
  axom::Array<IndexType> counts(ncells);
  axom::Array<IndexType> candidates;
  bvh.findPoints(offsets, counts, candidates, ncells, zip_test);

  spin::UniformGrid<IndexType, NDIMS> ug(lo, hi, res);

  for(IndexType i = 0; i < ncells; ++i)
  {
    PointDbl q = PointDbl {centroids[i][0], centroids[i][1], centroids[i][2]};
    const int donorCellIdx = ug.getBinIndex(q);
    EXPECT_EQ(counts[i], 1);
    EXPECT_EQ(donorCellIdx, candidates[offsets[i]]);
  }  // END for all cell centroids

  axom::deallocate(xs);
  axom::deallocate(ys);
  axom::deallocate(zs);
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
void check_find_points_zip2d()
{
  constexpr int NDIMS = 2;
  constexpr IndexType N = 4;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  using BoxType = typename primal::BoundingBox<FloatType, NDIMS>;
  using PointType = primal::Point<FloatType, NDIMS>;
  using PointDbl = primal::Point<double, NDIMS>;

  double lo[NDIMS] = {0.0, 0.0};
  double hi[NDIMS] = {3.0, 3.0};
  int res[NDIMS] = {N - 1, N - 1};

  mint::UniformMesh mesh(lo, hi, N, N);
  FloatType* centroids_raw =
    mesh.createField<FloatType>("centroid", mint::CELL_CENTERED, NDIMS);
  PointType* centroids = reinterpret_cast<PointType*>(centroids_raw);
  const IndexType ncells = mesh.getNumberOfCells();

  BoxType* aabbs = nullptr;
  generate_aabbs_and_centroids(&mesh, aabbs, centroids);

  FloatType* xs = axom::allocate<FloatType>(ncells);
  FloatType* ys = axom::allocate<FloatType>(ncells);
  FloatType* zs = nullptr;

  for(int icell = 0; icell < ncells; icell++)
  {
    xs[icell] = centroids[icell][0];
    ys[icell] = centroids[icell][1];
  }

  primal::ZipIndexable<PointType> zip_test {{xs, ys, zs}};

  // construct the BVH
  spin::BVH<NDIMS, ExecSpace, FloatType> bvh;
  bvh.setScaleFactor(1.0);  // i.e., no scaling
  bvh.initialize(aabbs, ncells);

  BoxType bounds = bvh.getBounds();

  for(int idim = 0; idim < NDIMS; ++idim)
  {
    EXPECT_DOUBLE_EQ(bounds.getMin()[idim], lo[idim]);
    EXPECT_DOUBLE_EQ(bounds.getMax()[idim], hi[idim]);
  }

  // traverse the BVH to find the candidates for all the centroids
  axom::Array<IndexType> offsets(ncells);
  axom::Array<IndexType> counts(ncells);
  axom::Array<IndexType> candidates;
  bvh.findPoints(offsets, counts, candidates, ncells, zip_test);

  spin::UniformGrid<IndexType, NDIMS> ug(lo, hi, res);

  for(IndexType i = 0; i < ncells; ++i)
  {
    PointDbl q = PointDbl {centroids[i][0], centroids[i][1]};
    const int donorCellIdx = ug.getBinIndex(q);
    EXPECT_EQ(counts[i], 1);
    EXPECT_EQ(donorCellIdx, candidates[offsets[i]]);
  }  // END for all cell centroids

  axom::deallocate(aabbs);

  axom::setDefaultAllocator(current_allocator);
}

//---------------------------------------------------------------------------

/*!
 * \brief Makes an array of 2D points for testing the BVH.
 *
 * \return points The array of points to use for testing.
 */
template <typename FloatType>
axom::Array<axom::primal::Point<FloatType, 2>> make_query_points_2d()
{
  int dims[] = {101, 101};
  FloatType x0 = -1.5, x1 = 1.5;
  FloatType y0 = -1.5, y1 = 1.5;
  axom::Array<FloatType> x(dims[0]), y(dims[1]);
  axom::numerics::linspace(x0, x1, x.data(), dims[0]);
  axom::numerics::linspace(y0, y1, y.data(), dims[1]);
  axom::Array<axom::primal::Point<FloatType, 2>> points;
  points.reserve(dims[0] * dims[1]);
  for(int j = 0; j < dims[1]; j++)
  {
    for(int i = 0; i < dims[0]; i++)
    {
      points.push_back(axom::primal::Point<FloatType, 2>({x[i], y[j]}));
    }
  }
  return points;
}

//---------------------------------------------------------------------------

/*!
 * \brief Tests inserting a single point into the BVH and querying a set of
 *        points against it, using a pattern seen in DistributedClosestPoint.
 *        The test demonstrates that the BVH returns an in range value to
 *        the checkMinDist lambda.
 *
 * \param [in] bvh The BVH object to use.
 * \param [in] points The points that will be inserted in the BVH as bboxes.
 * \param [in] query_pts The array of query points to use against the BVH.
 * \param [in] none Whether to add no points to the BVH. If this is the case,
 *                  none of the queries against the BVH will set the minElem
 *                  value and we use none so we know the expected result.
 */
template <typename BVHType, typename PointType>
void bvh_compute_point_distances_2d(BVHType& bvh,
                                    axom::Array<PointType>& points,
                                    axom::Array<PointType>& query_pts,
                                    bool none)
{
  EXPECT_TRUE(!query_pts.empty());

  using ExecSpace = typename BVHType::ExecSpaceType;
  using FloatType = typename PointType::CoordType;
  using BoxType = axom::primal::BoundingBox<FloatType, 2>;
  constexpr int INVALID_ELEMENT_INDEX = -1;
  constexpr int FIRST_ELEMENT_INDEX = 0;

  // Borrowed from DistibutedClosestPoint.
  struct MinCandidate
  {
    /// Squared distance to query point
    double minSqDist {numerics::floating_point_limits<double>::max()};
    /// Index within mesh of closest element
    int minElem {INVALID_ELEMENT_INDEX};
  };

  // Initialize the BVH with the points.
  axom::IndexType npts = none ? 0 : points.size();
  BoxType* bboxes =
    axom::allocate<BoxType>(npts ? npts : 1,  // do not allocate 0 elements
                            axom::execution_space<ExecSpace>::allocatorID());
  for(axom::IndexType i = 0; i < npts; i++) bboxes[i] = BoxType(points[i]);
  bvh.initialize(bboxes, npts);

  // Call the BVH like the DistributedClosestPoint does.
  auto it = bvh.getTraverser();

  // Make a results array to contain the data values we read.
  axom::Array<int> results(query_pts.size());
  axom::ArrayView<int> results_view(results.data(), results.size());

  // Loop over the query points.
  axom::ArrayView<PointType> pointsView(points.data(), points.size());
  axom::ArrayView<PointType> query_pts_view(query_pts.data(), query_pts.size());
  npts = query_pts.size();
  axom::for_all<ExecSpace>(
    npts,
    AXOM_LAMBDA(int32 idx) mutable {
      // Get the current query point.
      auto qpt = query_pts_view[idx];
      MinCandidate curr_min;

      auto checkMinDist = [&](int32 current_node, const int32* leaf_nodes) {
        int candidate_idx = leaf_nodes[current_node];
        const PointType candidate_pt = pointsView[candidate_idx];
        const double sq_dist = squared_distance(qpt, candidate_pt);

        if(sq_dist < curr_min.minSqDist)
        {
          curr_min.minSqDist = sq_dist;
          curr_min.minElem = candidate_idx;
        }
      };

      // Borrowed from DistributedClosestPoint.
      auto traversePredicate = [&](const PointType& p, const BoxType& bb) -> bool {
        auto sqDist = squared_distance(p, bb);
        return sqDist <= curr_min.minSqDist;
      };

      // Traverse the tree, searching for the point with minimum distance.
      it.traverse_tree(qpt, checkMinDist, traversePredicate);

      // Save the index of the minElem.
      results_view[idx] = curr_min.minElem;
    });

  // Make sure all of the indices we found are the expected value.
  int expected_idx = none ? INVALID_ELEMENT_INDEX : FIRST_ELEMENT_INDEX;
  for(axom::IndexType i = 0; i < results.size(); i++)
  {
    EXPECT_EQ(expected_idx, results[i]);
  }

  axom::deallocate(bboxes);
}

//------------------------------------------------------------------------------
template <typename ExecType, typename FloatType>
void check_0_or_1_bbox_2d()
{
  // Temporarily force the allocator for the ExecType so the query points as
  // well as any other important allocations are available in that space.
  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecType>::allocatorID());

  // Make a 1 element array that we'll access from a BVH callback lambda.
  using PointType = axom::primal::Point<FloatType, 2>;
  axom::Array<PointType> src_pts;
  PointType a({0.45, 0.8});
  src_pts.push_back(a);

  // Make query points.
  auto query_pts = make_query_points_2d<FloatType>();

  // 1 src point
  spin::BVH<2, ExecType, FloatType> bvh;
  bvh_compute_point_distances_2d(bvh, src_pts, query_pts, false);

  // 0 src points
  spin::BVH<2, ExecType, FloatType> bvh2;
  bvh_compute_point_distances_2d(bvh2, src_pts, query_pts, true);

  axom::setDefaultAllocator(current_allocator);
}

} /* end unnamed namespace */

//------------------------------------------------------------------------------
// UNIT TESTS
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
TEST(spin_bvh, construct3D_sequential_zip)
{
  check_build_bvh_zip3d<axom::SEQ_EXEC, double>();
  check_build_bvh_zip3d<axom::SEQ_EXEC, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, find_points_2d_sequential_zip)
{
  check_find_points_zip2d<axom::SEQ_EXEC, double>();
  check_find_points_zip2d<axom::SEQ_EXEC, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, find_points_3d_sequential_zip)
{
  check_find_points_zip3d<axom::SEQ_EXEC, double>();
  check_find_points_zip3d<axom::SEQ_EXEC, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, single_bbox_sequential)
{
  check_0_or_1_bbox_2d<axom::SEQ_EXEC, double>();
  check_0_or_1_bbox_2d<axom::SEQ_EXEC, float>();
}

//------------------------------------------------------------------------------
#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)

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

//------------------------------------------------------------------------------
TEST(spin_bvh, single_bbox_omp)
{
  check_0_or_1_bbox_2d<axom::OMP_EXEC, double>();
  check_0_or_1_bbox_2d<axom::OMP_EXEC, float>();
}

#endif

//------------------------------------------------------------------------------
#if defined(AXOM_USE_GPU) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)

TEST(spin_bvh, construct2D_device)
{
  constexpr int BLOCK_SIZE = 256;

  #if defined(__CUDACC__)
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;
  #elif defined(__HIPCC__)
  using exec = axom::HIP_EXEC<BLOCK_SIZE>;
  #else
  using exec = axom::SEQ_EXEC;
  #endif

  check_build_bvh2d<exec, double>();
  check_build_bvh2d<exec, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, construct3D_device)
{
  constexpr int BLOCK_SIZE = 256;

  #if defined(__CUDACC__)
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;
  #elif defined(__HIPCC__)
  using exec = axom::HIP_EXEC<BLOCK_SIZE>;
  #else
  using exec = axom::SEQ_EXEC;
  #endif

  check_build_bvh3d<exec, double>();
  check_build_bvh3d<exec, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, find_bounding_boxes_3d_device)
{
  constexpr int BLOCK_SIZE = 256;

  #if defined(__CUDACC__)
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;
  #elif defined(__HIPCC__)
  using exec = axom::HIP_EXEC<BLOCK_SIZE>;
  #else
  using exec = axom::SEQ_EXEC;
  #endif

  check_find_bounding_boxes3d<exec, double>();
  check_find_bounding_boxes3d<exec, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, find_bounding_boxes_2d_device)
{
  constexpr int BLOCK_SIZE = 256;

  #if defined(__CUDACC__)
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;
  #elif defined(__HIPCC__)
  using exec = axom::HIP_EXEC<BLOCK_SIZE>;
  #else
  using exec = axom::SEQ_EXEC;
  #endif

  check_find_bounding_boxes2d<exec, double>();
  check_find_bounding_boxes2d<exec, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, find_rays_3d_device)
{
  constexpr int BLOCK_SIZE = 256;

  #if defined(__CUDACC__)
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;
  #elif defined(__HIPCC__)
  using exec = axom::HIP_EXEC<BLOCK_SIZE>;
  #else
  using exec = axom::SEQ_EXEC;
  #endif

  check_find_rays3d<exec, double>();
  check_find_rays3d<exec, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, find_rays_2d_device)
{
  constexpr int BLOCK_SIZE = 256;

  #if defined(__CUDACC__)
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;
  #elif defined(__HIPCC__)
  using exec = axom::HIP_EXEC<BLOCK_SIZE>;
  #else
  using exec = axom::SEQ_EXEC;
  #endif

  check_find_rays2d<exec, double>();
  check_find_rays2d<exec, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, find_points_3d_device)
{
  constexpr int BLOCK_SIZE = 256;

  #if defined(__CUDACC__)
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;
  #elif defined(__HIPCC__)
  using exec = axom::HIP_EXEC<BLOCK_SIZE>;
  #else
  using exec = axom::SEQ_EXEC;
  #endif

  check_find_points3d<exec, double>();
  check_find_points3d<exec, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, find_points_2d_device)
{
  constexpr int BLOCK_SIZE = 256;

  #if defined(__CUDACC__)
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;
  #elif defined(__HIPCC__)
  using exec = axom::HIP_EXEC<BLOCK_SIZE>;
  #else
  using exec = axom::SEQ_EXEC;
  #endif

  check_find_points2d<exec, double>();
  check_find_points2d<exec, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, single_box2d_device)
{
  constexpr int BLOCK_SIZE = 256;

  #if defined(__CUDACC__)
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;
  #elif defined(__HIPCC__)
  using exec = axom::HIP_EXEC<BLOCK_SIZE>;
  #else
  using exec = axom::SEQ_EXEC;
  #endif

  check_single_box2d<exec, double>();
  check_single_box2d<exec, float>();
}

//------------------------------------------------------------------------------
TEST(spin_bvh, single_box3d_device)
{
  constexpr int BLOCK_SIZE = 256;

  #if defined(__CUDACC__)
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;
  #elif defined(__HIPCC__)
  using exec = axom::HIP_EXEC<BLOCK_SIZE>;
  #else
  using exec = axom::SEQ_EXEC;
  #endif

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
    rm.makeAllocator<umpire::strategy::QuickPool>("DEVICE_POOL",
                                                  allocator,
                                                  POOL_SIZE);

  const int allocID = pool_allocator.getId();

  // aliases & constants
  constexpr int NDIMS = 3;

  constexpr int BLOCK_SIZE = 256;

  #if defined(__CUDACC__)
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;
  #elif defined(__HIPCC__)
  using exec = axom::HIP_EXEC<BLOCK_SIZE>;
  #else
  using exec = axom::SEQ_EXEC;
  #endif

  constexpr int NUM_BOXES = 1;

  using FloatType = double;
  using BoxType = typename primal::BoundingBox<FloatType, NDIMS>;
  using PointType = typename primal::Point<FloatType, NDIMS>;

  // single bounding box in [0,1] x [0,1] x [0,1]
  BoxType* boxes = axom::allocate<BoxType>(1, allocID);
  axom::for_all<exec>(
    0,
    1,
    AXOM_LAMBDA(axom::IndexType idx) {
      boxes[idx] = {PointType(0.), PointType(1.)};
    });

  // construct a BVH with a single box
  spin::BVH<NDIMS, exec, FloatType> bvh;

  bvh.setScaleFactor(1.0);  // i.e., no scaling
  bvh.setAllocatorID(allocID);
  bvh.initialize(boxes, NUM_BOXES);
  EXPECT_EQ(bvh.getAllocatorID(), allocID);

  // run the find algorithm w/ the centroid of the bounding box as input.
  // Should return one and only one candidate that corresponds to the
  // single bounding box.
  PointType* centroid = axom::allocate<PointType>(NUM_BOXES, allocID);
  axom::for_all<exec>(
    0,
    1,
    AXOM_LAMBDA(axom::IndexType idx) {
      centroid[idx] = PointType {0.5, 0.5, 0.5};
    });

  axom::Array<IndexType> offsets(NUM_BOXES, NUM_BOXES, allocID);
  axom::Array<IndexType> counts(NUM_BOXES, NUM_BOXES, allocID);
  axom::Array<IndexType> candidates;
  bvh.findPoints(offsets, counts, candidates, NUM_BOXES, centroid);
  EXPECT_TRUE(candidates.size() > 0);

  // Ensure the BVH uses interally the supplied pool allocator
  EXPECT_EQ(candidates.getAllocatorID(), allocID);

  axom::deallocate(centroid);
  axom::deallocate(boxes);
}

//---------------------------------------------------------------------------
TEST(spin_bvh, single_bbox_device)
{
  constexpr int BLOCK_SIZE = 256;

  #if defined(__CUDACC__)
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;
  #elif defined(__HIPCC__)
  using exec = axom::HIP_EXEC<BLOCK_SIZE>;
  #else
  using exec = axom::SEQ_EXEC;
  #endif

  check_0_or_1_bbox_2d<exec, double>();
  check_0_or_1_bbox_2d<exec, float>();
}

#endif /* AXOM_USE_GPU && AXOM_USE_RAJA && AXOM_USE_UMPIRE */

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
