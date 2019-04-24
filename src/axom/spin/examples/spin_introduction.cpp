// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! \file spin_introduction.cpp
 *  \brief This example code is a demonstration of Axom's spin component.
 *
 *  This file shows how to use spin with primal to instantiate, populate,
 *  and use a spatial index.  Running the executable produced from
 *  this file will produce a collection of Asymptote source files.  When
 *  compiled, the Asymptote files produce the figures that accompany spin's
 *  Sphinx documentation.
 */

/* This example code contains snippets used in the spin Sphinx documentation.
 * They begin and end with comments such as
 *
 * prims_header_start
 * prims_header_end
 * clip_header_start
 * clip_header_end
 * closest_point_header_start
 * closest_point_header_end
 *
 * each prepended with an underscore.
 */

// Axom primitives
#include "axom/spin/Primitives.hpp"
#include "axom/primal/geometry/Triangle.hpp"

// Axom operations
#include "axom/primal/operators/compute_bounding_box.hpp"
#include "axom/primal/operators/intersect.hpp"

// C++ headers
#include <cmath> // do we need this?
#include <iostream>
#include <fstream>

#include "fmt/fmt.hpp"

// almost all our examples are in 3D
const int in3D = 3;
const int in2D = 2;

// primitives represented by doubles in 3D
using BoundingBoxType = axom::spin::BoundingBox<double, in3D> ;
using PointType = axom::spin::Point<double, in3D> ;
using TriangleType = axom::primal::Triangle<double, in3D> ;

// Axom spatial index
// _ugrid_triintersect_header_start
#include "axom/spin/UniformGrid.hpp"
// the UniformGrid will store ints ("thing" indexes) in 3D
using UniformGridType = axom::spin::UniformGrid<int, in3D> ;
// _ugrid_triintersect_header_end
// _bvhtree_header_start
#include "axom/spin/BVHTree.hpp"
// the BVHTree is in 2D, storing an index to 2D triangles
using BVHTree2DType = axom::spin::BVHTree<int, in2D> ;
// supporting classes
using BoundingBox2DType = axom::spin::BoundingBox<double, in2D> ;
using Point2DType = axom::spin::Point<double, in2D> ;
using Triangle2DType = axom::primal::Triangle<double, in2D> ;
// _bvhtree_header_end

// _naive_triintersect_start
void findTriIntersectionsNaively(
  std::vector<TriangleType> & tris,
  std::vector< std::pair<int, int> > & clashes
  )
{
  int tcount = tris.size();

  for (int i = 0 ; i < tcount ; ++i)
  {
    TriangleType & t1 = tris[i];
    for (int j = i + 1 ; j < tcount ; ++j)
    {
      TriangleType & t2 = tris[j];
      if (intersect(t1, t2))
      {
        clashes.push_back(std::make_pair(i, j));
      }
    }
  }
}
// _naive_triintersect_end

BoundingBoxType findBbox(std::vector<TriangleType> & tris)
{
  BoundingBoxType bbox;

  for (size_t i = 0 ; i < tris.size() ; ++i)
  {
    bbox.addPoint(tris[i][0]);
    bbox.addPoint(tris[i][1]);
    bbox.addPoint(tris[i][2]);
  }

  return bbox;
}

BoundingBoxType findBbox(TriangleType & tri)
{
  BoundingBoxType bbox;

  bbox.addPoint(tri[0]);
  bbox.addPoint(tri[1]);
  bbox.addPoint(tri[2]);

  return bbox;
}

BoundingBox2DType findBbox(Triangle2DType & tri)
{
  BoundingBox2DType bbox;

  bbox.addPoint(tri[0]);
  bbox.addPoint(tri[1]);
  bbox.addPoint(tri[2]);

  return bbox;
}

// _ugrid_build_start
BoundingBoxType findBbox(std::vector<TriangleType> & tris);
BoundingBoxType findBbox(TriangleType & tri);

UniformGridType * buildUniformGrid(std::vector<TriangleType> & tris)
{
  // Prepare to construct the UniformGrid.
  BoundingBoxType allbbox = findBbox(tris);
  const PointType & minBBPt = allbbox.getMin();
  const PointType & maxBBPt = allbbox.getMax();

  int tcount = tris.size();

  // The number of bins along one side of the UniformGrid.
  // This is a heuristic.
  int res = (int)(1 + std::pow(tcount, 1/3.));
  int ress[3] = {res, res, res};

  // Construct the UniformGrid with minimum point, maximum point,
  // and number of bins along each side.  Then insert the triangles.
  UniformGridType * ugrid =
    new UniformGridType(minBBPt.data(), maxBBPt.data(), ress);
  for (int i = 0 ; i < tcount ; ++i)
  {
    TriangleType & t1 = tris[i];
    BoundingBoxType bbox = findBbox(t1);
    ugrid->insert(bbox, i);
  }

  return ugrid;
}
// _ugrid_build_end

// _ugrid_candidate_start
void findNeighborCandidates(TriangleType & t1,
                            int i,
                            UniformGridType * ugrid,
                            std::vector<int> & neighbors)
{
  BoundingBoxType bbox = findBbox(t1);

  // Get all the bins t1 occupies
  const std::vector<int> bToCheck = ugrid->getBinsForBbox(bbox);
  size_t checkcount = bToCheck.size();

  // Load all the triangles in these bins whose indices are
  // greater than i into a vector.
  for (size_t curb = 0 ; curb < checkcount ; ++curb)
  {
    std::vector<int> ntlist = ugrid->getBinContents(bToCheck[curb]);
    for (size_t j = 0 ; j < ntlist.size() ; ++j)
    {
      if (ntlist[j] > i)
      {
        neighbors.push_back(ntlist[j]);
      }
    }
  }

  // Sort the neighboring triangles, and throw out duplicates.
  // This is not strictly necessary but saves some calls to intersect().
  std::sort(neighbors.begin(), neighbors.end());
  std::vector<int>::iterator jend =
    std::unique(neighbors.begin(), neighbors.end());
  neighbors.erase(jend, neighbors.end());
}
// _ugrid_candidate_end

// _ugrid_triintersect_start
void findTriIntersectionsAccel(
  std::vector<TriangleType> & tris,
  UniformGridType * ugrid,
  std::vector< std::pair<int, int> > & clashes)
{
  int tcount = tris.size();

  // For each triangle t1,
  for (int i = 0 ; i < tcount ; ++i)
  {
    TriangleType & t1 = tris[i];
    std::vector<int> neighbors;
    findNeighborCandidates(t1, i, ugrid, neighbors);

    // Test for intersection between t1 and each of its neighbors.
    int ncount = neighbors.size();
    for (int n = 0; n < ncount; ++n)
    {
      int j = neighbors[n];
      TriangleType & t2 = tris[j];
      if (axom::primal::intersect(t1, t2))
      {
        clashes.push_back(std::make_pair(i, j));
      }
    }
  }
}
// _ugrid_triintersect_end

void makeUniformGridTriangles(std::vector<TriangleType> & tris)
{
  PointType p[8];
  p[0] = PointType::make_point(0.3,  0.93, 0.03);
  p[1] = PointType::make_point(0.1,  0.85, 0.01);
  p[2] = PointType::make_point(0.3,  0.78, 0.03);
  p[3] = PointType::make_point(0.18, 0.36, 0.018);
  p[4] = PointType::make_point(0.8,  0.58, 0.08);
  p[5] = PointType::make_point(0.6,  0.5, 0.06);
  p[6] = PointType::make_point(0.55, 0.42, 0.055);
  p[7] = PointType::make_point(0.61, 0.1, 0.061);

  TriangleType t0(p[0], p[1], p[2]);    tris.push_back(t0);
  TriangleType t1(p[2], p[1], p[3]);    tris.push_back(t1);
  TriangleType t2(p[2], p[3], p[6]);    tris.push_back(t2);
  TriangleType t3(p[6], p[3], p[7]);    tris.push_back(t3);
  TriangleType t4(p[4], p[2], p[6]);    tris.push_back(t4);
  TriangleType t5(p[4], p[5], p[7]);    tris.push_back(t5);
}

void printPairs(std::string title,
                std::vector< std::pair<int, int> > & clashes)
{
  int ccount = clashes.size();
  std::cout << ccount << title << std::endl;
  for (int i = 0 ; i < ccount ; ++i)
  {
    std::cout << clashes[i].first << "   " << clashes[i].second << std::endl;
  }
}

void driveUniformGrid()
{
  std::vector<TriangleType> tris;
  makeUniformGridTriangles(tris);

  std::vector< std::pair<int, int> > naiveclashes;
  findTriIntersectionsNaively(tris, naiveclashes);

  UniformGridType * ugrid = buildUniformGrid(tris);
  std::vector< std::pair<int, int> > accelclashes;
  findTriIntersectionsAccel(tris, ugrid, accelclashes);
  delete ugrid;

  std::cout << "----- driveUniformGrid -----" << std::endl;
  printPairs(" clashes found by naive algorithm:", naiveclashes);
  printPairs(" clashes found by accelerated algorithm:", accelclashes);
  std::cout << std::endl;
}

// _bvhtree_build_start
BoundingBox2DType findBbox(Triangle2DType & tri);

BVHTree2DType* buildBVHTree(std::vector<Triangle2DType> & tris)
{
  // Initialize BVHTree with the triangles
  const int MaxBinFill = 1;
  const int MaxLevels = 4;
  int tricount = tris.size();
  BVHTree2DType* tree  = new BVHTree2DType( tricount, MaxLevels );

  for ( int i=0 ; i < tricount ; ++i )
  {
    tree->insert( findBbox( tris[i] ), i );
  }

  // Build bounding volume hierarchy
  tree->build( MaxBinFill );

  return tree;
}
// _bvhtree_build_end

// _bvhtree_candidate_start
void findCandidateBVHTreeBins(BVHTree2DType * tree,
                              Point2DType ppoint,
                              std::vector<int> & candidates)
{
  // Which triangles does the probe point intersect?
  // Get the candidate bins
  std::vector<int> bins;
  tree->find(ppoint, bins);
  size_t nbins = bins.size();

  // for each candidate bin,
  for (size_t curb = 0 ; curb < nbins ; ++curb)
  {
    // get its size and object array
    int bcount = tree->getBucketNumObjects(bins[curb]);
    const int* ary = tree->getBucketObjectArray(bins[curb]);

    // For each object in the current bin,
    for (int j = 0 ; j < bcount ; ++j)
    {
      // find the tree's internal object ID
      int treeObjID = ary[j];
      // and use it to retrieve the triangle's ID.
      int triID = tree->getObjectData(treeObjID);

      // Then store the ID in the candidates list.
      candidates.push_back(triID);
    }
  }

  // Sort the candidate triangles, and throw out duplicates.
  // This is not strictly necessary but saves some calls to checkInTriangle().
  std::sort(candidates.begin(), candidates.end());
  std::vector<int>::iterator jend =
    std::unique(candidates.begin(), candidates.end());
  candidates.erase(jend, candidates.end());
}
// _bvhtree_candidate_end

// _bvhtree_cand_int_start
void findIntersectionsWithCandidates(std::vector<Triangle2DType> & tris,
                                     std::vector<int> & candidates,
                                     Point2DType ppoint,
                                     std::vector<int> & intersections)
{
  // Test if ppoint lands in any of its neighbor triangles.
  int csize = candidates.size();
  for (int i = 0; i < csize; ++i)
  {
    Triangle2DType & t = tris[candidates[i]];
    if (t.checkInTriangle(ppoint))
    {
      intersections.push_back(candidates[i]);
    }
  }
}
// _bvhtree_cand_int_end

void makeBVHTreeTriangles(std::vector<Triangle2DType> & tris)
{
  Point2DType p[19];
  p[0] = Point2DType::make_point(.13, .88);
  p[1] = Point2DType::make_point(.26, .87);
  p[2] = Point2DType::make_point(.11, .77);
  p[3] = Point2DType::make_point(.18, .78);
  p[4] = Point2DType::make_point(.13, .74);
  p[5] = Point2DType::make_point(.37, .75);
  p[6] = Point2DType::make_point(.12, .61);
  p[7] = Point2DType::make_point(.25, .51);
  p[8] = Point2DType::make_point(.11, .44);
  p[9] = Point2DType::make_point(.26, .40);
  p[10] = Point2DType::make_point(.12, .25);
  p[11] = Point2DType::make_point(.85, .38);
  p[12] = Point2DType::make_point(.94, .37);
  p[13] = Point2DType::make_point(.84, .26);
  p[14] = Point2DType::make_point(.92, .27);
  p[15] = Point2DType::make_point(.96, .28);
  p[16] = Point2DType::make_point(.84, .16);
  p[17] = Point2DType::make_point(.92, .16);
  p[18] = Point2DType::make_point(.93, .09);

  Triangle2DType t0(p[1], p[0], p[3]);      tris.push_back(t0);
  Triangle2DType t1(p[0], p[2], p[3]);      tris.push_back(t1);
  Triangle2DType t2(p[1], p[3], p[5]);      tris.push_back(t2);
  Triangle2DType t3(p[3], p[2], p[4]);      tris.push_back(t3);
  Triangle2DType t4(p[3], p[4], p[5]);      tris.push_back(t4);
  Triangle2DType t5(p[2], p[6], p[4]);      tris.push_back(t5);
  Triangle2DType t6(p[5], p[4], p[6]);      tris.push_back(t6);
  Triangle2DType t7(p[5], p[6], p[7]);      tris.push_back(t7);
  Triangle2DType t8(p[7], p[6], p[8]);      tris.push_back(t8);
  Triangle2DType t9(p[7], p[8], p[9]);      tris.push_back(t9);
  Triangle2DType t10(p[9], p[8], p[10]);    tris.push_back(t10);
  Triangle2DType t11(p[11], p[13], p[14]);  tris.push_back(t11);
  Triangle2DType t12(p[12], p[11], p[14]);  tris.push_back(t12);
  Triangle2DType t13(p[12], p[14], p[15]);  tris.push_back(t13);
  Triangle2DType t14(p[14], p[13], p[16]);  tris.push_back(t14);
  Triangle2DType t15(p[14], p[16], p[17]);  tris.push_back(t15);
  Triangle2DType t16(p[16], p[18], p[17]);  tris.push_back(t16);
}

void driveBVHTree()
{
  std::vector<Triangle2DType> tris;
  makeBVHTreeTriangles(tris);
  Point2DType ppoint = Point2DType::make_point(0.45, 0.25);
  std::vector<int> intersections, candidates;

  BVHTree2DType* tree = buildBVHTree(tris);
  findCandidateBVHTreeBins(tree, ppoint, candidates);
  findIntersectionsWithCandidates(tris, candidates, ppoint, intersections);
  tree->writeVtkFile("BVHTree.out.vtk");

  std::cout << "----- driveBVHTree -----" << std::endl;
  std::cout << "Point " << ppoint << " hit the following " <<
    intersections.size() << " triangles (0 expected):" << std::endl;
  for (size_t i = 0 ; i < intersections.size() ; ++i)
  {
    std::cout << intersections[i] << std::endl;
  }

  delete tree;
}

int main(int argc, char** argv)
{

  // Deal with unused variables
  AXOM_DEBUG_VAR(argc);
  AXOM_DEBUG_VAR(argv);

  driveUniformGrid();
  driveBVHTree();

  return 0;
}
