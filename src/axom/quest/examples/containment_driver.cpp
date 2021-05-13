// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file containment_driver.cpp
 * \brief Basic demo of point containment acceleration structure over surfaces.
 */

// Axom includes
#include "axom/core.hpp"
#include "axom/primal.hpp"
#include "axom/quest.hpp"
#include "axom/spin.hpp"
#include "axom/mint.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"

#include "fmt/fmt.hpp"
#include "CLI11/CLI11.hpp"

// C/C++ includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <fstream>
#include <iomanip>  // for setprecision()

namespace mint = axom::mint;
namespace primal = axom::primal;
namespace quest = axom::quest;
namespace slic = axom::slic;

using UMesh = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;

using TriVertIndices = primal::Point<axom::IndexType, 3>;
using SpaceTriangle = primal::Triangle<double, 3>;

using Octree3D = quest::InOutOctree<3>;

using GeometricBoundingBox = Octree3D::GeometricBoundingBox;
using SpacePt = Octree3D::SpacePt;
using SpaceVector = Octree3D::SpaceVector;
using GridPt = Octree3D::GridPt;
using BlockIndex = Octree3D::BlockIndex;

//------------------------------------------------------------------------------

/** Computes the bounding box of the surface mesh */
GeometricBoundingBox compute_bounds(mint::Mesh* mesh)
{
  SLIC_ASSERT(mesh != nullptr);

  GeometricBoundingBox meshBB;
  SpacePt pt;

  for(int i = 0; i < mesh->getNumberOfNodes(); ++i)
  {
    mesh->getNode(i, pt.data());
    meshBB.addPoint(pt);
  }

  SLIC_ASSERT(meshBB.isValid());

  return meshBB;
}

/**
 * Query the inOutOctree using uniform grid of resolution \a gridRes
 * in region defined by bounding box \a queryBounds
 */
void testContainmentOnRegularGrid(const Octree3D& inOutOctree,
                                  const GeometricBoundingBox& queryBounds,
                                  const GridPt& gridRes,
                                  bool isBatched,
                                  bool shouldOutputMeshes)
{
  const double* low = queryBounds.getMin().data();
  const double* high = queryBounds.getMax().data();
  mint::UniformMesh* umesh =
    new mint::UniformMesh(low, high, gridRes[0], gridRes[1], gridRes[2]);

  const int nnodes = umesh->getNumberOfNodes();
  int* containment = umesh->createField<int>("containment", mint::NODE_CENTERED);

  SLIC_ASSERT(containment != nullptr);

  axom::utilities::Timer timer;
  if(!isBatched)
  {
    timer.start();
    for(int inode = 0; inode < nnodes; ++inode)
    {
      primal::Point<double, 3> pt;
      umesh->getNode(inode, pt.data());

      containment[inode] = inOutOctree.within(pt) ? 1 : 0;
    }
    timer.stop();
  }
  else
  {
    timer.start();

    // Allocate space for the coordinate arrays
    double* x = axom::allocate<double>(nnodes);
    double* y = axom::allocate<double>(nnodes);
    double* z = axom::allocate<double>(nnodes);

// Determine an appropriate execution policy
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)
    using ExecPolicy = axom::OMP_EXEC;
#else
    using ExecPolicy = axom::SEQ_EXEC;
#endif

    // Fill the coordinate arrays
    mint::for_all_nodes<ExecPolicy, mint::xargs::xyz>(
      umesh,
      AXOM_LAMBDA(axom::IndexType idx, double xx, double yy, double zz) {
        x[idx] = xx;
        y[idx] = yy;
        z[idx] = zz;
      });

    // Loop through the points using ExecPolicy
    axom::for_all<ExecPolicy>(0, nnodes, [&](axom::IndexType idx) {
      const bool inside = inOutOctree.within(SpacePt {x[idx], y[idx], z[idx]});
      containment[idx] = inside ? 1 : 0;
    });

    // Deallocate the coordinate arrays
    axom::deallocate(x);
    axom::deallocate(y);
    axom::deallocate(z);

    timer.stop();
  }

  SLIC_INFO(
    fmt::format("\tQuerying {}^3 containment field "
                "took {} seconds (@ {} queries per second)",
                gridRes,
                timer.elapsed(),
                nnodes / timer.elapsed()));

  if(shouldOutputMeshes)
  {
    std::stringstream sstr;
    const auto& res = gridRes;

    const bool resAllSame = (res[0] == res[1] && res[1] == res[2]);
    std::string resStr = resAllSame
      ? fmt::format("{}", res[0])
      : fmt::format("{}_{}_{}", res[0], res[1], res[2]);

    sstr << "gridContainment_" << resStr << ".vtk";
    mint::write_vtk(umesh, sstr.str());
  }

  delete umesh;
}

/**
 * \brief Extracts the vertex indices of cell cellIndex from the mesh
 */
TriVertIndices getTriangleVertIndices(mint::Mesh* mesh, axom::IndexType cellIndex)
{
  SLIC_ASSERT(mesh != nullptr);
  SLIC_ASSERT(cellIndex >= 0 && cellIndex < mesh->getNumberOfCells());

  TriVertIndices tvInd;
  mesh->getCellNodeIDs(cellIndex, tvInd.data());
  return tvInd;
}

/**
 * \brief Extracts the positions of a traingle's vertices from the mesh
 * \return The triangle vertex positions in a SpaceTriangle instance
 */
SpaceTriangle getMeshTriangle(mint::Mesh* mesh, const TriVertIndices& vertIndices)
{
  SLIC_ASSERT(mesh != nullptr);

  SpaceTriangle tri;
  for(int i = 0; i < 3; ++i) mesh->getNode(vertIndices[i], tri[i].data());

  return tri;
}

/**
 * \brief Computes some statistics about the surface mesh.
 *
 * Specifically, computes histograms (and ranges) of the edge lengths and
 * triangle areas on a logarithmic scale and logs the results
 */
void print_surface_stats(mint::Mesh* mesh)
{
  SLIC_ASSERT(mesh != nullptr);

  SpacePt pt;

  using MinMaxRange = primal::BoundingBox<double, 1>;
  using LengthType = MinMaxRange::PointType;

  MinMaxRange meshEdgeLenRange;
  MinMaxRange meshTriAreaRange;
  const int nCells = mesh->getNumberOfCells();
  using TriIdxSet = std::set<axom::IndexType>;
  TriIdxSet badTriangles;

  // simple binning based on the exponent
  using LogHistogram = std::map<int, int>;
  LogHistogram edgeLenHist;  // Create histogram of edge lengths (log scale)
  LogHistogram areaHist;     // Create histogram of triangle areas (log scale)

  using LogRangeMap = std::map<int, MinMaxRange>;
  LogRangeMap edgeLenRangeMap;  // Tracks range of edge lengths at each scale
  LogRangeMap areaRangeMap;     // Tracks range of triangle areas at each scale

  using TriVertIndices = primal::Point<axom::IndexType, 3>;
  int expBase2;

  // Traverse mesh triangles and bin the edge lengths and areas
  for(int i = 0; i < nCells; ++i)
  {
    // Get the indices and positions of the triangle's three vertices
    TriVertIndices vertIndices = getTriangleVertIndices(mesh, i);
    SpaceTriangle tri = getMeshTriangle(mesh, vertIndices);

    // Compute edge stats -- note edges are double counted
    for(int j = 0; j < 3; ++j)
    {
      double len = SpaceVector(tri[j], tri[(j + 1) % 3]).norm();
      if(axom::utilities::isNearlyEqual(len, 0.))
      {
        badTriangles.insert(i);
      }
      else
      {
        LengthType edgeLen(len);
        meshEdgeLenRange.addPoint(edgeLen);
        std::frexp(len, &expBase2);
        edgeLenHist[expBase2]++;
        edgeLenRangeMap[expBase2].addPoint(edgeLen);
      }
    }

    // Compute triangle area stats
    double area = tri.area();
    if(axom::utilities::isNearlyEqual(area, 0.))
    {
      badTriangles.insert(i);
    }
    else
    {
      LengthType triArea(area);
      meshTriAreaRange.addPoint(triArea);
      std::frexp(area, &expBase2);
      areaHist[expBase2]++;
      areaRangeMap[expBase2].addPoint(triArea);
    }
  }

  // Log the results
  const int nVerts = mesh->getNumberOfNodes();
  SLIC_INFO(
    fmt::format("Mesh has {} vertices  and {} triangles.", nVerts, nCells));

  SLIC_INFO("Edge length range: " << meshEdgeLenRange);
  SLIC_INFO("Triangle area range is: " << meshTriAreaRange);

  fmt::memory_buffer edgeHistStr;
  fmt::format_to(edgeHistStr, "Edge length histogram (lg-arithmic): ");
  for(LogHistogram::const_iterator it = edgeLenHist.begin();
      it != edgeLenHist.end();
      ++it)
  {
    fmt::format_to(edgeHistStr,
                   "\n\texp: {}\tcount: {}\tRange: {}",
                   it->first,
                   it->second / 2,
                   edgeLenRangeMap[it->first]);
  }
  SLIC_DEBUG(fmt::to_string(edgeHistStr));

  fmt::memory_buffer triHistStr;
  fmt::format_to(triHistStr, "Triangle areas histogram (lg-arithmic): ");
  for(LogHistogram::const_iterator it = areaHist.begin(); it != areaHist.end();
      ++it)
  {
    fmt::format_to(triHistStr,
                   "\n\texp: {}\tcount: {}\tRange: {}",
                   it->first,
                   it->second,
                   areaRangeMap[it->first]);
  }
  SLIC_DEBUG(fmt::to_string(triHistStr));

  if(!badTriangles.empty())
  {
    fmt::memory_buffer badTriStr;
    fmt::format_to(badTriStr,
                   "The following triangle(s) have zero area/edge lengths:");
    for(TriIdxSet::const_iterator it = badTriangles.begin();
        it != badTriangles.end();
        ++it)
    {
      fmt::format_to(badTriStr, "\n\tTriangle {}", *it);
      TriVertIndices vertIndices;
      mesh->getCellNodeIDs(*it, vertIndices.data());

      SpacePt vertPos;
      for(int j = 0; j < 3; ++j)
      {
        mesh->getNode(vertIndices[j], vertPos.data());
        fmt::format_to(badTriStr,
                       "\n\t\t vId: {} @ position: {}",
                       vertIndices[j],
                       vertPos);
      }
    }
    SLIC_DEBUG(fmt::to_string(badTriStr));
  }
}

/**
 * \brief Finds the octree leaf containing the given query point,
 * and optionally refines the leaf
 */
void refineAndPrint(Octree3D& octree,
                    const SpacePt& queryPt,
                    bool shouldRefine = true)
{
  BlockIndex leafBlock = octree.findLeafBlock(queryPt);

  if(shouldRefine)
  {
    octree.refineLeaf(leafBlock);
    leafBlock = octree.findLeafBlock(queryPt);
  }

  GeometricBoundingBox blockBB = octree.blockBoundingBox(leafBlock);
  bool containsPt = blockBB.contains(queryPt);

  SLIC_INFO(
    fmt::format("\t(gridPt: {}; lev: {}) with bounds {} {} query point.",
                leafBlock.pt(),
                leafBlock.level(),
                blockBB,
                (containsPt ? " contains " : "does not contain ")));
}

/** Struct to parse and store the input parameters */
struct Input
{
public:
  std::string stlFile;
  int maxQueryLevel {7};
  std::vector<double> queryBoxMins;
  std::vector<double> queryBoxMaxs;

private:
  bool m_verboseOutput {false};
  bool m_hasUserQueryBox {false};
  bool m_use_batched_query {false};

public:
  Input()
  {
// set default stl file, when the 'data' submodule is available
#ifdef AXOM_DATA_DIR
    {
      namespace fs = axom::utilities::filesystem;
      const auto dir = fs::joinPath(AXOM_DATA_DIR, "quest");
      stlFile = fs::joinPath(dir, "plane_simp.stl");
    }
#endif

// increase default for max query resolution in release builds
#ifndef AXOM_DEBUG
    maxQueryLevel += 2;
#endif
  }

  bool isVerbose() const { return m_verboseOutput; }

  bool useBatchedQuery() const { return m_use_batched_query; }

  bool hasUserQueryBox() const { return m_hasUserQueryBox; }

  void parse(int argc, char** argv, CLI::App& app)
  {
    app.add_option("stlFile", stlFile, "Path to input mesh")->check(CLI::ExistingFile);

    app
      .add_flag("-v,--verbose",
                m_verboseOutput,
                "Enable/disable verbose output, "
                "including outputting generated containment grids.")
      ->capture_default_str();

    app
      .add_option(
        "-l,--levels",
        maxQueryLevel,
        "Max query resolution. \n"
        "Will query uniform grids at levels 1 through the provided level")
      ->capture_default_str()
      ->check(CLI::PositiveNumber);

    // Optional bounding box for query region
    auto* minbb =
      app.add_option("--min", queryBoxMins, "Min bounds for query box (x,y,z)")
        ->expected(3);
    auto* maxbb =
      app.add_option("--max", queryBoxMaxs, "Max bounds for query box (x,y,z)")
        ->expected(3);
    minbb->needs(maxbb);
    maxbb->needs(minbb);

    app
      .add_flag("--batched",
                m_use_batched_query,
                "uses a single batched query on all points instead of many "
                "individual queries")
      ->capture_default_str();

    app.get_formatter()->column_width(35);

    // could throw an exception
    app.parse(argc, argv);

    slic::setLoggingMsgLevel(m_verboseOutput ? slic::message::Debug
                                             : slic::message::Info);

    m_hasUserQueryBox = app.count("--min") == 3 && app.count("--max") == 3;
  }
};

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  axom::slic::SimpleLogger logger;  // create & initialize logger
  // slic::debug::checksAreErrors = true;

  // Set up and parse command line arguments
  Input params;
  CLI::App app {"Driver for In/Out surface containment query"};

  try
  {
    params.parse(argc, argv, app);
  }
  catch(const CLI::ParseError& e)
  {
    return app.exit(e);
  }

  // Load mesh file
  mint::Mesh* surface_mesh = nullptr;
  {
    SLIC_INFO(fmt::format("{:*^80}", " Loading the mesh "));
    SLIC_INFO("Reading file: " << params.stlFile << "...");

    quest::STLReader* reader = new quest::STLReader();
    reader->setFileName(params.stlFile);
    reader->read();

    // Create surface mesh
    surface_mesh = new UMesh(3, mint::TRIANGLE);
    reader->getMesh(static_cast<UMesh*>(surface_mesh));

    SLIC_INFO("Mesh has " << surface_mesh->getNumberOfNodes() << " nodes and "
                          << surface_mesh->getNumberOfCells() << " cells.");

    delete reader;
    reader = nullptr;
  }
  SLIC_ASSERT(surface_mesh != nullptr);

  // Compute mesh bounding box and log some stats about the surface
  GeometricBoundingBox meshBB = compute_bounds(surface_mesh);
  SLIC_INFO("Mesh bounding box: " << meshBB);
  print_surface_stats(surface_mesh);

  // Create octree over mesh's bounding box
  SLIC_INFO(fmt::format("{:*^80}", " Generating the octree "));
  Octree3D octree(meshBB, surface_mesh);
  octree.generateIndex();

  print_surface_stats(surface_mesh);
  mint::write_vtk(surface_mesh, "meldedTriMesh.vtk");

  SLIC_INFO(fmt::format("{:*^80}", " Querying the octree "));

  // Define the query bounding box
  GeometricBoundingBox queryBB;
  if(params.hasUserQueryBox())
  {
    queryBB.addPoint(SpacePt(params.queryBoxMins.data(), 3));
    queryBB.addPoint(SpacePt(params.queryBoxMaxs.data(), 3));
  }
  else
  {
    queryBB = octree.boundingBox();
  }
  SLIC_INFO("Bounding box for query points: " << queryBB);

  // Query the mesh
  for(int i = 1; i < params.maxQueryLevel; ++i)
  {
    int res = 1 << i;
    testContainmentOnRegularGrid(octree,
                                 queryBB,
                                 GridPt::make_point(res, res, res),
                                 params.useBatchedQuery(),
                                 params.isVerbose());
  }

  if(!params.isVerbose())
  {
    slic::setLoggingMsgLevel(slic::message::Warning);
  }

  // Other octree operations
  //-- find leaf block of a given query point at various levels of resolution
  SLIC_INFO(fmt::format("{:*^80}", " Other octree operations "));
  double alpha = 2. / 3.;
  SpacePt queryPt = SpacePt::lerp(meshBB.getMin(), meshBB.getMax(), alpha);

  SLIC_INFO("Finding associated grid point for query point: " << queryPt);
  for(int lev = 0; lev < octree.maxLeafLevel(); ++lev)
  {
    GridPt gridPt = octree.findGridCellAtLevel(queryPt, lev);
    SLIC_INFO(
      fmt::format("  {1} @ level {0}\n\t[max gridPt: {2}; spacing: {3};\n\t "
                  "bounding box {4}]",
                  lev,
                  gridPt,
                  octree.maxGridCellAtLevel(lev),
                  octree.spacingAtLevel(lev),
                  octree.blockBoundingBox(gridPt, lev)));
  }

  SLIC_INFO("Recursively refining around query point: " << queryPt);
  refineAndPrint(octree, queryPt, false);
  //for(int i=0; i< octree.maxInternalLevel(); ++i)
  //  refineAndPrint(octree, queryPt);

  // Reclaim memory
  if(surface_mesh != nullptr)
  {
    delete surface_mesh;
    surface_mesh = nullptr;
  }

  return 0;
}
