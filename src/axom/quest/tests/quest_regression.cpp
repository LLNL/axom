// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file
 * This file generates and runs regression tests for the Quest signed distance
 * and point containment queries. It generates the signed distance
 * representation and/or the InOutOctree representation and queries it over a
 * uniform grid of a given resolution and bounding box using the quest C
 * interface.
 *
 * The baseline files are stored in a Sidre datastore
 * with the following structure:
 *  /mesh_name          (string: name of mesh file, without paths)
 *  /mesh_bounding_box  (a doubles: min_x, min_y, min_z, max_x, max_y, max_z)
 *  /query_resolution   (3 ints: i, j, k of query grid)
 *  /octree_containment (ints: one for each query point, with value 0 or 1)
 *  /bvh_containment    (ints: one for each query point, with value 0 or 1)
 *  /bvh_distance       (doubles: one for each query point,
 *                       value is min signed distance from the associated point
 *                       to the surface).
 *
 *  The 'octree_containment' field is generated by default, or when the
 *  '--containment' command line option is present.  Similarly for the
 *  'bvh_containment' and 'bvh_distance' and the '--distance' command line
 *  option.
 *
 * Usage: run with "-h" or "--help" to see command line options
 */

// axom includes
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/sidre.hpp"
#include "axom/mint.hpp"

#include "axom/quest/interface/inout.hpp"
#include "axom/quest/interface/signed_distance.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

// MPI includes
#include "mpi.h"

#ifdef AXOM_USE_OPENMP
  #include "omp.h"
#endif

// C/C++ includes
#include <cmath>
#include <string>

constexpr int DIM = 3;
// Default resolution of query grid
constexpr int DEFAULT_RESOLUTION = 32;
// Max number of disagreeing entries to show when comparing results
constexpr int MAX_RESULTS = 10;

namespace mint = axom::mint;
namespace primal = axom::primal;
namespace quest = axom::quest;
namespace sidre = axom::sidre;
namespace utilities = axom::utilities;

using SpaceBoundingBox = primal::BoundingBox<double, DIM>;
using SpacePt = primal::Point<double, DIM>;
using SpaceVec = primal::Vector<double, DIM>;
using GridPt = primal::Point<int, DIM>;

/// Simple structure to hold the command line arguments
struct Input
{
  Input() = default;

  ~Input()
  {
    if(queryMesh != nullptr)
    {
      delete queryMesh;
      queryMesh = nullptr;
    }
  }

  std::string meshName;
  std::string baselineRoot;

  std::vector<int> m_resolution;
  std::vector<double> m_queryBoxMins;
  std::vector<double> m_queryBoxMaxs;

  SpaceBoundingBox meshBoundingBox;
  GridPt queryResolution {DEFAULT_RESOLUTION};

  mint::UniformMesh* queryMesh {nullptr};

  bool testDistance {true};
  bool testContainment {true};

  bool hasBaseline() const { return !baselineRoot.empty(); }
  bool hasBoundingBox() const { return meshBoundingBox != SpaceBoundingBox(); }
  bool hasQueryMesh() const { return queryMesh != nullptr; }

  /// Parses the command line options
  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    app.add_option("-m,--mesh", meshName, "Surface mesh file (STL files are currently supported)")
      ->required()
      ->check(axom::CLI::ExistingFile);

    app
      .add_flag("--distance,!--no-distance",
                testDistance,
                "Indicates whether to test the signed distance")
      ->capture_default_str();
    app
      .add_flag("--containment,!--no-containment",
                testContainment,
                "Indicates whether to test the point containment")
      ->capture_default_str();

    // Note: Baselines comparisons only supported when Axom is configured with hdf5
    // Users must supply either a baseline, or both the query resolution and bounding box
#ifdef AXOM_USE_HDF5
    app
      .add_option("-b,--baseline",
                  baselineRoot,
                  "root file of baseline, a sidre rootfile.\n"
                  "Note: Only supported when Axom configured with hdf5")
      ->check(axom::CLI::ExistingFile);
#endif

    // user can supply 1 or 3 values for resolution the following is not quite right
    auto* res = app
                  .add_option("-r,--resolution",
                              m_resolution,
                              "Resolution of the sample grid.\n"
                              "Note: Only used when baseline not supplied.")
                  ->expected(-1);

    // Optional bounding box for query region
    auto* minbb = app
                    .add_option("--min",
                                m_queryBoxMins,
                                "Min bounds for query box (x,y,z).\n"
                                "Note: Only used when baseline not supplied.")
                    ->expected(3);
    auto* maxbb = app
                    .add_option("--max",
                                m_queryBoxMaxs,
                                "Max bounds for query box (x,y,z)\n"
                                "Note: Only used when baseline not supplied.")
                    ->expected(3);

    // Add some requirements
    minbb->needs(maxbb);
    minbb->needs(res);
    maxbb->needs(minbb);
    maxbb->needs(res);
    res->needs(minbb);
    res->needs(maxbb);

    app.get_formatter()->column_width(45);

    // could throw an exception
    app.parse(argc, argv);

    // Check if resolution or bounding box values were provided in addition
    // to baseline. If so, inform user that these will be overridden by baseline
    bool hasRes = !m_resolution.empty();
    bool hasBBox = !m_queryBoxMins.empty();
    if(!baselineRoot.empty())
    {
      if(hasRes || hasBBox)
      {
        SLIC_INFO("Baseline mesh will override values for resolution and bounding box");
      }
    }
    else if(!hasRes || !hasBBox)
    {
      throw axom::CLI::Error("IncorrectConstruction",
                             "Resolution and bounding box required when baseline is not supplied");
    }
    else
    {
      const int dim = 3;
      queryResolution =
        (m_resolution.size() == 1) ? GridPt(m_resolution[0]) : GridPt(m_resolution.data(), dim);
      meshBoundingBox.addPoint(SpacePt(m_queryBoxMins.data(), dim));
      meshBoundingBox.addPoint(SpacePt(m_queryBoxMaxs.data(), dim));
    }

    if(!testContainment && !testDistance)
    {
      throw axom::CLI::Error("IncorrectConstruction",
                             "At least one of {--distance; --containment} must be enabled.");
    }
  }
};

/// Loads the baseline dataset into the given sidre group
void loadBaselineData(sidre::Group* grp, Input& args)
{
  sidre::IOManager reader(MPI_COMM_WORLD);
  reader.read(grp, args.baselineRoot, "sidre_hdf5");

  /// Check that the required fields are present

  if(!grp->hasView("mesh_name"))
  {
    SLIC_ERROR("Baseline must include a 'mesh_name' view");
  }

  // Check for bounding box, and load into the args instance
  if(!grp->hasView("mesh_bounding_box"))
  {
    SLIC_ERROR("Baseline must include a 'mesh_bounding_box' view");
  }
  else
  {
    sidre::View* view = grp->getView("mesh_bounding_box");
    if(view->getNumElements() != 6)
    {
      SLIC_ERROR("Bounding box must contain six doubles");
    }

    double* data = view->getData();
    args.meshBoundingBox = SpaceBoundingBox(SpacePt(data, 3), SpacePt(data + 3, 3));
  }

  // Check for query grid resolution, and load into the args instance
  if(!grp->hasView("query_resolution"))
  {
    SLIC_ERROR("Baseline must include a 'query_resolution' view");
  }
  else
  {
    sidre::View* view = grp->getView("query_resolution");
    if(view->getNumElements() != 3)
    {
      SLIC_ERROR("Query resolution must contain three ints");
    }

    int* data = view->getData();
    args.queryResolution = GridPt(data, 3);
  }

  // Optionally check for the InOutOctree point containment data
  if(args.testContainment)
  {
    if(!grp->hasView("octree_containment"))
    {
      SLIC_ERROR("Requested containment, but baseline "
                 << "does not have a 'octree_containment' view");
    }
    else
    {
      SLIC_ASSERT_MSG(grp->getView("octree_containment")->getTypeID() == sidre::INT_ID,
                      "Type of 'octree_containment' view must be int (SIDRE_INT_ID)");
    }
  }

  // Optionally check for the Signed distance point containment and distance data
  if(args.testDistance)
  {
    if(!grp->hasView("bvh_distance"))
    {
      SLIC_ERROR("Requested distance, but baseline "
                 << "does not have a 'bvh_distance' view");
    }
    else
    {
      SLIC_ASSERT_MSG(grp->getView("bvh_distance")->getTypeID() == sidre::DOUBLE_ID,
                      "Type of 'bvh_distance' view must be double (SIDRE_DOUBLE_ID)");
    }

    if(!grp->hasView("bvh_containment"))
    {
      SLIC_ERROR("Requested distance, but baseline does not "
                 << " have a 'bvh_containment' view");
    }
    else
    {
      SLIC_ASSERT_MSG(grp->getView("bvh_containment")->getTypeID() == sidre::INT_ID,
                      "Type of 'bvh_containment' view must be int (SIDRE_INT_ID)");
    }
  }
}

/**
 * \brief Generates a mint Uniform mesh with the given bounding box and resolution
 * \note Allocates a UniformMesh instance, which must be deleted by the user
 */
mint::UniformMesh* createQueryMesh(const SpaceBoundingBox& bb, const GridPt& res)
{
  const double* low = bb.getMin().data();
  const double* high = bb.getMax().data();

  return new mint::UniformMesh(low, high, res[0] + 1, res[1] + 1, res[2] + 1);
}

/**
 * Runs the InOutOctree point containment queries and adds results as scalar
 * field on uniform mesh
 */
void runContainmentQueries(Input& clargs)
{
  SLIC_INFO(axom::fmt::format("Initializing InOutOctree over mesh '{}'...", clargs.meshName));
  utilities::Timer buildTimer(true);

  quest::inout_init(clargs.meshName, MPI_COMM_WORLD);

  buildTimer.stop();
  SLIC_INFO(axom::fmt::format("Initialization took {} seconds.", buildTimer.elapsed()));

  SpacePt bbMin;
  SpacePt bbMax;
  quest::inout_mesh_min_bounds(bbMin.data());
  quest::inout_mesh_max_bounds(bbMax.data());

  if(!clargs.hasBoundingBox())
  {
    clargs.meshBoundingBox = SpaceBoundingBox(bbMin, bbMax);
    clargs.meshBoundingBox.scale(1.5);
  }

  if(!clargs.hasQueryMesh())
  {
    clargs.queryMesh = createQueryMesh(clargs.meshBoundingBox, clargs.queryResolution);
  }

  SLIC_INFO("Mesh bounding box is: " << SpaceBoundingBox(bbMin, bbMax));
  SLIC_INFO("Query bounding box is: " << clargs.meshBoundingBox);

#ifdef AXOM_USE_OPENMP
  #pragma omp parallel
  #pragma omp master
  SLIC_INFO(
    axom::fmt::format("Querying InOutOctree on uniform grid "
                      "of resolution {} using {} threads",
                      clargs.queryResolution,
                      omp_get_num_threads()));
#else
  SLIC_INFO(axom::fmt::format("Querying InOutOctree on uniform grid of resolution {}",
                              clargs.queryResolution));
#endif

  // Add a scalar field for the containment queries
  SLIC_ASSERT(clargs.queryMesh != nullptr);
  mint::UniformMesh* umesh = clargs.queryMesh;
  const axom::IndexType nnodes = umesh->getNumberOfNodes();

  int* containment = umesh->createField<int>("octree_containment", mint::NODE_CENTERED);
  SLIC_ASSERT(containment != nullptr);

  double* xcoords = new double[nnodes];
  double* ycoords = new double[nnodes];
  double* zcoords = new double[nnodes];
  utilities::Timer fillTimer(true);

#pragma omp parallel for schedule(static)
  for(int inode = 0; inode < nnodes; ++inode)
  {
    axom::IndexType i, j, k;
    umesh->getNodeGridIndex(inode, i, j, k);

    xcoords[inode] = umesh->evaluateCoordinate(i, mint::X_COORDINATE);
    ycoords[inode] = umesh->evaluateCoordinate(j, mint::Y_COORDINATE);
    zcoords[inode] = umesh->evaluateCoordinate(k, mint::Z_COORDINATE);
  }
  fillTimer.stop();

  utilities::Timer queryTimer(true);
  quest::inout_evaluate(xcoords, ycoords, zcoords, nnodes, containment);
  queryTimer.stop();

  SLIC_INFO(axom::fmt::format("Filling coordinates array took {} seconds", fillTimer.elapsed()));
  SLIC_INFO(
    axom::fmt::format("Querying {}^3 containment field (InOutOctree) "
                      "took {} seconds (@ {} queries per second)",
                      clargs.queryResolution,
                      queryTimer.elapsed(),
                      nnodes / queryTimer.elapsed()));

  delete[] xcoords;
  delete[] ycoords;
  delete[] zcoords;

  quest::inout_finalize();
}

/**
 * Runs the SignedDistance point containment and distance queries and adds
 * results as scalar field on uniform mesh
 */
void runDistanceQueries(Input& clargs)
{
  SLIC_INFO(axom::fmt::format("Initializing linear BVH over mesh '{}'...", clargs.meshName));
  utilities::Timer buildTimer(true);

  quest::signed_distance_init(clargs.meshName, MPI_COMM_WORLD);

  buildTimer.stop();

  SLIC_INFO(axom::fmt::format("Initialization took {} seconds.", buildTimer.elapsed()));

  SpacePt bbMin, bbMax;
  quest::signed_distance_get_mesh_bounds(bbMin.data(), bbMax.data());

  if(!clargs.hasBoundingBox())
  {
    clargs.meshBoundingBox = SpaceBoundingBox(bbMin, bbMax);
    clargs.meshBoundingBox.scale(1.5);
  }

  if(!clargs.hasQueryMesh())
  {
    clargs.queryMesh = createQueryMesh(clargs.meshBoundingBox, clargs.queryResolution);
  }

  SLIC_INFO("Mesh bounding box is: " << SpaceBoundingBox(bbMin, bbMax));
  SLIC_INFO("Query bounding box is: " << clargs.meshBoundingBox);

#ifdef AXOM_USE_OPENMP
  #pragma omp parallel
  #pragma omp master
  SLIC_INFO(
    axom::fmt::format("Querying BVH on uniform grid "
                      "of resolution {} using {} threads",
                      clargs.queryResolution,
                      omp_get_num_threads()));
#else
  SLIC_INFO(
    axom::fmt::format("Querying BVH on uniform grid of resolution {}", clargs.queryResolution));
#endif

  // Add a scalar field for the containment queries
  SLIC_ASSERT(clargs.queryMesh != nullptr);
  axom::mint::UniformMesh* umesh = clargs.queryMesh;
  const int nnodes = umesh->getNumberOfNodes();

  int* containment = umesh->createField<int>("bvh_containment", axom::mint::NODE_CENTERED);
  double* distance = umesh->createField<double>("bvh_distance", axom::mint::NODE_CENTERED);
  SLIC_ASSERT(containment != nullptr);
  SLIC_ASSERT(distance != nullptr);

  double* xcoords = new double[nnodes];
  double* ycoords = new double[nnodes];
  double* zcoords = new double[nnodes];
  utilities::Timer fillTimer(true);

#pragma omp parallel for schedule(static)
  for(int inode = 0; inode < nnodes; ++inode)
  {
    axom::IndexType i, j, k;
    umesh->getNodeGridIndex(inode, i, j, k);

    xcoords[inode] = umesh->evaluateCoordinate(i, axom::mint::X_COORDINATE);
    ycoords[inode] = umesh->evaluateCoordinate(j, axom::mint::Y_COORDINATE);
    zcoords[inode] = umesh->evaluateCoordinate(k, axom::mint::Z_COORDINATE);
  }
  fillTimer.stop();

  utilities::Timer distanceTimer(true);
  quest::signed_distance_evaluate(xcoords, ycoords, zcoords, nnodes, distance);
  distanceTimer.stop();

  for(int inode = 0; inode < nnodes; ++inode)
  {
    containment[inode] = (std::signbit(distance[inode]) != 0) ? 1 : 0;
  }

  SLIC_INFO(axom::fmt::format("Filling coordinates array took {} seconds", fillTimer.elapsed()));
  SLIC_INFO(
    axom::fmt::format("Querying {}^3 signed distance field (BVH) "
                      "took {} seconds (@ {} queries per second)",
                      clargs.queryResolution,
                      distanceTimer.elapsed(),
                      nnodes / distanceTimer.elapsed()));

  delete[] xcoords;
  delete[] ycoords;
  delete[] zcoords;

  quest::signed_distance_finalize();
}

/**
 * \brief Function to compare the results from the InOutOctree and
 * SignedDistance queries
 * \return True if all results agree, False otherwise.
 * \note When there are differences, the first few are logged
 */
bool compareDistanceAndContainment(Input& clargs)
{
  SLIC_ASSERT(clargs.hasQueryMesh());

  bool passed = true;

  mint::UniformMesh* umesh = clargs.queryMesh;
  const int nnodes = umesh->getNumberOfNodes();

  if(!clargs.testContainment)
  {
    SLIC_INFO("Cannot compare signed distance and InOutOctree "
              << "-- InOutOctree was not generated");
  }

  if(!clargs.testDistance)
  {
    SLIC_INFO("Cannot compare signed distance and InOutOctree "
              << "-- Signed distance was not generated");
  }

  if(clargs.testContainment && clargs.testDistance)
  {
    // compare containment results of the two approaches
    int diffCount = 0;
    axom::fmt::memory_buffer out;

    int* oct_containment = umesh->getFieldPtr<int>("octree_containment", mint::NODE_CENTERED);
    int* bvh_containment = umesh->getFieldPtr<int>("bvh_containment", mint::NODE_CENTERED);
    double* bvh_distance = umesh->getFieldPtr<double>("bvh_distance", mint::NODE_CENTERED);

    for(int inode = 0; inode < nnodes; ++inode)
    {
      const int oct_c = oct_containment[inode];
      const int bvh_c = bvh_containment[inode];

      if(bvh_c != oct_c)
      {
        // allow for differences between the two approaches really close to the boundary
        const double bvh_d = bvh_distance[inode];
        if(!axom::utilities::isNearlyEqual(bvh_d, 0.))
        {
          if(diffCount < MAX_RESULTS)
          {
            primal::Point<double, 3> pt;
            umesh->getNode(inode, pt.data());

            axom::fmt::format_to(std::back_inserter(out),
                                 "\n  Disagreement on sample {} @ {}.  "
                                 "Signed distance: {} ({}) -- InOutOctree: {} ",
                                 inode,
                                 pt,
                                 bvh_d,
                                 bvh_c ? "inside" : "outside",
                                 oct_c ? "inside" : "outside");
          }
          ++diffCount;
        }
      }
    }

    if(diffCount != 0)
    {
      passed = false;
      SLIC_INFO(
        axom::fmt::format("** Disagreement between SignedDistance and InOutOctree queries.\n"
                          "There were {} differences.\n Showing first {} results -- {}",
                          diffCount,
                          std::min(diffCount, MAX_RESULTS),
                          std::string(out.begin(), out.end())));
    }
  }

  return passed;
}

/**
 * \brief Function to compare the results from the InOutOctree and SignedDistance queries
 * \return True if all results agree, False otherwise.
 * \note When there are differences, the first few are logged
 */
bool compareToBaselineResults(axom::sidre::Group* grp, Input& clargs)
{
  SLIC_ASSERT(grp != nullptr);
  SLIC_ASSERT(clargs.hasQueryMesh());

  bool passed = true;

  mint::UniformMesh* umesh = clargs.queryMesh;
  const int nnodes = umesh->getNumberOfNodes();

  int* base_inout_containment = nullptr;
  int* exp_inout_containment = nullptr;
  int* base_sd_containment = nullptr;
  int* exp_sd_containment = nullptr;
  double* base_sd_distance = nullptr;
  double* exp_sd_distance = nullptr;

  // grab pointers to data arrays when appropriate
  if(clargs.testContainment)
  {
    base_inout_containment = grp->getView("octree_containment")->getArray();
    exp_inout_containment = umesh->getFieldPtr<int>("octree_containment", mint::NODE_CENTERED);
  }

  if(clargs.testDistance)
  {
    base_sd_containment = grp->getView("bvh_containment")->getArray();
    exp_sd_containment = umesh->getFieldPtr<int>("bvh_containment", mint::NODE_CENTERED);

    base_sd_distance = grp->getView("bvh_distance")->getArray();
    exp_sd_distance = umesh->getFieldPtr<double>("bvh_distance", mint::NODE_CENTERED);
  }

  if(clargs.testContainment)
  {
    int diffCount = 0;
    axom::fmt::memory_buffer out;

    for(int inode = 0; inode < nnodes; ++inode)
    {
      const int expected = base_inout_containment[inode];
      const int actual = exp_inout_containment[inode];
      if(expected != actual)
      {
        if(diffCount < MAX_RESULTS)
        {
          primal::Point<double, 3> pt;
          umesh->getNode(inode, pt.data());

          axom::fmt::format_to(std::back_inserter(out),
                               "\n  Disagreement on sample {} @ {}.  Expected {}, got {}",
                               inode,
                               pt,
                               expected ? "inside" : "outside",
                               actual ? "inside" : "outside");

          if(clargs.testDistance)
          {
            axom::fmt::format_to(std::back_inserter(out),
                                 " -- distance from query point to surface is {} ",
                                 base_sd_distance[inode]);
          }
        }
        ++diffCount;
      }
    }

    if(diffCount != 0)
    {
      passed = false;
      SLIC_INFO(
        axom::fmt::format("** Containment test failed.  There were {} differences. "
                          "Showing first {} -- {}",
                          diffCount,
                          std::min(diffCount, MAX_RESULTS),
                          std::string(out.begin(), out.end())));
    }
  }

  if(clargs.testDistance)
  {
    int diffCount = 0;
    axom::fmt::memory_buffer out;

    for(int inode = 0; inode < nnodes; ++inode)
    {
      const int expected_c = base_sd_containment[inode];
      const int actual_c = exp_sd_containment[inode];
      const double expected_d = base_sd_distance[inode];
      const double actual_d = exp_sd_distance[inode];
      if(expected_c != actual_c || !utilities::isNearlyEqual(expected_d, actual_d))
      {
        if(diffCount < MAX_RESULTS)
        {
          primal::Point<double, 3> pt;
          umesh->getNode(inode, pt.data());

          axom::fmt::format_to(std::back_inserter(out),
                               "\n  Disagreement on sample {} @ {}. Expected {} ({}), got {} ({})",
                               inode,
                               pt,
                               expected_d,
                               expected_c ? "inside" : "outside",
                               actual_d,
                               actual_c ? "inside" : "outside");
        }
        ++diffCount;
      }
    }

    if(diffCount != 0)
    {
      passed = false;
      SLIC_INFO(
        axom::fmt::format("** Distance test failed.  There were {} differences. "
                          "Showing first {} -- {}",
                          diffCount,
                          std::min(diffCount, MAX_RESULTS),
                          std::string(out.begin(), out.end())));
    }
  }

  return passed;
}

/**
 * \brief Saves the current results as a new baseline
 * \note The baseline will be output as
 *       a sidre root file ./<mesh>_<res>_baseline.root
 *       and a corresponding folder ./<mesh>_<res>_baseline/
 *       (both in the same directory)
 */
void saveBaseline(axom::sidre::Group* grp, Input& clargs)
{
  SLIC_ASSERT(grp != nullptr);
  SLIC_ASSERT(clargs.hasQueryMesh());

  std::string fullMeshName = clargs.meshName;
  std::size_t found = fullMeshName.find_last_of("/");
  std::string meshName = fullMeshName.substr(found + 1);

  found = meshName.find_last_of(".");
  std::string meshNameNoExt = meshName.substr(0, found);

  grp->createViewString("mesh_name", meshName);

  sidre::View* view = nullptr;

  view = grp->createView("mesh_bounding_box", sidre::DOUBLE_ID, 6)->allocate();
  double* bb = view->getArray();
  clargs.meshBoundingBox.getMin().to_array(bb);
  clargs.meshBoundingBox.getMax().to_array(bb + 3);

  view = grp->createView("query_resolution", sidre::INT_ID, 3)->allocate();
  clargs.queryResolution.to_array(view->getArray());

  axom::mint::UniformMesh* umesh = clargs.queryMesh;
  const int nnodes = umesh->getNumberOfNodes();
  if(clargs.testContainment)
  {
    int* oct_containment = umesh->getFieldPtr<int>("octree_containment", mint::NODE_CENTERED);
    view = grp->createView("octree_containment", sidre::INT_ID, nnodes)->allocate();
    int* contData = view->getArray();
    std::copy(oct_containment, oct_containment + nnodes, contData);
  }

  if(clargs.testDistance)
  {
    int* bvh_containment = umesh->getFieldPtr<int>("bvh_containment", mint::NODE_CENTERED);
    view = grp->createView("bvh_containment", sidre::INT_ID, nnodes)->allocate();
    int* contData = view->getArray();
    std::copy(bvh_containment, bvh_containment + nnodes, contData);

    double* bvh_distance = umesh->getFieldPtr<double>("bvh_distance", mint::NODE_CENTERED);
    view = grp->createView("bvh_distance", sidre::DOUBLE_ID, nnodes)->allocate();

    double* distData = view->getArray();
    std::copy(bvh_distance, bvh_distance + nnodes, distData);
  }

  const GridPt& res = clargs.queryResolution;
  bool resAllSame = (res[0] == res[1] && res[1] == res[2]);
  std::string resStr = resAllSame ? axom::fmt::format("{}", res[0])
                                  : axom::fmt::format("{}_{}_{}", res[0], res[1], res[2]);

  std::string outfile = axom::fmt::format("{}_{}_{}", meshNameNoExt, resStr, "baseline");
  std::string protocol = "sidre_hdf5";
  sidre::IOManager writer(MPI_COMM_WORLD);
  writer.write(grp, 1, outfile, protocol);
  SLIC_INFO(axom::fmt::format("** Saved baseline file '{}' using '{}' protocol.", outfile, protocol));
}

/// Runs regression test for quest containment and signed distance queries
int main(int argc, char** argv)
{
  bool allTestsPassed = true;

  // initialize the problem
  MPI_Init(&argc, &argv);

  axom::slic::SimpleLogger logger;
  sidre::DataStore ds;

  // parse the command arguments
  Input args;
  axom::CLI::App app {"Regression tester for point containment and signed distance tests"};

  try
  {
    args.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    return app.exit(e);
  }

#ifdef AXOM_USE_HDF5
  // load the baseline file for comparisons and additional test parameters
  // This is currently only supported when hdf5 is enabled.
  if(args.hasBaseline())
  {
    loadBaselineData(ds.getRoot(), args);
  }
#endif

  // run the containment queries
  if(args.testContainment)
  {
    SLIC_INFO("Running containment queries");
    runContainmentQueries(args);
    SLIC_INFO("--");
  }

  // run the distance queries
  if(args.testDistance)
  {
    SLIC_INFO("Running distance queries");
    runDistanceQueries(args);
    SLIC_INFO("--");
  }

  // Compare signs of current results on SignedDistance and InOutOctree
  bool methodsAgree = true;
  if(args.testContainment && args.testDistance)
  {
    SLIC_INFO("Comparing results from containment and distance queries");

    methodsAgree = compareDistanceAndContainment(args);

    SLIC_INFO("** Methods " << (methodsAgree ? "agree" : "do not agree"));

    allTestsPassed = allTestsPassed && methodsAgree;

    SLIC_INFO("--");
  }

#ifdef AXOM_USE_HDF5
  // compare current results to baselines or generate new baselines
  // This is currently only supported when hdf5 is enabled.
  bool baselinePassed = true;
  if(args.hasBaseline())
  {
    SLIC_INFO("Comparing results to baselines");
    baselinePassed = compareToBaselineResults(ds.getRoot(), args);

    SLIC_INFO("** Baseline tests " << (baselinePassed ? "passed" : "failed"));
    allTestsPassed = allTestsPassed && baselinePassed;
  }
  else
  {
    SLIC_INFO("Saving results as new baseline.");
    saveBaseline(ds.getRoot(), args);
  }
  SLIC_INFO("--");
#endif

  // finalize
  MPI_Finalize();
  return (allTestsPassed) ? 0 : 1;
}
