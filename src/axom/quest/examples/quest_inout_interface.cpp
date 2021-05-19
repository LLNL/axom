// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file quest_inout_interface.cpp
 *
 * \brief Simple example that exercises the C-style quest_inout interface
 * which determines if a point is contained within an enclosed volume
 * defined by a surface mesh.
 */

// Axom includes
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"

// _quest_inout_interface_include_start
#include "axom/quest/interface/inout.hpp"

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#endif
// _quest_inout_interface_include_end

// C/C++ includes
#include <vector>

namespace slic = axom::slic;
namespace primal = axom::primal;
namespace quest = axom::quest;

typedef primal::Point<double, 3> PointType;
typedef primal::BoundingBox<double, 3> BoxType;
typedef std::vector<PointType> CoordsVec;

//------------------------------------------------------------------------------

/*!
 * \brief Utility function to generate a random set of query points
 * within a given bounding box
 *
 * \param [out] queryPoints The generated set of query points
 * \param [in] bbox The bounding box for the query points
 * \param [in] numPoints The number of points to generate
 */
void generateQueryPoints(CoordsVec& queryPoints, BoxType const& bbox, int numPoints)
{
  queryPoints.clear();
  queryPoints.reserve(numPoints);
  for(int i = 0; i < numPoints; ++i)
  {
    double x = axom::utilities::random_real(bbox.getMin()[0], bbox.getMax()[0]);
    double y = axom::utilities::random_real(bbox.getMin()[1], bbox.getMax()[1]);
    double z = axom::utilities::random_real(bbox.getMin()[2], bbox.getMax()[2]);
    queryPoints.push_back(PointType::make_point(x, y, z));
  }
}

/*!
 * \brief Utility function to initialize the logger
 */
void initializeLogger()
{
  // Initialize Logger
  slic::initialize();
  slic::setLoggingMsgLevel(axom::slic::message::Info);

  slic::LogStream* logStream;

#ifdef AXOM_USE_MPI
  std::string fmt = "[<RANK>][<LEVEL>]: <MESSAGE>\n";
  #ifdef AXOM_USE_LUMBERJACK
  const int RLIMIT = 8;
  logStream = new slic::LumberjackStream(&std::cout, MPI_COMM_WORLD, RLIMIT, fmt);
  #else
  logStream = new slic::SynchronizedStream(&std::cout, MPI_COMM_WORLD, fmt);
  #endif
#else
  std::string fmt = "[<LEVEL>]: <MESSAGE>\n";
  logStream = new slic::GenericOutputStream(&std::cout, fmt);
#endif  // AXOM_USE_MPI

  slic::addStreamToAllMsgLevels(logStream);
}

/*!
 * \brief Utility function to finalize the logger
 */
void finalizeLogger()
{
  if(slic::isInitialized())
  {
    slic::flushStreams();
    slic::finalize();
  }
}

/*!
 * \brief Utility function for a clean abort that finalizes
 * quest_inout, slic and mpi
 */
void cleanAbort()
{
  if(quest::inout_initialized())
  {
    quest::inout_finalize();
  }

  finalizeLogger();

  axom::utilities::processAbort();
}

/*!
 * \brief A simple example illustrating the use of the Quest C-Style interface.
 *
 * \note To run the example
 * \verbatim
 *
 *   [mpirun -np N] ./quest_inout_interface_ex <stl_file>
 *
 * \endverbatim
 */
int main(int argc, char** argv)
{
  // -- Initialize MPI
#ifdef AXOM_USE_MPI
  MPI_Init(&argc, &argv);
#endif

  // -- Initialize logger
  initializeLogger();

  // -- Parse command line options
  if(argc != 2)
  {
#ifdef AXOM_USE_MPI
    SLIC_WARNING("Usage: [mpirun -np N] ./quest_inout_interface_ex <stl_file>");
#else
    SLIC_WARNING("Usage: ./quest_inout_interface_ex <stl_file>");
#endif
    cleanAbort();
  }

  std::string fileName = std::string(argv[1]);
  bool isVerbose = false;
  const int npoints = 100000;
  double weldThresh = 1E-9;

  int rc = quest::QUEST_INOUT_SUCCESS;

  // _quest_inout_interface_parameters_start
  // -- Set quest_inout parameters
  rc = quest::inout_set_verbose(isVerbose);
  if(rc != quest::QUEST_INOUT_SUCCESS)
  {
    cleanAbort();
  }

  rc = quest::inout_set_vertex_weld_threshold(weldThresh);
  if(rc != quest::QUEST_INOUT_SUCCESS)
  {
    cleanAbort();
  }
  // _quest_inout_interface_parameters_end

  // -- Initialize quest_inout
  axom::utilities::Timer timer;
  SLIC_INFO("Initializing quest_inout...");
  timer.start();
  {
// _quest_inout_interface_init_start
#ifdef AXOM_USE_MPI
    rc = quest::inout_init(fileName, MPI_COMM_WORLD);
#else
    rc = quest::inout_init(fileName);
#endif
    // _quest_inout_interface_init_end
  }
  timer.stop();
  if(rc != quest::QUEST_INOUT_SUCCESS)
  {
    cleanAbort();
  }
  else
  {
    SLIC_INFO("  initialization took " << timer.elapsed() << " seconds.");
  }

  // -- Query mesh bounding box and center of mass
  double bbMin[3], bbMax[3], cMass[3];
  rc = quest::inout_mesh_min_bounds(bbMin);
  if(rc != quest::QUEST_INOUT_SUCCESS)
  {
    cleanAbort();
  }

  rc = quest::inout_mesh_max_bounds(bbMax);
  if(rc != quest::QUEST_INOUT_SUCCESS)
  {
    cleanAbort();
  }

  rc = quest::inout_mesh_center_of_mass(cMass);
  if(rc != quest::QUEST_INOUT_SUCCESS)
  {
    cleanAbort();
  }

  BoxType bbox = BoxType(PointType(bbMin), PointType(bbMax));
  SLIC_INFO("Mesh bounding box: " << bbox);
  SLIC_INFO("Mesh center of mass: " << PointType(cMass));

  // -- Generate query points
  CoordsVec queryPoints;
  generateQueryPoints(queryPoints, bbox, npoints);

  // -- Run the queries
  int numInside = 0;
  SLIC_INFO("Querying mesh with " << npoints << " query points...");
  timer.start();
  for(auto& pt : queryPoints)
  {
    // _quest_inout_interface_test_start
    const double x = pt[0];
    const double y = pt[1];
    const double z = pt[2];

    const bool ins = quest::inout_evaluate(x, y, z);
    numInside += ins ? 1 : 0;
    // _quest_inout_interface_test_end
  }
  timer.stop();

  // -- Output some query statistics
  {
    SLIC_INFO("  queries took " << timer.elapsed() << " seconds.");

    const double rate = queryPoints.size() / timer.elapsed();
    SLIC_INFO("  query rate: " << rate << " queries per second.");

    const double intPercent =
      (100 * numInside) / static_cast<double>(queryPoints.size());
    SLIC_INFO("  " << numInside << " of " << queryPoints.size() << " ("
                   << intPercent << "%) "
                   << " of the query points were contained in the surface.");
  }

  // -- Finalize quest_inout
  // _quest_inout_interface_finalize_start
  quest::inout_finalize();
  // _quest_inout_interface_finalize_end

  // -- Finalize logger
  finalizeLogger();

  // -- Finalize MPI
#ifdef AXOM_USE_MPI
  MPI_Finalize();
#endif

  return 0;
}
