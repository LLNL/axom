// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core/utilities/Timer.hpp"
#include "axom/quest/PointInCell.hpp"

#include "axom/slic.hpp"

#include "axom/CLI11.hpp"
#include "axom/fmt.hpp"

#ifdef AXOM_USE_MFEM
  #include "axom/quest/detail/PointInCellMeshWrapper_mfem.hpp"
#else
  #error "Quest's PointInCell on mfem meshes requires mfem library."
#endif

// namespace aliases
namespace primal = axom::primal;
namespace quest = axom::quest;
namespace slic = axom::slic;
namespace utilities = axom::utilities;

enum class ExecPolicy
{
  CPU,
  OpenMP,
  GPU
};

const std::map<std::string, ExecPolicy> validExecPolicies {
  {"seq", ExecPolicy::CPU},
#ifdef AXOM_USE_OPENMP
  {"omp", ExecPolicy::OpenMP},
#endif
#ifdef AXOM_USE_CUDA
  {"gpu", ExecPolicy::GPU}
#endif
};

void initialize_logger();
void finalize_logger();

template <int NDIMS>
primal::Point<double, NDIMS> get_rand_pt(
  const primal::BoundingBox<double, NDIMS>& bbox)
{
  primal::Point<double, NDIMS> rnd_pt;
  for(int i = 0; i < NDIMS; i++)
  {
    rnd_pt[i] = axom::utilities::random_real(bbox.getMin()[i], bbox.getMax()[i]);
  }
  return rnd_pt;
}

template <typename ExecSpace>
void benchmark_point_in_cell(mfem::Mesh& mesh, int npts)
{
  if(mesh.SpaceDimension() == 2)
  {
    using PointType = primal::Point<double, 2>;
    using BoxType = primal::BoundingBox<double, 2>;
    using mesh_tag = axom::quest::quest_point_in_cell_mfem_tag;
    using IndexType = typename quest::PointInCellTraits<mesh_tag>::IndexType;

    BoxType meshBb;
    {
      mfem::Vector meshMin, meshMax;
      mesh.GetBoundingBox(meshMin, meshMax);
      meshBb.addPoint(PointType(meshMin.GetData()));
      meshBb.addPoint(PointType(meshMax.GetData()));
    }

    PointType* pts = axom::allocate<PointType>(npts);

    utilities::Timer timeInitRandPts(true);
    // Generate random points
    for(int i = 0; i < npts; i++)
    {
      pts[i] = get_rand_pt(meshBb);
    }
    SLIC_INFO(axom::fmt::format("Constructed {} random points in {} s.",
                                npts,
                                timeInitRandPts.elapsed()));

    primal::Point<int, 2> bins(25);

    utilities::Timer timeInitQuery(true);
    // Create PointInCell query object
    quest::PointInCell<mesh_tag, ExecSpace> query(&mesh, bins.data());
    SLIC_INFO(axom::fmt::format("Initialized point-in-cell query in {} s.",
                                timeInitQuery.elapsed()));

    IndexType* outCellIds = axom::allocate<IndexType>(npts);
    // Run query
    utilities::Timer timeRunQuery(true);
    query.locatePoints(npts, pts, outCellIds);
    double time = timeRunQuery.elapsed();
    SLIC_INFO(
      axom::fmt::format("Ran query on {} points in {} s -- rate: {} q/s",
                        npts,
                        time,
                        npts / time));

    axom::deallocate(pts);
  }
  else if(mesh.SpaceDimension() == 3)
  {
    using PointType = primal::Point<double, 3>;
    using BoxType = primal::BoundingBox<double, 3>;
    using mesh_tag = axom::quest::quest_point_in_cell_mfem_tag;
    using IndexType = typename quest::PointInCellTraits<mesh_tag>::IndexType;

    BoxType meshBb;
    {
      mfem::Vector meshMin, meshMax;
      mesh.GetBoundingBox(meshMin, meshMax);
      meshBb.addPoint(PointType(meshMin.GetData()));
      meshBb.addPoint(PointType(meshMax.GetData()));
    }

    PointType* pts = axom::allocate<PointType>(npts);

    // Generate random points

    utilities::Timer timeInitRandPts(true);
    // Generate random points
    for(int i = 0; i < npts; i++)
    {
      pts[i] = get_rand_pt(meshBb);
    }
    SLIC_INFO(axom::fmt::format("Constructed {} random points in {} s.",
                                npts,
                                timeInitRandPts.elapsed()));

    primal::Point<int, 3> bins(25);

    utilities::Timer timeInitQuery(true);
    // Create PointInCell query object
    quest::PointInCell<mesh_tag, ExecSpace> query(&mesh, bins.data());
    SLIC_INFO(axom::fmt::format("Initialized point-in-cell query in {} s.",
                                timeInitQuery.elapsed()));

    IndexType* outCellIds = axom::allocate<IndexType>(npts);
    // Run query
    utilities::Timer timeRunQuery(true);
    query.locatePoints(npts, pts, outCellIds);
    double time = timeRunQuery.elapsed();
    SLIC_INFO(
      axom::fmt::format("Ran query on {} points in {} s -- rate: {} q/s",
                        npts,
                        time,
                        npts / time));

    axom::deallocate(pts);
  }
};

struct Arguments
{
  std::string file_name;
  int num_rand_pts {10000};
  ExecPolicy exec_space {ExecPolicy::CPU};

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    app
      .add_option("-f,--file", this->file_name, "specifies the input mesh file")
      ->check(axom::CLI::ExistingFile)
      ->required();

    app.add_option("-n,--num_pts",
                   this->num_rand_pts,
                   "the number of points to query");

    std::string pol_info =
      "Sets execution space of the SignedDistance query.\n";
    pol_info += "Set to \'seq\' to use sequential execution policy.";
#ifdef AXOM_USE_OPENMP
    pol_info += "\nSet to \'omp\' to use an OpenMP execution policy.";
#endif
#ifdef AXOM_USE_CUDA
    pol_info += "\nSet to \'gpu\' to use a GPU execution policy.";
#endif
    app.add_option("-e, --exec_space", this->exec_space, pol_info)
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(validExecPolicies));

    app.get_formatter()->column_width(40);

    // could throw an exception
    app.parse(argc, argv);

    slic::flushStreams();
  }
};

int main(int argc, char** argv)
{
  // STEP 1: initialize the logger
  initialize_logger();

  // STEP 2: parse command line arguments
  Arguments args;
  axom::CLI::App app {"Driver for signed distance query"};

  try
  {
    args.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    int retval = -1;
    retval = app.exit(e);
    finalize_logger();
    return retval;
  }
#ifdef AXOM_USE_CUDA
  if(args.exec_space == ExecPolicy::GPU)
  {
    using GPUExec = axom::CUDA_EXEC<256>;
    axom::setDefaultAllocator(axom::execution_space<GPUExec>::allocatorID());
  }
#endif

  // Open MFEM mesh
  utilities::Timer ctorMeshTimer(true);
  mfem::Mesh testMesh(args.file_name.c_str());
  SLIC_INFO(axom::fmt::format("Initialized MFEM mesh {} in {} s.",
                              args.file_name,
                              ctorMeshTimer.elapsed()));

  switch(args.exec_space)
  {
  case ExecPolicy::CPU:
    benchmark_point_in_cell<axom::SEQ_EXEC>(testMesh, args.num_rand_pts);
    break;
#ifdef AXOM_USE_OPENMP
  case ExecPolicy::OpenMP:
    benchmark_point_in_cell<axom::OMP_EXEC>(testMesh, args.num_rand_pts);
    break;
#endif
#ifdef AXOM_USE_CUDA
  case ExecPolicy::GPU:
    benchmark_point_in_cell<axom::CUDA_EXEC<256>>(testMesh, args.num_rand_pts);
    break;
#endif
  default:
    SLIC_ERROR("Unsupported execution space.");
    return 1;
  }
  finalize_logger();
  return 0;
}

//------------------------------------------------------------------------------
void initialize_logger()
{
  // initialize logger
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Info);

  // setup the logstreams
  std::string fmt = "";
  slic::LogStream* logStream = nullptr;

  fmt = "[<LEVEL>]: <MESSAGE>\n";
  logStream = new slic::GenericOutputStream(&std::cout, fmt);

  // register stream objects with the logger
  slic::addStreamToAllMsgLevels(logStream);
}

//------------------------------------------------------------------------------
void finalize_logger()
{
  slic::flushStreams();
  slic::finalize();
}
