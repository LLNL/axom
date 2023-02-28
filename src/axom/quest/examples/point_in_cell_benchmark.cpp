// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/quest/PointInCell.hpp"

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

namespace
{
#ifdef AXOM_USE_GPU
  #ifdef AXOM_USE_HIP
using GPUExec = axom::HIP_EXEC<256>;
  #else
using GPUExec = axom::CUDA_EXEC<256>;
  #endif
#endif

enum class ExecPolicy
{
  CPU,
  OpenMP,
  GPU
};

const std::map<std::string, ExecPolicy> validExecPolicies {
  {"seq", ExecPolicy::CPU},
#ifdef AXOM_USE_RAJA
  #ifdef AXOM_USE_OPENMP
  {"omp", ExecPolicy::OpenMP},
  #endif
  #ifdef AXOM_USE_GPU
  {"gpu", ExecPolicy::GPU}
  #endif
#endif
};

}  // namespace

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

template <typename ExecSpace, int DIM>
void benchmark_point_in_cell(mfem::Mesh& mesh, int npts, int nbins)
{
  using PointType = primal::Point<double, DIM>;
  using BoxType = primal::BoundingBox<double, DIM>;
  using mesh_tag = axom::quest::quest_point_in_cell_mfem_tag;
  using IndexType = typename quest::PointInCellTraits<mesh_tag>::IndexType;

  BoxType meshBb;
  {
    mfem::Vector meshMin, meshMax;
    mesh.GetBoundingBox(meshMin, meshMax);
    meshBb.addPoint(PointType(meshMin.GetData()));
    meshBb.addPoint(PointType(meshMax.GetData()));
  }

  axom::Array<PointType> pts(npts, npts, axom::getDefaultAllocatorID());

  utilities::Timer timeInitRandPts(true);
  // Generate random points
  for(int i = 0; i < npts; i++)
  {
    pts[i] = get_rand_pt(meshBb);
  }
  SLIC_INFO(axom::fmt::format("Constructed {} random points in {} s.",
                              npts,
                              timeInitRandPts.elapsed()));

  primal::Point<int, DIM> bins(nbins);

  utilities::Timer timeInitQuery(true);
  quest::PointInCell<mesh_tag, ExecSpace> query(&mesh, bins.data());
  SLIC_INFO(axom::fmt::format("Initialized point-in-cell query in {} s.",
                              timeInitQuery.elapsed()));

  axom::Array<IndexType> outCellIds(npts, npts, axom::getDefaultAllocatorID());

  // Run query
  utilities::Timer timeRunQuery(true);
  query.locatePoints(pts.view(), outCellIds.data());
  double time = timeRunQuery.elapsed();
  SLIC_INFO(axom::fmt::format("Ran query on {} points in {} s -- rate: {} q/s",
                              npts,
                              time,
                              npts / time));
}

template <typename ExecSpace>
void benchmark_point_in_cell(mfem::Mesh& mesh, int npts, int nbins)
{
  switch(mesh.SpaceDimension())
  {
  case 2:
    benchmark_point_in_cell<ExecSpace, 2>(mesh, npts, nbins);
    break;
  case 3:
    benchmark_point_in_cell<ExecSpace, 3>(mesh, npts, nbins);
    break;
  default:
    SLIC_ERROR("Unsupported dimension (" << mesh.SpaceDimension() << ")");
    break;
  }
}

struct Arguments
{
  std::string file_name;
  int num_rand_pts {10000};
  int num_bins {25};
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

    app.add_option("-b,--num-bins",
                   this->num_bins,
                   "the number of bins to construct for each dimension");

    std::string pol_info = "Sets execution space of the PointInCell query.\n";
    pol_info += "Set to 'seq' to use sequential execution policy.";
#ifdef AXOM_USE_RAJA
  #ifdef AXOM_USE_OPENMP
    pol_info += "\nSet to 'omp' to use an OpenMP execution policy.";
  #endif
  #ifdef AXOM_USE_GPU
    pol_info += "\nSet to 'gpu' to use a GPU execution policy.";
  #endif
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
  slic::SimpleLogger logger(slic::message::Info);

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
    return retval;
  }

#ifdef AXOM_USE_GPU
  if(args.exec_space == ExecPolicy::GPU)
  {
    axom::setDefaultAllocator(axom::execution_space<GPUExec>::allocatorID());
  }
#endif

  // Open MFEM mesh
  utilities::Timer ctorMeshTimer(true);
  mfem::Mesh testMesh(args.file_name.c_str());
  testMesh.EnsureNodes();
  SLIC_INFO(
    axom::fmt::format("Initialized MFEM mesh '{}' in {} s. \n"
                      "\tThe mesh has {} {}d elements and its order is {}.",
                      args.file_name,
                      ctorMeshTimer.elapsed(),
                      testMesh.GetNE(),
                      testMesh.Dimension(),
                      testMesh.GetNodalFESpace()->FEColl()->GetOrder()));

  switch(args.exec_space)
  {
  case ExecPolicy::CPU:
    benchmark_point_in_cell<axom::SEQ_EXEC>(testMesh,
                                            args.num_rand_pts,
                                            args.num_bins);
    break;
#ifdef AXOM_USE_RAJA
  #ifdef AXOM_USE_OPENMP
  case ExecPolicy::OpenMP:
    benchmark_point_in_cell<axom::OMP_EXEC>(testMesh,
                                            args.num_rand_pts,
                                            args.num_bins);
    break;
  #endif
  #ifdef AXOM_USE_GPU
  case ExecPolicy::GPU:
    benchmark_point_in_cell<GPUExec>(testMesh, args.num_rand_pts, args.num_bins);
    break;
  #endif
#endif
  default:
    SLIC_ERROR("Unsupported execution space.");
    return 1;
  }

  return 0;
}
