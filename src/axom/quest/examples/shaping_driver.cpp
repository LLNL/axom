// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file shaping_driver.cpp
 * \brief Driver for shaping material volume fractions onto a simulation mesh
 */

// Axom includes
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/sidre.hpp"
#include "axom/primal.hpp"
#include "axom/spin.hpp"
#include "axom/klee.hpp"
#include "axom/mint.hpp"
#include "axom/quest.hpp"

#include "fmt/fmt.hpp"
#include "fmt/locale.h"
#include "CLI11/CLI11.hpp"

#include "mfem.hpp"

// C/C++ includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <fstream>
#include <iomanip>  // for setprecision()

#include "mpi.h"

// NOTE: Axom must be configured with AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION compiler define for the klee containment driver

namespace mint = axom::mint;
namespace quest = axom::quest;
namespace slic = axom::slic;
namespace klee = axom::klee;

using UMesh = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;

using Octree3D = quest::InOutOctree<3>;
using GeometricBoundingBox = Octree3D::GeometricBoundingBox;
using SpacePt = Octree3D::SpacePt;

using QFunctionCollection = mfem::NamedFieldsMap<mfem::QuadratureFunction>;
using DenseTensorCollection = mfem::NamedFieldsMap<mfem::DenseTensor>;

using VolFracSampling = quest::shaping::VolFracSampling;

//------------------------------------------------------------------------------

/** Struct to parse and store the input parameters */
struct Input
{
public:
  std::string shapeFile;
  klee::ShapeSet shapeSet;

  int maxQueryLevel {5};
  int quadratureOrder {5};
  int outputOrder {2};
  VolFracSampling vfSampling {VolFracSampling::SAMPLE_AT_QPTS};

  std::vector<double> queryBoxMins;
  std::vector<double> queryBoxMaxs;

private:
  bool m_verboseOutput {false};
  bool m_hasUserQueryBox {false};
  GeometricBoundingBox m_problemBoundingBox;

public:
  Input()
  {
// increase default for max query resolution in release builds
#ifndef AXOM_DEBUG
    maxQueryLevel += 2;
#endif
  }

  bool isVerbose() const { return m_verboseOutput; }

  bool hasUserQueryBox() const { return m_hasUserQueryBox; }

  GeometricBoundingBox problemBoundingBox() { return m_problemBoundingBox; }

  // HACK to allow 2D querying.
  // For now, use a 2D mesh if the supplied 3rd parameter in the bounding boxes are 0
  // TODO: Generalize to better support for 2D querying by allowing the user
  //       to more easily supply the 2D bounding box
  //       Alternatively, we will be using the "slice" operation, when supported by klee/quest
  int meshDimension() const
  {
    if(hasUserQueryBox())
    {
      return (queryBoxMins[2] == 0. && queryBoxMaxs[2] == 0) ? 2 : 3;
    }

    return 3;
  }

  void parse(int argc, char** argv, CLI::App& app)
  {
    app.add_option("shapeFile", shapeFile, "Path to input shape file")
      ->check(CLI::ExistingFile)
      ->required();

    app
      .add_flag("-v,--verbose", m_verboseOutput, "Enable/disable verbose output")
      ->capture_default_str();

    app
      .add_option(
        "-l,--levels",
        maxQueryLevel,
        "Max query resolution. \n"
        "Will query uniform grids at levels 1 through the provided level")
      ->capture_default_str()
      ->check(CLI::PositiveNumber);

    app
      .add_option("-o,--order", outputOrder, "order of the output grid function")
      ->capture_default_str()
      ->check(CLI::PositiveNumber);

    app
      .add_option("-q,--quadrature-order",
                  quadratureOrder,
                  "Quadrature order for sampling the inout field. \n"
                  "Determines number of samples per element in determining "
                  "volume fraction field")
      ->capture_default_str()
      ->check(CLI::PositiveNumber);

    std::map<std::string, VolFracSampling> vfsamplingMap {
      {"qpts", VolFracSampling::SAMPLE_AT_QPTS},
      {"dofs", VolFracSampling::SAMPLE_AT_DOFS}};
    app
      .add_option("-s,--sampling-type",
                  vfSampling,
                  "Sampling strategy. \n"
                  "Sampling either at quadrature points or collocated with "
                  "degrees of freedom")
      ->capture_default_str()
      ->transform(CLI::CheckedTransformer(vfsamplingMap, CLI::ignore_case));

    // Optional bounding box for query region
    // Note: We're using CLI11@1.8; CLI11@1.9 allows expected(min,max) to simplify 2D and 3D processing
    auto* minbb =
      app.add_option("--min", queryBoxMins, "Min bounds for query box (x,y,z)")
        ->expected(3);
    auto* maxbb =
      app
        .add_option("--max",
                    queryBoxMaxs,
                    "Max bounds for query box (x,y,z). \n"
                    "Query mesh will be 2D if min[2] == max[2] == 0")
        ->expected(3);
    minbb->needs(maxbb);
    maxbb->needs(minbb);

    app.get_formatter()->column_width(35);

    // could throw an exception
    app.parse(argc, argv);

    slic::setLoggingMsgLevel(m_verboseOutput ? slic::message::Debug
                                             : slic::message::Info);

    m_hasUserQueryBox = app.count("--min") == 3 && app.count("--max") == 3;
  }

  void initializeProblemBoundingBox()
  {
    m_problemBoundingBox.clear();

    if(this->hasUserQueryBox())
    {
      m_problemBoundingBox.addPoint(SpacePt(queryBoxMins.data(), 3));
      m_problemBoundingBox.addPoint(SpacePt(queryBoxMaxs.data(), 3));
    }
    else
    {
      // Hack -- If not supplied by user, create problem bounding box by reading all STL meshes
      for(const auto& s : shapeSet.getShapes())
      {
        // SLIC_INFO(fmt::format("Reading shape '{}' of material '{}'",
        //                       s.getName(),
        //                       s.getMaterial()));

        // Initial assumption is that all shapes are provided as STL files
        auto& geom = s.getGeometry();
        if(geom.getFormat() == "stl")
        {
          quest::STLReader reader;
          reader.setFileName(shapeSet.resolvePath(geom.getPath()));
          reader.read();

          UMesh surface_mesh(3, mint::TRIANGLE);
          reader.getMesh(&surface_mesh);

          // NOTE: We're not yet applying transformations!

          GeometricBoundingBox bb = quest::shaping::compute_bounds(surface_mesh);
          m_problemBoundingBox.addBox(bb);
        }
      }
    }
  }
};

void initializeMesh(Input& params, axom::sidre::MFEMSidreDataCollection* dc)
{
  // Create a background mesh
  // Generate an mfem Cartesian mesh, scaled to the bounding box range
  const int res = 1 << params.maxQueryLevel;

  const int dim = params.meshDimension();

  auto bbox = params.problemBoundingBox();
  auto range = bbox.range();
  auto low = bbox.getMin();
  mfem::Mesh* mesh = nullptr;

  switch(dim)
  {
  case 2:
    mesh = new mfem::Mesh(res,
                          res,
                          mfem::Element::QUADRILATERAL,
                          false,
                          range[0],
                          range[1]);
    break;
  case 3:
    mesh = new mfem::Mesh(res,
                          res,
                          res,
                          mfem::Element::HEXAHEDRON,
                          false,
                          range[0],
                          range[1],
                          range[2]);
    break;
  default:
    SLIC_ERROR("Only 2D and 3D meshes are currently supported.");
    break;
  }

  // Offset to the mesh to lie w/in the bounding box
  for(int i = 0; i < mesh->GetNV(); ++i)
  {
    double* v = mesh->GetVertex(i);
    for(int d = 0; d < dim; ++d)
    {
      v[d] += low[d];
    }
  }

  mesh->EnsureNodes();
  dc->SetMeshNodesName("positions");

#ifdef MFEM_USE_MPI
  dc->SetMesh(MPI_COMM_WORLD, mesh);
#endif
}

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

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

  // Load shape file and extract info
  params.shapeSet = klee::readShapeSet(params.shapeFile);
  const klee::Dimensions shapeDim = params.shapeSet.getDimensions();

  params.initializeProblemBoundingBox();

  auto bbox = params.problemBoundingBox();
  SLIC_INFO("Mesh bounding box: " << bbox);

  // Add the mesh to an mfem data collection
  quest::MFEMShaping shaper(params.shapeSet);

  initializeMesh(params, shaper.getDC());
  const int queryDim = params.meshDimension();
  const int sampleOrder = params.quadratureOrder;
  const int outputOrder = params.outputOrder;

  // Sample the InOut quadrature field for each shape using an InOut octree
  // Assumptions: Each shape has a unique name
  SLIC_INFO(fmt::format("{:=^80}", "Sampling InOut fields for shapes"));
  for(const auto& s : params.shapeSet.getShapes())
  {
    shaper.loadShape(s);

    auto* surface_mesh = shaper.getSurfaceMesh();

    // Compute mesh bounding box and log some stats about the surface
    shaper.prepareShapeQuery(shapeDim, s);

    shaper.runShapeQuery(params.vfSampling, sampleOrder, outputOrder);

    shaper.finalizeShapeQuery();
  }

  auto& inoutQFuncs = shaper.getInoutQFuncs();

  // Apply replacement rules to the quadrature points
  // Assumptions: The replacement rules have been validated, yielding a valid DAG for the replacements
  SLIC_INFO(
    fmt::format("{:=^80}", "Applying replacement rules over the shapes"));

  // generate a map from materials to shape names
  std::map<std::string, std::vector<std::string>> materialsToShapes;
  {
    for(const auto& s : params.shapeSet.getShapes())
    {
      materialsToShapes[s.getMaterial()].push_back(s.getName());
    }

    // generate a map from materials to set of materials that it is replaced by
    std::map<std::string, std::set<std::string>> replaced_by;
    for(const auto& s : params.shapeSet.getShapes())
    {
      auto& material_s = s.getMaterial();
      for(const auto& t : params.shapeSet.getShapes())
      {
        auto& material_t = t.getMaterial();
        if(material_s != material_t && s.replaces(material_t))
        {
          replaced_by[material_t].insert(material_s);
        }
      }
    }

    SLIC_INFO("Replacement rules:");
    for(const auto& s : params.shapeSet.getShapes())
    {
      SLIC_INFO(fmt::format("Shape '{}' of material {} is replaced by: [{}]",
                            s.getName(),
                            s.getMaterial(),
                            fmt::join(replaced_by[s.getMaterial()], " ")));
    }

    // Merge all shapes of a given material into a single material QFunc
    for(const auto& kv : materialsToShapes)
    {
      auto& mat = kv.first;
      quest::shaping::mergeQFuncs(mat, kv.second, inoutQFuncs);
    }

    // Merge all shapes of a given material into a single material QFunc
    for(const auto& kv : replaced_by)
    {
      auto& mat = kv.first;
      quest::shaping::replaceMaterials(mat, kv.second, inoutQFuncs);
    }
  }
  // Generate the volume fractions from the InOut quadrature fields
  SLIC_INFO(
    fmt::format("{:=^80}", "Generating volume fraction fields for materials"));
  for(const auto& kv : materialsToShapes)
  {
    const std::string shapeName = kv.first;
    const int outputOrder = params.outputOrder;

    SLIC_INFO(fmt::format("Generating volume fraction fields for '{}' shape",
                          shapeName));

    switch(params.vfSampling)
    {
    case VolFracSampling::SAMPLE_AT_QPTS:
      quest::shaping::computeVolumeFractions(shapeName,
                                             shaper.getDC(),
                                             inoutQFuncs,
                                             outputOrder);
      break;
    case VolFracSampling::SAMPLE_AT_DOFS:
      /* no-op for now */
      break;
    }
  }

// Save meshes and fields
#ifdef MFEM_USE_MPI
  shaper.getDC()->Save();
#endif

  MPI_Finalize();

  return 0;
}
