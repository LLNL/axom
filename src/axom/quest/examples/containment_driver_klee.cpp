// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file containment_driver.cpp
 * \brief Basic demo of point containment acceleration structure over surfaces.
 */

// Axom includes
#include "axom/core.hpp"
#include "axom/primal.hpp"
#include "axom/sidre.hpp"
#include "axom/spin.hpp"
#include "axom/quest.hpp"
#include "axom/mint.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/klee.hpp"

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
namespace primal = axom::primal;
namespace quest = axom::quest;
namespace slic = axom::slic;
namespace klee = axom::klee;

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
enum VolFracSampling
{
  SAMPLE_AT_DOFS,
  SAMPLE_AT_QPTS
};

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

void FCT_project(mfem::DenseMatrix& M,
                 mfem::DenseMatrixInverse& M_inv,
                 mfem::Vector& m,
                 mfem::Vector& x,  // indicators
                 double y_min,
                 double y_max,
                 mfem::Vector& xy)  // indicators * rho
{
  // [IN]  - M, M_inv, m, x, y_min, y_max
  // [OUT] - xy

  using namespace mfem;

  const int s = M.Size();

  xy.SetSize(s);

  // Compute the high-order projection in xy
  M_inv.Mult(m, xy);

  // Q0 solutions can't be adjusted conservatively. It's what it is.
  if(xy.Size() == 1)
  {
    return;
  }

  // Compute the lumped mass matrix in ML
  Vector ML(s);
  M.GetRowSums(ML);

  //Ensure dot product is done on the CPU
  double dMLX(0);
  for(int i = 0; i < x.Size(); ++i)
  {
    dMLX += ML(i) * x(i);
  }

  const double y_avg = m.Sum() / dMLX;

#ifdef AXOM_DEBUG
  SLIC_WARNING_IF(!(y_min < y_avg + 1e-12 && y_avg < y_max + 1e-12),
                  fmt::format("Average ({}) is out of bounds [{},{}]: ",
                              y_avg,
                              y_min - 1e-12,
                              y_max + 1e-12));
#endif

  Vector z(s);
  Vector beta(s);
  Vector Mxy(s);
  M.Mult(xy, Mxy);
  for(int i = 0; i < s; i++)
  {
    // Some different options for beta:
    //beta(i) = 1.0;
    beta(i) = ML(i) * x(i);
    //beta(i) = ML(i)*(x(i) + 1e-14);
    //beta(i) = ML(i);
    //beta(i) = Mxy(i);

    // The low order flux correction
    z(i) = m(i) - ML(i) * x(i) * y_avg;
  }

  // Make beta_i sum to 1
  beta /= beta.Sum();

  DenseMatrix F(s);
  for(int i = 1; i < s; i++)
  {
    for(int j = 0; j < i; j++)
    {
      F(i, j) = M(i, j) * (xy(i) - xy(j)) + (beta(j) * z(i) - beta(i) * z(j));
    }
  }

  Vector gp(s), gm(s);
  gp = 0.0;
  gm = 0.0;
  for(int i = 1; i < s; i++)
  {
    for(int j = 0; j < i; j++)
    {
      double fij = F(i, j);
      if(fij >= 0.0)
      {
        gp(i) += fij;
        gm(j) -= fij;
      }
      else
      {
        gm(i) += fij;
        gp(j) -= fij;
      }
    }
  }

  for(int i = 0; i < s; i++)
  {
    xy(i) = x(i) * y_avg;
  }

  for(int i = 0; i < s; i++)
  {
    double mi = ML(i), xyLi = xy(i);
    double rp = std::max(mi * (x(i) * y_max - xyLi), 0.0);
    double rm = std::min(mi * (x(i) * y_min - xyLi), 0.0);
    double sp = gp(i), sm = gm(i);

    gp(i) = (rp < sp) ? rp / sp : 1.0;
    gm(i) = (rm > sm) ? rm / sm : 1.0;
  }

  for(int i = 1; i < s; i++)
  {
    for(int j = 0; j < i; j++)
    {
      double fij = F(i, j), aij;

      if(fij >= 0.0)
      {
        aij = std::min(gp(i), gm(j));
      }
      else
      {
        aij = std::min(gm(i), gp(j));
      }

      fij *= aij;
      xy(i) += fij / ML(i);
      xy(j) -= fij / ML(j);
    }
  }
}

void generate_volume_fractions_baseline(mfem::DataCollection* dc,
                                        mfem::QuadratureFunction* inout,
                                        const std::string& name,  // vol_frac
                                        int /*order*/)
{
  const int order = inout->GetSpace()->GetElementIntRule(0).GetOrder();

  mfem::Mesh* mesh = dc->GetMesh();
  const int dim = mesh->Dimension();
  const int NE = mesh->GetNE();

  std::cout << fmt::format("Mesh has dim {} and {} elements", dim, NE)
            << std::endl;

  mfem::L2_FECollection* fec =
    new mfem::L2_FECollection(order, dim, mfem::BasisType::Positive);
  mfem::FiniteElementSpace* fes = new mfem::FiniteElementSpace(mesh, fec);
  mfem::GridFunction* volFrac = new mfem::GridFunction(fes);
  volFrac->MakeOwner(fec);
  dc->RegisterField(name, volFrac);

  (*volFrac) = (*inout);
}

void generate_volume_fractions(mfem::DataCollection* dc,
                               mfem::QuadratureFunction* inout,
                               const std::string& name,  // vol_frac
                               int order)
{
  const int sampleOrder = inout->GetSpace()->GetElementIntRule(0).GetOrder();
  const int sampleNQ = inout->GetSpace()->GetElementIntRule(0).GetNPoints();
  const int sampleSZ = inout->GetSpace()->GetSize();

  SLIC_INFO(fmt::format(std::locale("en_US.UTF-8"),
                        "In generate_volume_fractions: sample order {} | "
                        "sample num qpts {} |  total samples {:L}",
                        sampleOrder,
                        sampleNQ,
                        sampleSZ));

  mfem::Mesh* mesh = dc->GetMesh();
  const int dim = mesh->Dimension();
  const int NE = mesh->GetNE();

  SLIC_INFO(fmt::format(std::locale("en_US.UTF-8"),
                        "Mesh has dim {} and {:L} elements",
                        dim,
                        NE));

  mfem::L2_FECollection* fec =
    new mfem::L2_FECollection(order, dim, mfem::BasisType::Positive);
  mfem::FiniteElementSpace* fes = new mfem::FiniteElementSpace(mesh, fec);
  mfem::GridFunction* volFrac = new mfem::GridFunction(fes);
  volFrac->MakeOwner(fec);
  dc->RegisterField(name, volFrac);

  axom::utilities::Timer timer(true);
  {
    mfem::MassIntegrator mass_integrator;  // use the default for exact integration; lower for approximate

    mfem::QuadratureFunctionCoefficient qfc(*inout);
    mfem::DomainLFIntegrator rhs(qfc);

    // assume all elts are the same
    const auto& ir = inout->GetSpace()->GetElementIntRule(0);
    rhs.SetIntRule(&ir);

    mfem::DenseMatrix m;
    mfem::DenseMatrixInverse mInv;
    mfem::Vector b, x;
    mfem::Array<int> dofs;
    mfem::Vector one;
    const double minY = 0.;
    const double maxY = 1.;

    for(int i = 0; i < NE; ++i)
    {
      auto* T = mesh->GetElementTransformation(i);
      auto* el = fes->GetFE(i);

      mass_integrator.AssembleElementMatrix(*el, *T, m);
      rhs.AssembleRHSElementVect(*el, *T, b);
      mInv.Factor(m);

      // Use FCT limiting -- similar to Remap
      // Limit the function to be between 0 and 1
      // Q: Is there a better way limiting algorithm for this?
      if(one.Size() != b.Size())
      {
        one.SetSize(b.Size());
        one = 1.0;
      }
      FCT_project(m, mInv, b, one, minY, maxY, x);

      fes->GetElementDofs(i, dofs);
      volFrac->SetSubVector(dofs, x);
    }
  }
  timer.stop();
  SLIC_INFO(fmt::format(std::locale("en_US.UTF-8"),
                        "\t Generating volume fractions took {:.3f} seconds (@ "
                        "{:L} dofs processed per second)",
                        timer.elapsed(),
                        static_cast<int>(fes->GetNDofs() / timer.elapsed())));

  // // Alt strategy for Q0: take a simple average of samples
  // const int nq = sampleNQ;
  // double sc = 1. / nq;
  // for(int i = 0; i < NE; ++i)
  // {
  //   int sum = 0;
  //   for(int k = 0; k < nq; ++k)
  //   {
  //     sum += (*inout)(i * nq + k);
  //   }
  //   (*volFrac)[i] = sum * sc;
  // }
}

/**
 * Compute volume fractions function for shape on a grid of resolution \a gridRes
 * in region defined by bounding box \a queryBounds
 */
void computeVolumeFractionsBaseline(int id,
                                    const Octree3D& inOutOctree,
                                    mfem::DataCollection* dc,
                                    int AXOM_NOT_USED(sampleRes),
                                    int outputOrder)
{
  // Step 1 -- generate a QField w/ the spatial coordinates
  mfem::Mesh* mesh = dc->GetMesh();
  const int NE = mesh->GetNE();
  const int dim = mesh->Dimension();

  if(NE < 1)
  {
    SLIC_WARNING("Mesh has no elements!");
    return;
  }

  mfem::L2_FECollection* coll =
    new mfem::L2_FECollection(outputOrder, dim, mfem::BasisType::Positive);
  mfem::FiniteElementSpace* fes = new mfem::FiniteElementSpace(mesh, coll);
  mfem::GridFunction* volFrac = new mfem::GridFunction(fes);
  volFrac->MakeOwner(coll);
  dc->RegisterField(fmt::format("vol_frac_{:03}", id), volFrac);

  auto* fe = fes->GetFE(0);
  auto& ir = fe->GetNodes();

  // Assume all elements have the same integration rule
  const int nq = ir.GetNPoints();
  const auto* geomFactors =
    mesh->GetGeometricFactors(ir, mfem::GeometricFactors::COORDINATES);

  mfem::DenseTensor pos_coef(dim, nq, NE);

  // Rearrange positions into quadrature function
  {
    for(int i = 0; i < NE; ++i)
    {
      for(int j = 0; j < dim; ++j)
      {
        for(int k = 0; k < nq; ++k)
        {
          pos_coef(j, k, i) = geomFactors->X((i * nq * dim) + (j * nq) + k);
        }
      }
    }
  }

  // Step 2 -- sample the in/out field at each point -- store directly in volFrac grid function
  mfem::Vector res(nq);
  mfem::Array<int> dofs;
  for(int i = 0; i < NE; ++i)
  {
    mfem::DenseMatrix& m = pos_coef(i);
    for(int p = 0; p < nq; ++p)
    {
      const axom::primal::Point3D pt(m.GetColumn(p), dim);
      const bool in = inOutOctree.within(pt);
      res(p) = in ? 1. : 0.;
    }

    fes->GetElementDofs(i, dofs);
    volFrac->SetSubVector(dofs, res);
  }
}

/**
 * Compute volume fractions function for shape on a grid of resolution \a gridRes
 * in region defined by bounding box \a queryBounds
 */
void computeVolumeFractions(int id,
                            const Octree3D& inOutOctree,
                            mfem::DataCollection* dc,
                            int sampleRes,
                            int outputOrder)
{
  // Step 1 -- generate a QField w/ the spatial coordinates
  mfem::Mesh* mesh = dc->GetMesh();
  const int NE = mesh->GetNE();
  const int dim = mesh->Dimension();

  if(NE < 1)
  {
    SLIC_WARNING("Mesh has no elements!");
    return;
  }

  // convert requested samples into a compatible polynomial order
  // that will use that many samples: 2n-1 and 2n-2 will work
  // NOTE: Might be different for simplices
  const int sampleOrder = 2 * sampleRes - 1;
  mfem::QuadratureSpace sp(mesh, sampleOrder);

  // TODO: Should the samples be along a uniform grid
  //       instead of Guassian quadrature?
  //       This would need quadrature weights for the uniform
  //       samples -- Newton-Cotes ?
  //       With uniform points, we could do HO polynomial fitting
  //       Using 0s and 1s is non-oscillatory in Bernstein basis

  // Assume all elements have the same integration rule
  const auto& ir = sp.GetElementIntRule(0);
  const int nq = ir.GetNPoints();
  const auto* geomFactors =
    mesh->GetGeometricFactors(ir, mfem::GeometricFactors::COORDINATES);

  mfem::QuadratureFunction pos_coef(&sp, dim);

  // Rearrange positions into quadrature function
  {
    for(int i = 0; i < NE; ++i)
    {
      for(int j = 0; j < dim; ++j)
      {
        for(int k = 0; k < nq; ++k)
        {
          //X has dims nqpts x sdim x ne
          pos_coef((i * nq * dim) + (k * dim) + j) =
            geomFactors->X((i * nq * dim) + (j * nq) + k);
        }
      }
    }
  }

  const int vdim = 1;
  mfem::QuadratureFunction inout(&sp, vdim);

  // Step 2 -- sample the in/out field at each point -- store in QField 'inout'
  mfem::DenseMatrix m;
  mfem::Vector res;

  axom::utilities::Timer timer(true);
  for(int i = 0; i < NE; ++i)
  {
    pos_coef.GetElementValues(i, m);
    inout.GetElementValues(i, res);

    for(int p = 0; p < nq; ++p)
    {
      const axom::primal::Point3D pt(m.GetColumn(p), dim);
      const bool in = inOutOctree.within(pt);
      res(p) = in ? 1. : 0.;

      // SLIC_INFO(fmt::format("[{},{}] Pt: {}, In: {}", i,p,pt, (in? "yes" : "no") ));
    }
  }
  timer.stop();
  SLIC_INFO(fmt::format(
    std::locale("en_US.UTF-8"),
    "\t Sampling inout field took {} seconds (@ {:L} queries per second)",
    timer.elapsed(),
    static_cast<int>((NE * nq) / timer.elapsed())));

  // Step 3 -- project QField onto volume fractions field
  generate_volume_fractions(dc,
                            &inout,
                            fmt::format("vol_frac_{:03}", id),
                            outputOrder);
}

/** Struct to parse and store the input parameters */
struct Input
{
public:
  std::string shapeFile;
  klee::ShapeSet shapeSet;

  int maxQueryLevel {7};
  int quadratureOrder {5};
  int outputOrder {2};
  VolFracSampling vfSampling {SAMPLE_AT_QPTS};

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
    auto* minbb =
      app.add_option("--min", queryBoxMins, "Min bounds for query box (x,y,z)")
        ->expected(3);
    auto* maxbb =
      app.add_option("--max", queryBoxMaxs, "Max bounds for query box (x,y,z)")
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
          reader.setFileName(geom.getPath());
          reader.read();

          UMesh surface_mesh(3, mint::TRIANGLE);
          reader.getMesh(&surface_mesh);

          // NOTE: We're not yet applying transformations!

          GeometricBoundingBox bb = compute_bounds(&surface_mesh);
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

  auto bbox = params.problemBoundingBox();
  auto range = bbox.range();
  auto low = bbox.getMin();
  mfem::Mesh* mesh = new mfem::Mesh(res,
                                    res,
                                    res,
                                    mfem::Element::HEXAHEDRON,
                                    false,
                                    range[0],
                                    range[1],
                                    range[2]);

  // Offset to the mesh to lie w/in the bounding box
  for(int i = 0; i < mesh->GetNV(); ++i)
  {
    double* v = mesh->GetVertex(i);
    v[0] += low[0];
    v[1] += low[1];
    v[2] += low[2];
  }

  dc->SetMesh(MPI_COMM_WORLD, mesh);
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
  std::fstream file(params.shapeFile);
  params.shapeSet = klee::readShapeSet(file);

  params.initializeProblemBoundingBox();

  auto bbox = params.problemBoundingBox();
  SLIC_INFO("Mesh bounding box: " << bbox);

  // Add the mesh to an mfem data collection
  axom::sidre::MFEMSidreDataCollection dc("shaping", nullptr, true);

  initializeMesh(params, &dc);

  // for each shape, generate a volume fraction field
  int id = 0;
  for(const auto& s : params.shapeSet.getShapes())
  {
    std::string stlPath = s.getGeometry().getPath();
    std::string outMsg = fmt::format(" Loading mesh {} ", id);
    SLIC_INFO(fmt::format("{:*^80}", outMsg));
    SLIC_INFO("Reading file: " << stlPath << "...");

    quest::STLReader reader;
    reader.setFileName(stlPath);
    reader.read();

    // Create surface mesh
    mint::Mesh* surface_mesh = new UMesh(3, mint::TRIANGLE);
    reader.getMesh(static_cast<UMesh*>(surface_mesh));

    SLIC_INFO("Mesh has " << surface_mesh->getNumberOfNodes() << " nodes and "
                          << surface_mesh->getNumberOfCells() << " cells.");

    // Compute mesh bounding box and log some stats about the surface
    GeometricBoundingBox meshBB = compute_bounds(surface_mesh);
    SLIC_INFO("Mesh bounding box: " << meshBB);

    // Create octree over mesh's bounding box
    SLIC_INFO(fmt::format("{:*^80}", " Generating the octree "));
    Octree3D octree(meshBB, surface_mesh);
    octree.generateIndex();

    {
      const int nVerts = surface_mesh->getNumberOfNodes();
      const int nCells = surface_mesh->getNumberOfCells();

      SLIC_INFO(fmt::format(
        "After welding, surface mesh has {} vertices  and {} triangles.",
        nVerts,
        nCells));
      mint::write_vtk(surface_mesh, fmt::format("meldedTriMesh_{:03}.vtk", id));
    }

    SLIC_INFO(fmt::format("{:*^80}", " Querying the octree "));

    // Generate volume fractions of the desired output order and sampling resolution
    const int sampleOrder = params.quadratureOrder;
    const int outputOrder = params.outputOrder;
    switch(params.vfSampling)
    {
    case SAMPLE_AT_QPTS:
      computeVolumeFractions(id, octree, &dc, sampleOrder, outputOrder);
      break;
    case SAMPLE_AT_DOFS:
      computeVolumeFractionsBaseline(id, octree, &dc, sampleOrder, outputOrder);
      break;
    }

    delete surface_mesh;
    surface_mesh = nullptr;

    ++id;
  }

  // Save meshes and fields
  dc.Save();

  MPI_Finalize();

  return 0;
}
