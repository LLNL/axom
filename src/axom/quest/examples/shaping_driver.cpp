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

using QFunctionCollection = mfem::NamedFieldsMap<mfem::QuadratureFunction>;
using DenseTensorCollection = mfem::NamedFieldsMap<mfem::DenseTensor>;

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

/**
 * Compute volume fractions function for shape on a grid of resolution \a gridRes
 * in region defined by bounding box \a queryBounds
 */
void computeVolumeFractionsBaseline(const std::string shapeName,
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
  auto volFracName = fmt::format("vol_frac_{}", shapeName);
  dc->RegisterField(volFracName, volFrac);

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

void generatePositionsQFunction(mfem::Mesh* mesh,
                                QFunctionCollection& inoutQFuncs,
                                int sampleRes)
{
  SLIC_ASSERT(mesh != nullptr);
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
  mfem::QuadratureSpace* sp = new mfem::QuadratureSpace(mesh, sampleOrder);

  // TODO: Should the samples be along a uniform grid
  //       instead of Guassian quadrature?
  //       This would need quadrature weights for the uniform
  //       samples -- Newton-Cotes ?
  //       With uniform points, we could do HO polynomial fitting
  //       Using 0s and 1s is non-oscillatory in Bernstein basis

  // Assume all elements have the same integration rule
  const auto& ir = sp->GetElementIntRule(0);
  const int nq = ir.GetNPoints();
  const auto* geomFactors =
    mesh->GetGeometricFactors(ir, mfem::GeometricFactors::COORDINATES);

  mfem::QuadratureFunction* pos_coef = new mfem::QuadratureFunction(sp, dim);
  pos_coef->SetOwnsSpace(true);

  // Rearrange positions into quadrature function
  {
    for(int i = 0; i < NE; ++i)
    {
      for(int j = 0; j < dim; ++j)
      {
        for(int k = 0; k < nq; ++k)
        {
          //X has dims nqpts x sdim x ne
          (*pos_coef)((i * nq * dim) + (k * dim) + j) =
            geomFactors->X((i * nq * dim) + (j * nq) + k);
        }
      }
    }
  }

  // register positions with the QFunction collection, which wil handle its deletion
  inoutQFuncs.Register("positions", pos_coef, true);
}

void sampleInOutField(const std::string& shapeName,
                      const Octree3D& inOutOctree,
                      mfem::DataCollection* dc,
                      QFunctionCollection& inoutQFuncs,
                      int sampleRes)
{
  auto* mesh = dc->GetMesh();
  SLIC_ASSERT(mesh != nullptr);
  const int NE = mesh->GetNE();
  const int dim = mesh->Dimension();

  // Generate a Quadrature Function with the geometric positions, if not already available
  if(!inoutQFuncs.Has("positions"))
  {
    generatePositionsQFunction(mesh, inoutQFuncs, sampleRes);
  }

  // Access the positions QFunc and associated QuadratureSpace
  mfem::QuadratureFunction* pos_coef = inoutQFuncs.Get("positions");
  mfem::QuadratureSpace* sp = pos_coef->GetSpace();
  const int nq = sp->GetElementIntRule(0).GetNPoints();

  // Sample the in/out field at each point
  // store in QField which we register with the QFunc collection
  const std::string inoutName = fmt::format("inout_{}", shapeName);
  const int vdim = 1;
  auto* inout = new mfem::QuadratureFunction(sp, vdim);
  inoutQFuncs.Register(inoutName, inout, true);

  mfem::DenseMatrix m;
  mfem::Vector res;

  axom::utilities::Timer timer(true);
  for(int i = 0; i < NE; ++i)
  {
    pos_coef->GetElementValues(i, m);
    inout->GetElementValues(i, res);

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
    "\t Sampling inout field '{}' took {} seconds (@ {:L} queries per second)",
    inoutName,
    timer.elapsed(),
    static_cast<int>((NE * nq) / timer.elapsed())));
}

/**
 * Utility function to take the union of inouts on all shapes for a given material
 * 
 * Note: Registers the new QFunction with \a inoutQFuncs
 */
void mergeQFuncs(const std::string& material,
                 const std::vector<std::string>& shapes,
                 QFunctionCollection& inoutQFuncs)
{
  if(shapes.empty())
  {
    return;
  }

  std::vector<std::pair<std::string, mfem::QuadratureFunction*>> shapeQFuncs;
  for(auto& s : shapes)
  {
    std::string name = fmt::format("inout_{}", s);
    auto* q = inoutQFuncs.Get(name);
    SLIC_ASSERT(q != nullptr);
    shapeQFuncs.push_back({name, q});
  }

  // initialize material Q function from first shape
  auto* firstQFunc = shapeQFuncs[0].second;
  auto* mat_inout = new mfem::QuadratureFunction(*firstQFunc);
  inoutQFuncs.Register(fmt::format("mat_inout_{}", material), mat_inout, true);

  const int SZ = mat_inout->Size();
  double* mdata = mat_inout->GetData();

  // add in Q functions from other shapes
  const int nShapes = shapes.size();
  for(int i = 1; i < nShapes; ++i)
  {
    auto& pair = shapeQFuncs[i];
    double* sdata = pair.second->GetData();
    for(int j = 0; j < SZ; ++j)
    {
      mdata[j] = (mdata[j] + sdata[j] > 0) ? 1 : 0;
    }
  }
}

/**
 * Utility function to zero out inout quadrature points when a material is replaced
 */
void replaceMaterials(const std::string& material,
                      const std::set<std::string>& replacements,
                      QFunctionCollection& inoutQFuncs)
{
  if(replacements.empty())
  {
    return;
  }

  auto* matQFunc = inoutQFuncs.Get(fmt::format("mat_inout_{}", material));
  SLIC_ASSERT(matQFunc != nullptr);

  const int SZ = matQFunc->Size();
  double* mData = matQFunc->GetData();

  for(auto& m : replacements)
  {
    std::string name = fmt::format("mat_inout_{}", m);
    auto* q = inoutQFuncs.Get(name);
    SLIC_ASSERT(q != nullptr);
    double* rData = q->GetData();

    for(int j = 0; j < SZ; ++j)
    {
      mData[j] = rData[j] > 0 ? 0 : mData[j];
    }
  }
}

/**
 * Compute volume fractions function for shape on a grid of resolution \a gridRes
 * in region defined by bounding box \a queryBounds
 */
void computeVolumeFractions(const std::string& shapeName,
                            mfem::DataCollection* dc,
                            QFunctionCollection& inoutQFuncs,
                            int outputOrder)
{
  auto inoutName = fmt::format("mat_inout_{}", shapeName);
  auto volFracName = fmt::format("vol_frac_{}", shapeName);

  // Grab a pointer to the inout samples QFunc
  mfem::QuadratureFunction* inout = inoutQFuncs.Get(inoutName);

  const int sampleOrder = inout->GetSpace()->GetElementIntRule(0).GetOrder();
  const int sampleNQ = inout->GetSpace()->GetElementIntRule(0).GetNPoints();
  const int sampleSZ = inout->GetSpace()->GetSize();
  SLIC_INFO(fmt::format(std::locale("en_US.UTF-8"),
                        "In computeVolumeFractions(): sample order {} | "
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

  // Project QField onto volume fractions field

  mfem::L2_FECollection* fec =
    new mfem::L2_FECollection(outputOrder, dim, mfem::BasisType::Positive);
  mfem::FiniteElementSpace* fes = new mfem::FiniteElementSpace(mesh, fec);
  mfem::GridFunction* volFrac = new mfem::GridFunction(fes);
  volFrac->MakeOwner(fec);
  dc->RegisterField(volFracName, volFrac);

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
  SLIC_INFO(
    fmt::format(std::locale("en_US.UTF-8"),
                "\t Generating volume fractions '{}' took {:.3f} seconds (@ "
                "{:L} dofs processed per second)",
                volFracName,
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
          reader.setFileName(shapeSet.resolvePath(geom.getPath()));
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

  mesh->EnsureNodes();
  dc->SetMeshNodesName("positions");
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
  params.shapeSet = klee::readShapeSet(params.shapeFile);

  params.initializeProblemBoundingBox();

  auto bbox = params.problemBoundingBox();
  SLIC_INFO("Mesh bounding box: " << bbox);

  // Add the mesh to an mfem data collection
  axom::sidre::MFEMSidreDataCollection dc("shaping", nullptr, true);

  initializeMesh(params, &dc);

  // TODO: Use MfemSidreDataCollection QFuncs for this when we upgrade to post mfem@4.3
  QFunctionCollection inoutQFuncs;
  DenseTensorCollection inoutDofs;

  // Sample the InOut quadrature field for each shape using an InOut octree
  // Assumptions: Each shape has a unique name
  SLIC_INFO(fmt::format("{:=^80}", "Sampling InOut fields for shapes"));
  for(const auto& s : params.shapeSet.getShapes())
  {
    SLIC_ASSERT(s.getGeometry().getFormat() == "stl");

    const std::string shapeName = s.getName();

    std::string stlPath = params.shapeSet.resolvePath(s.getGeometry().getPath());
    std::string outMsg = fmt::format(" Loading mesh '{}' ", shapeName);
    SLIC_INFO(fmt::format("{:-^80}", outMsg));
    SLIC_INFO("Reading file: " << stlPath << "...");

    quest::STLReader reader;
    reader.setFileName(stlPath);
    reader.read();

    // Create surface mesh
    mint::Mesh* surface_mesh = new UMesh(3, mint::TRIANGLE);
    reader.getMesh(static_cast<UMesh*>(surface_mesh));

    // Compute mesh bounding box and log some stats about the surface
    GeometricBoundingBox meshBB = compute_bounds(surface_mesh);
    SLIC_INFO("Mesh bounding box: " << meshBB);

    // Create octree over mesh's bounding box
    SLIC_INFO(fmt::format("{:-^80}", " Generating the octree "));
    Octree3D octree(meshBB, surface_mesh);
    octree.generateIndex();

    {
      const int nVerts = surface_mesh->getNumberOfNodes();
      const int nCells = surface_mesh->getNumberOfCells();

      SLIC_INFO(fmt::format(
        "After welding, surface mesh has {} vertices  and {} triangles.",
        nVerts,
        nCells));
      mint::write_vtk(surface_mesh,
                      fmt::format("meldedTriMesh_{}.vtk", shapeName));
    }

    SLIC_INFO(fmt::format("{:-^80}", " Querying the octree "));

    // Sample the InOut field at the mesh quadrature points
    const int sampleOrder = params.quadratureOrder;
    const int outputOrder = params.outputOrder;
    switch(params.vfSampling)
    {
    case SAMPLE_AT_QPTS:
      sampleInOutField(shapeName, octree, &dc, inoutQFuncs, sampleOrder);
      break;
    case SAMPLE_AT_DOFS:
      computeVolumeFractionsBaseline(shapeName,
                                     octree,
                                     &dc,
                                     sampleOrder,
                                     outputOrder);
      break;
    }

    delete surface_mesh;
    surface_mesh = nullptr;
  }

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
      mergeQFuncs(mat, kv.second, inoutQFuncs);
    }

    // Merge all shapes of a given material into a single material QFunc
    for(const auto& kv : replaced_by)
    {
      auto& mat = kv.first;
      replaceMaterials(mat, kv.second, inoutQFuncs);
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
    case SAMPLE_AT_QPTS:
      computeVolumeFractions(shapeName, &dc, inoutQFuncs, outputOrder);
      break;
    case SAMPLE_AT_DOFS:
      /* no-op for now */
      break;
    }
  }

  // Save meshes and fields
  dc.Save();

  MPI_Finalize();

  return 0;
}
