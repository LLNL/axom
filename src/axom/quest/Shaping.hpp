// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file Shaping.hpp
 *
 * \brief Helper class for shaping queries
 */

#ifndef AXOM_QUEST_SHAPING__HPP_
#define AXOM_QUEST_SHAPING__HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"
#include "axom/spin.hpp"
#include "axom/klee.hpp"

#include "axom/quest/InOutOctree.hpp"
#include "axom/quest/interface/internal/QuestHelpers.hpp"

#include "mfem.hpp"

#include "fmt/fmt.hpp"
#include "fmt/locale.h"

namespace axom
{
namespace quest
{
namespace shaping
{
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

enum class VolFracSampling : int
{
  SAMPLE_AT_DOFS,
  SAMPLE_AT_QPTS
};

/** Computes the bounding box of the surface mesh */
GeometricBoundingBox compute_bounds(const mint::Mesh& mesh)
{
  GeometricBoundingBox meshBB;
  SpacePt pt;

  for(int i = 0; i < mesh.getNumberOfNodes(); ++i)
  {
    mesh.getNode(i, pt.data());
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
  // Note: indexing F(i,j) where  0 <= j < i < s
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

}  // end namespace shaping

namespace unused
{
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
}  // end namespace unused

template <int NDIMS>
class SamplingShaper
{
public:
  static constexpr int DIM = NDIMS;
  using InOutOctreeType = quest::InOutOctree<DIM>;

  using GeometricBoundingBox = typename InOutOctreeType::GeometricBoundingBox;
  using SpacePt = typename InOutOctreeType::SpacePt;
  using SpaceVector = typename InOutOctreeType::SpaceVector;
  using GridPt = typename InOutOctreeType::GridPt;
  using BlockIndex = typename InOutOctreeType::BlockIndex;

public:
  SamplingShaper(const std::string& shapeName, mint::Mesh* surfaceMesh)
    : m_shapeName(shapeName)
    , m_surfaceMesh(surfaceMesh)
  { }

  ~SamplingShaper() { delete m_octree; }

  mint::Mesh* getSurfaceMesh() const { return m_surfaceMesh; }

  /// Computes the bounding box of the surface mesh
  void computeBounds()
  {
    SLIC_ASSERT(m_surfaceMesh != nullptr);

    m_bbox.clear();
    SpacePt pt;

    for(int i = 0; i < m_surfaceMesh->getNumberOfNodes(); ++i)
    {
      m_surfaceMesh->getNode(i, pt.data());
      m_bbox.addPoint(pt);
    }

    SLIC_ASSERT(m_bbox.isValid());

    SLIC_INFO("Mesh bounding box: " << m_bbox);
  }

  void initSpatialIndex()
  {
    // Create octree over mesh's bounding box
    m_octree = new InOutOctreeType(m_bbox, m_surfaceMesh);
    m_octree->generateIndex();
  }

  void sampleInOutField(mfem::DataCollection* dc,
                        shaping::QFunctionCollection& inoutQFuncs,
                        int sampleRes)
  {
    auto* mesh = dc->GetMesh();
    SLIC_ASSERT(mesh != nullptr);
    const int NE = mesh->GetNE();
    const int dim = mesh->Dimension();

    // Generate a Quadrature Function with the geometric positions, if not already available
    if(!inoutQFuncs.Has("positions"))
    {
      shaping::generatePositionsQFunction(mesh, inoutQFuncs, sampleRes);
    }

    // Access the positions QFunc and associated QuadratureSpace
    mfem::QuadratureFunction* pos_coef = inoutQFuncs.Get("positions");
    mfem::QuadratureSpace* sp = pos_coef->GetSpace();
    const int nq = sp->GetElementIntRule(0).GetNPoints();

    // Sample the in/out field at each point
    // store in QField which we register with the QFunc collection
    const std::string inoutName = fmt::format("inout_{}", m_shapeName);
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
        const SpacePt pt(m.GetColumn(p), dim);
        const bool in = m_octree->within(pt);
        res(p) = in ? 1. : 0.;

        // SLIC_INFO(fmt::format("[{},{}] Pt: {}, In: {}", i,p,pt, (in? "yes" : "no") ));
      }
    }
    timer.stop();

    SLIC_INFO(fmt::format(std::locale("en_US.UTF-8"),
                          "\t Sampling inout field '{}' took {} seconds (@ "
                          "{:L} queries per second)",
                          inoutName,
                          timer.elapsed(),
                          static_cast<int>((NE * nq) / timer.elapsed())));
  }

  /**
  * Compute volume fractions function for shape on a grid of resolution \a gridRes
  * in region defined by bounding box \a queryBounds
  */
  void computeVolumeFractionsBaseline(mfem::DataCollection* dc,
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
    auto volFracName = fmt::format("vol_frac_{}", m_shapeName);
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
        const SpacePt pt(m.GetColumn(p), dim);
        const bool in = m_octree->within(pt);
        res(p) = in ? 1. : 0.;
      }

      fes->GetElementDofs(i, dofs);
      volFrac->SetSubVector(dofs, res);
    }
  }

private:
  std::string m_shapeName;
  GeometricBoundingBox m_bbox;
  mint::Mesh* m_surfaceMesh {nullptr};
  InOutOctreeType* m_octree {nullptr};
};

/// Helper class to use an MFEM sampling-based shaper
class MFEMShaping
{
public:
  MFEMShaping(const klee::ShapeSet& shapeSet, sidre::MFEMSidreDataCollection* dc)
    : m_shapeSet(shapeSet)
    , m_dc(dc)
  {
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
    m_comm = m_dc->GetComm();
#endif
  }

  void setSamplesPerKnotSpan(int nSamples)
  {
    using axom::utilities::clampLower;
    SLIC_WARNING_IF(
      nSamples < 1,
      fmt::format(
        "Samples per knot span must be at least 1. Provided value was {}",
        nSamples));

    m_samplesPerKnotSpan = clampLower(nSamples, 1);
  }

  void setVertexWeldThreshold(double threshold)
  {
    SLIC_WARNING_IF(
      threshold <= 0.,
      fmt::format(
        "Vertex weld threshold should be positive Provided value was {}",
        threshold));

    m_vertexWeldThreshold = threshold;
  }

  sidre::MFEMSidreDataCollection* getDC() { return m_dc; }
  mint::Mesh* getSurfaceMesh() const { return m_surfaceMesh; }

  /// Loads the shape from file into m_surfaceMesh
  void loadShape(const klee::Shape& shape)
  {
    using axom::utilities::string::endsWith;
    SLIC_ASSERT(shape.getGeometry().getFormat() == "stl" ||
                shape.getGeometry().getFormat() == "c2c");

    const std::string shapeName = shape.getName();
    std::string outMsg = fmt::format(" Loading shape '{}' ", shapeName);
    SLIC_INFO(fmt::format("{:-^80}", outMsg));

    std::string shapePath = m_shapeSet.resolvePath(shape.getGeometry().getPath());
    SLIC_INFO("Reading file: " << shapePath << "...");

    if(endsWith(shapePath, ".stl"))
    {
      quest::internal::read_stl_mesh(shapePath, m_surfaceMesh, m_comm);
    }
#ifdef AXOM_USE_C2C
    else if(endsWith(shapePath, ".contour"))
    {
      quest::internal::read_c2c_mesh(shapePath,
                                     m_samplesPerKnotSpan,
                                     m_vertexWeldThreshold,
                                     m_surfaceMesh,
                                     m_comm);
    }
#endif
    else
    {
      SLIC_ERROR(
        fmt::format("Unsupported filetype for this Axom configuration. "
                    "Provided file was '{}'",
                    shapePath));
    }
  }

  void applyTransforms(const klee::Shape& shape)
  {
    // TODO: Implement this as a set of affine transforms to vertices of mesh
    AXOM_UNUSED_VAR(shape);
  }

  void prepareShapeQuery(klee::Dimensions shapeDimension, const klee::Shape& shape)
  {
    const auto& shapeName = shape.getName();

    SLIC_INFO(fmt::format("{:-^80}", " Generating the octree "));
    switch(shapeDimension)
    {
    case klee::Dimensions::Two:
      m_samplingShaper2D = new SamplingShaper<2>(shapeName, m_surfaceMesh);
      m_samplingShaper2D->computeBounds();
      m_samplingShaper2D->initSpatialIndex();
      m_surfaceMesh = m_samplingShaper2D->getSurfaceMesh();
      break;

    case klee::Dimensions::Three:
      m_samplingShaper3D = new SamplingShaper<3>(shapeName, m_surfaceMesh);
      m_samplingShaper3D->computeBounds();
      m_samplingShaper3D->initSpatialIndex();
      m_surfaceMesh = m_samplingShaper3D->getSurfaceMesh();
      break;

    default:
      SLIC_ERROR(
        "Shaping dimension must be 2 or 3, but requested dimension was "
        << static_cast<int>(shapeDimension));
      break;
    }

    // Check that one of sampling shapers (2D or 3D) is null and the other is not
    SLIC_ASSERT((m_samplingShaper2D == nullptr && m_samplingShaper3D != nullptr) ||
                (m_samplingShaper3D == nullptr && m_samplingShaper2D != nullptr));

    // Output some logging info and dump the mesh
    {
      const int nVerts = m_surfaceMesh->getNumberOfNodes();
      const int nCells = m_surfaceMesh->getNumberOfCells();

      SLIC_INFO(fmt::format(
        "After welding, surface mesh has {} vertices  and {} triangles.",
        nVerts,
        nCells));
      mint::write_vtk(m_surfaceMesh,
                      fmt::format("meldedTriMesh_{}.vtk", shapeName));
    }
  }

  // Handles 2D or 3D shaping, based on the template and associated parameter
  template <typename SamplingShaperType>
  void runShapeQueryImpl(SamplingShaperType* shaper,
                         shaping::VolFracSampling vfSampling,
                         int samplingOrder,
                         int outputOrder)
  {
    // Sample the InOut field at the mesh quadrature points
    switch(vfSampling)
    {
    case shaping::VolFracSampling::SAMPLE_AT_QPTS:
      shaper->sampleInOutField(m_dc, m_inoutQFuncs, samplingOrder);
      break;
    case shaping::VolFracSampling::SAMPLE_AT_DOFS:
      shaper->computeVolumeFractionsBaseline(m_dc, samplingOrder, outputOrder);
      break;
    }
  }

  void runShapeQuery(shaping::VolFracSampling vfSampling,
                     int samplingOrder,
                     int outputOrder)
  {
    SLIC_INFO(fmt::format("{:-^80}", " Querying the octree "));

    switch(getShapeDimension())
    {
    case klee::Dimensions::Two:
      runShapeQueryImpl(m_samplingShaper2D, vfSampling, samplingOrder, outputOrder);
      break;
    case klee::Dimensions::Three:
      runShapeQueryImpl(m_samplingShaper3D, vfSampling, samplingOrder, outputOrder);
      break;
    }
  }

  void finalizeShapeQuery()
  {
    delete m_samplingShaper2D;
    m_samplingShaper2D = nullptr;

    delete m_samplingShaper3D;
    m_samplingShaper3D = nullptr;

    delete m_surfaceMesh;
    m_surfaceMesh = nullptr;
  }

  shaping::QFunctionCollection& getInoutQFuncs() { return m_inoutQFuncs; }
  // shaping::DenseTensorCollection getDenseTensorCollection() const { return m_inoutDofs; }

private:
  klee::Dimensions getShapeDimension() const
  {
    const bool has2D = (m_samplingShaper2D != nullptr);
    const bool has3D = (m_samplingShaper3D != nullptr);
    SLIC_ERROR_IF(!(has2D || has3D), "Shape not initialized");
    SLIC_ERROR_IF(has2D && has3D, "Cannot have concurrent 2D and 3D shapes");

    return has2D ? klee::Dimensions::Two : klee::Dimensions::Three;
  }

public:
  const klee::ShapeSet& m_shapeSet;
  sidre::MFEMSidreDataCollection* m_dc;

  // TODO: Use MfemSidreDataCollection QFuncs for this when we upgrade to post mfem@4.3
  shaping::QFunctionCollection m_inoutQFuncs;
  shaping::DenseTensorCollection m_inoutDofs;

  mint::Mesh* m_surfaceMesh {nullptr};
  SamplingShaper<2>* m_samplingShaper2D {nullptr};
  SamplingShaper<3>* m_samplingShaper3D {nullptr};

  int m_samplesPerKnotSpan {25};
  double m_vertexWeldThreshold {1e-9};

  MPI_Comm m_comm {MPI_COMM_SELF};
};

}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_QUEST_SHAPING__HPP_
