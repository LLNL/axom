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

#include "mfem.hpp"

#include "fmt/fmt.hpp"
#include "fmt/locale.h"

namespace axom
{
namespace quest
{
namespace shaping
{
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

}  // end namespace shaping
}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_QUEST_SHAPING__HPP_
