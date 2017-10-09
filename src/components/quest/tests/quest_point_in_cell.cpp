/*
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */



#include "gtest/gtest.h"

#include "axom_utils/Timer.hpp"

#include "mint/CurvilinearMesh.hpp"

#include "primal/Point.hpp"
#include "primal/BoundingBox.hpp"

#include "quest/ImplicitGrid.hpp"
#include "quest/PointInCell.hpp"

#include "quest_test_utilities.hpp"

#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

#include "fmt/format.h"

#include <sstream>
#include <vector>
#include <cmath>    // for pow
#include <cstdlib>  // for srand

namespace
{
  const unsigned int SRAND_SEED = 42; //105;

  const bool outputMeshMFEM = true;
  const bool outputMeshVTK = false;
}

namespace types2D
{
  const int DIM = 2;
  const int ELT_MULT_FAC = 4;

  const double EPS = 1e-8;

#ifdef AXOM_DEBUG
  const int NREFINE = 5;
  const int NUM_TEST_PTS = 10000;
#else
  const int NREFINE = 5;
  const int NUM_TEST_PTS = 100000;
#endif

  typedef axom::primal::BoundingBox<double, DIM> BBox;
  typedef axom::primal::Point<double,DIM> SpacePt;
  typedef axom::primal::Vector<double,DIM> SpaceVec;
  typedef axom::primal::Point<int,DIM> GridCell;
  typedef axom::quest::PointInCell<DIM> PointInCellType;
}

std::string meshPrefix =
    "MFEM mesh v1.0"  "\n\n"
    "dimension"         "\n"
    "2"               "\n\n"
    "elements"          "\n"
    "1"                 "\n"
    "1 3 0 1 2 3"     "\n\n"
    "boundary"          "\n"
    "0"               "\n\n";

std::string lowOrderVerts =
    "vertices"       "\n"
    "4"              "\n"
    "2"              "\n"
    "0     -{0}"     "\n"
    "{0}      0"     "\n"
    "0      {0}"     "\n"
    "-{0}     0"     "\n";

std::string highOrderNodes =
    "vertices"                      "\n"
    "4"                           "\n\n"
    "nodes"                         "\n"
    "FiniteElementSpace"            "\n"
    "FiniteElementCollection: {2}"  "\n"
    "VDim: 2"                       "\n"
    "Ordering: 1"                   "\n"
    "   0  -{0}"                    "\n"
    " {0}    0"                     "\n"
    "   0  {0}"                     "\n"
    "-{0}    0"                     "\n"
    " {1} -{1}"                     "\n"
    " {1}  {1}"                     "\n"
    "-{1}  {1}"                     "\n"
    "-{1} -{1}"                     "\n"
    "  0     0"                     "\n";


std::string c_shapedNodes =
    "vertices"                      "\n"
    "4"                           "\n\n"
    "nodes"                         "\n"
    "FiniteElementSpace"            "\n"
    "FiniteElementCollection: {0}"  "\n"
    "VDim: 2"                       "\n"
    "Ordering: 1"                   "\n"
    "0 0"                           "\n"
    "0 2"                           "\n"
    "0 6"                           "\n"
    "0 8"                           "\n"
    "0 1"                           "\n"
    "-6 4"                           "\n"
    "0 7"                           "\n"
    "-8 4"                           "\n"
    "-7 4"                           "\n";

std::vector< axom::primal::Point<double, 2> > generate2DTestPoints(double val)
{
  typedef axom::primal::Point<double, 2> SpacePt2D;
  const int SZ = types2D::NUM_TEST_PTS;

  std::vector<SpacePt2D> pts;
  pts.reserve( SZ + 10 );

  pts.push_back( SpacePt2D::zero() );
  pts.push_back( SpacePt2D::ones() );
  pts.push_back( SpacePt2D::make_point(.1, .1) );
  pts.push_back( SpacePt2D::make_point(.4, .4) );
  pts.push_back( SpacePt2D::make_point(val, 0) );
  pts.push_back( SpacePt2D::make_point(val  +0.005, 0) );


  for(int i=0; i< SZ; ++i)
  {
    pts.push_back( axom::quest::utilities::randomSpacePt<2>(0, 1.25 * val));
  }

  return pts;
}

std::vector< axom::primal::Point<double, 2> > generate2DIsoParTestPoints(int res)
{
  typedef axom::primal::Point<double, 2> SpacePt2D;

  std::vector<SpacePt2D> pts;

  for(int i=1; i < res; ++i)
    for(int j=1; j < res; ++j)
    {
      SpacePt2D pt = SpacePt2D::make_point( static_cast<double>(i)/res,
                                            static_cast<double>(j)/res);

      SpacePt2D off = axom::quest::utilities::randomSpacePt<2>( -types2D::EPS, +types2D::EPS);

      pts.push_back( pt.array() + off.array() );
    }

  return pts;
}


struct ExpectedValueFlat
{
  static const int DIM = 2;
  typedef typename axom::primal::Point<double, DIM> SpacePt;
  typedef typename axom::primal::Vector<double, DIM> SpaceVec;

  ExpectedValueFlat(double radius) : m_radius(radius) {}

  double radius() const { return m_radius; }

  bool canTestPoint(const SpacePt& pt)
  {
    return true;
  }

  bool expectedInMesh(const SpacePt& pt)
  {
    // We expect the point to be in the mesh if its L1 norm is less than or equal to m_radus
    double normL1 = axom::utilities::abs(pt[0]) +  axom::utilities::abs(pt[1]);
    return normL1 < m_radius || axom::utilities::isNearlyEqual(normL1, m_radius);
  }

  double m_radius;
};

// We have an inner radius DIAG_VAL where the point should be found in the mesh
// and an outer radius VERT_VAL, where the point should not be found
// there is an annulus between the two radii,
// where it is not easy to be sure analytically if the point should be in or out
struct ExpectedValueCurved
{
  static const int DIM = 2;
  typedef typename axom::primal::Point<double, DIM> SpacePt;
  typedef typename axom::primal::Vector<double, DIM> SpaceVec;

  ExpectedValueCurved(double radius) : m_radius(radius) {}

  double radius() const { return m_radius; }

  bool canTestPoint(const SpacePt& pt) const
  {
    double normL2 = SpaceVec( pt.array() ).norm();
    return !axom::utilities::isNearlyEqual(normL2, m_radius, .05 * m_radius);
  }

  bool expectedInMesh(const SpacePt& pt) const
  {
    double norm = SpaceVec( pt.array() ).norm();
    return (norm < m_radius) ? true : false;
  }

  double m_radius;
};


/*!
 * Compute the number of expected elements in uniform refinement with
 * a given refinement factor and level of resolution
 */
int expectedNumElts(int refinementFactor, int refinementLevel)
{
  double refFac = static_cast<double>(refinementFactor);
  double refLev = static_cast<double>(refinementLevel);

  return static_cast<int>(std::pow(refFac,refLev) );
}

/**
 * Jitter all degrees of freedom (dofs) of the mesh's nodal grid function
 * Implementation borrowed from mfem's mesh-explorer mini-app
 */
template<int DIM>
void jitterNodalValues(mfem::Mesh* mesh, double dx)
{
  mfem::GridFunction* nodes = mesh->GetNodes();
  if(nodes == AXOM_NULLPTR)
    return;

  mfem::FiniteElementSpace *fespace = nodes->FESpace();
  mfem::GridFunction rdm(fespace);
  rdm.Randomize( SRAND_SEED );
  rdm -= 0.5; // shift to random values in [-0.5,0.5]
  rdm *= dx;

  // compute minimal local mesh size
  mfem::Vector h0(fespace->GetNDofs());
  h0 = std::numeric_limits<double>::infinity();
  {
     mfem::Array<int> dofs;
     for (int i = 0; i < fespace->GetNE(); i++)
     {
        fespace->GetElementDofs(i, dofs);
        for (int j = 0; j < dofs.Size(); j++)
        {
           h0(dofs[j]) = std::min(h0(dofs[j]), mesh->GetElementSize(i));
        }
     }
  }

  // scale the random values to be of order of the local mesh size
  for (int i = 0; i < fespace->GetNDofs(); i++)
  {
     for (int d = 0; d < DIM; d++)
     {
        rdm(fespace->DofToVDof(i,d)) *= h0(i)/2.;
     }
  }

  // don't perturb the boundary
  bool keepBdry = true;
  if (keepBdry)
  {
     mfem::Array<int> vdofs;
     for (int i = 0; i < fespace->GetNBE(); i++)
     {
        fespace->GetBdrElementVDofs(i, vdofs);
        for (int j = 0; j < vdofs.Size(); j++)
        {
           rdm(vdofs[j]) = 0.0;
        }
     }
  }

  // Finally, add the displacements
  *nodes +=rdm;
}

template<typename ExpectedValueFunctor>
void testRandomPointsOnMesh(mfem::Mesh* mesh, ExpectedValueFunctor exp, const std::string& meshTypeStr)
{
  using namespace types2D;

  // output mesh in mfem and vtk format
  bool isCurved = mesh->GetNodes() != AXOM_NULLPTR;
  std::string filename = fmt::format("quest_point_in_cell_{}_quad",meshTypeStr);

  if(outputMeshMFEM)
  {
    mfem::VisItDataCollection dataCol(filename, mesh);
    if(mesh->GetNodes())
      dataCol.RegisterField("nodes", mesh->GetNodes());
    dataCol.Save();
  }
  if(outputMeshVTK)
  {
    std::ofstream mfem_stream(filename + ".vtk");
    mesh->PrintVTK(mfem_stream);
  }

  // Add mesh to the grid
  axom::utilities::Timer constructTimer(true);
  PointInCellType spatialIndex(mesh, GridCell(25));
  SLIC_INFO(fmt::format("Constructing index over {} quad mesh with {} elems took {} s",
      meshTypeStr,
      mesh->GetNE(),
      constructTimer.elapsed() ) );

  std::vector<SpacePt> pts = generate2DTestPoints( exp.radius() );
  SpacePt isoPar;
  int numCheckedPoints = 0;
  int numInverseXforms = 0;

  axom::utilities::Timer queryTimer(true);
  for( std::vector<SpacePt>::iterator it= pts.begin(); it != pts.end(); ++it)
  {
    SpacePt& queryPoint = *it;

    // Try to find the point
    int idx = spatialIndex.locatePoint( queryPoint, &isoPar );
    bool isInMesh = (idx != PointInCellType::NO_CELL);

    // Check if result matches our expectations (our simple model
    // supports verifying many, but not all query points)
    if( exp.canTestPoint(queryPoint) )
    {
      ++numCheckedPoints;

      bool expectedInMesh = exp.expectedInMesh(queryPoint);
      EXPECT_EQ( expectedInMesh, isInMesh) << "Point " << *it;
    }

    // Check if the transform's inverse gives us back our original point
    if(isInMesh)
    {
      ++numInverseXforms;

      SpacePt untransformPt = spatialIndex.findInSpace(idx, isoPar);

      for(int i=0; i< DIM; ++i)
      {
        EXPECT_NEAR(queryPoint[i],untransformPt[i], EPS);
      }
    }
  }


  SLIC_INFO(fmt::format("Querying {} pts on {} quad mesh took {} s -- rate: {} q/s (includes {} inverse xforms)",
      pts.size(),
      meshTypeStr,
      queryTimer.elapsed(),
      pts.size() / queryTimer.elapsed(),
      numInverseXforms) );

  if(isCurved)
  {
    SLIC_INFO(fmt::format("On {} mesh, verified plausibility in {} of {} cases ({:.1f}%)",
        meshTypeStr, numCheckedPoints, pts.size(), 100* static_cast<double>(numCheckedPoints)/ pts.size() ));
  }
}

void testIsoGridPointsOnMesh(mfem::Mesh* mesh, const std::string& meshTypeStr)
{
  using namespace types2D;

  std::string filename = fmt::format("quest_point_in_cell_{}_quad",meshTypeStr);

  // Add mesh to the grid
  axom::utilities::Timer constructTimer(true);
  PointInCellType spatialIndex(mesh, GridCell(25));
  SLIC_INFO(fmt::format("Constructing index over {} quad mesh with {} elems took {} s",
      meshTypeStr,
      mesh->GetNE(),
      constructTimer.elapsed() ) );


  // Test that a fixed set of isoparametric coords on each cell maps to the correct place.
  std::vector<SpacePt> pts = generate2DIsoParTestPoints(10);
  axom::utilities::Timer queryTimer2(true);
  SpacePt foundIsoPar;
  for(int eltId=0; eltId< mesh->GetNE(); ++eltId)
  {
    for( std::vector<SpacePt>::iterator it= pts.begin(); it != pts.end(); ++it)
    {
      SpacePt& isoparCenter = *it;
      SpacePt spacePt = spatialIndex.findInSpace(eltId, isoparCenter);

      int foundCellId = spatialIndex.locatePoint(spacePt, &foundIsoPar);

      // Check that the cell ids agree
      EXPECT_EQ( eltId, foundCellId)
          << fmt::format("For element {} -- computed space point {} from isoPar {} -- found isoPar is {}",
              eltId, spacePt, isoparCenter, foundIsoPar);

      // Check that the isoparametric coordinates agree
      for(int i=0; i< DIM; ++i)
      {
        EXPECT_NEAR(isoparCenter[i],foundIsoPar[i], EPS);
      }

    }
  }

  SLIC_INFO(fmt::format("Verifying {} pts on {} quad mesh took {} s -- rate: {} q/s",
      pts.size() * mesh->GetNE(),
      meshTypeStr,
      queryTimer2.elapsed(),
      pts.size() * mesh->GetNE() / queryTimer2.elapsed()  ) );

}

TEST( quest_point_in_cell, simple_bilinear_mesh )
{
  using namespace types2D;
  const double VAL = 0.5;

  SLIC_INFO("Testing that we can create a simple mfem mesh. "
      << "No point-in-cell queries");


  // Create a simple mfem mesh over a single quad
  //   -- a diamond with verts at +-(VAL,0) and +-(0, VAL)
  // This is a circle with radius VAL in the L1 norm
  std::stringstream sstr;
  sstr << meshPrefix
       << fmt::format(lowOrderVerts, VAL);
  mfem::Mesh mesh(sstr);

  // Add a bilinear gridfunction
  mfem::FiniteElementCollection* fec = mfem::FiniteElementCollection::New("Linear");
  mfem::FiniteElementSpace *fes = new mfem::FiniteElementSpace(&mesh, fec, 1);
  mfem::GridFunction gf(fes);
  gf(0) = gf(2) = 1.;
  gf(1) = gf(3) = 0.;

  std::string filename = fmt::format("simple_mesh");
  {
    mfem::VisItDataCollection dataCol(filename, &mesh);
    dataCol.RegisterField("data", &gf);
    dataCol.Save();
  }

}

TEST( quest_point_in_cell, pic_flat_quad )
{
  using namespace types2D;
  const double VAL = 0.5;

  SLIC_INFO("Construction and querying PointInCell structure"
      << " over flat quad mesh with " << NREFINE << " refinement levels.");


  // Create a simple mfem mesh over a single quad
  //   -- a diamond with verts at +-(VAL,0) and +-(0, VAL)
  // This is a circle with radius VAL in the L1 norm
  std::stringstream sstr;
  sstr << meshPrefix
       << fmt::format(lowOrderVerts, VAL);

  // Refine the mesh several times
  mfem::Mesh mesh(sstr);
  for(int i=0; i< NREFINE; ++i)
    mesh.UniformRefinement();

  // Sanity checks
  int expectedNE = expectedNumElts(ELT_MULT_FAC, NREFINE);
  EXPECT_EQ( expectedNE, mesh.GetNE() );
  EXPECT_EQ(AXOM_NULLPTR, mesh.GetNodes());

  SCOPED_TRACE("point_in_cell_flat_refined");
  std::string meshTypeStr = fmt::format("flat_{}", NREFINE);;
  testRandomPointsOnMesh(&mesh, ExpectedValueFlat(VAL), meshTypeStr);
  testIsoGridPointsOnMesh(&mesh, meshTypeStr);
}

TEST( quest_point_in_cell, pic_curved_quad )
{
  using namespace types2D;

  const double VERT_VAL = 0.5;
  const double DIAG_VAL = VERT_VAL * std::sqrt(2.)/2.;
  const std::string feCollName = "Quadratic";

  SLIC_INFO("Construction and querying PointInCell structure"
      << " over quadratic quad mesh with 0 and " << NREFINE << " refinement levels.");


  // Create a simple high order mfem mesh over a single quad
  //  -- a diamond with verts at +-(VAL,0) and +-(0, VAL)
  //  -- and nodes along diagonals (+- VAL2, +-VAL2}
  std::stringstream sstr;
  sstr << meshPrefix
       << fmt::format(highOrderNodes, VERT_VAL, DIAG_VAL, feCollName);

  mfem::Mesh mesh(sstr);

  // Test single element mesh
  {
    SCOPED_TRACE("point_in_cell_curved_level_single_elt");
    std::string meshTypeStr = fmt::format("{}_{}", feCollName, 1);
    testRandomPointsOnMesh(&mesh, ExpectedValueCurved(VERT_VAL), meshTypeStr);
    testIsoGridPointsOnMesh(&mesh, meshTypeStr);
  }

  for(int i=0; i< NREFINE; ++i)
    mesh.UniformRefinement();

  int expectedNE = expectedNumElts(ELT_MULT_FAC, NREFINE);
  EXPECT_EQ( expectedNE, mesh.GetNE() );
  EXPECT_NE(AXOM_NULLPTR, mesh.GetNodes());
  EXPECT_EQ(feCollName, mesh.GetNodalFESpace()->FEColl()->Name());

  // Test refined mesh
  {
    SCOPED_TRACE("point_in_cell_curved_refined");
    std::string meshTypeStr = fmt::format("{}_{}", feCollName, NREFINE);
    testRandomPointsOnMesh(&mesh, ExpectedValueCurved(VERT_VAL), meshTypeStr);
    testIsoGridPointsOnMesh(&mesh, meshTypeStr);
  }
}

TEST( quest_point_in_cell, pic_curved_quad_jittered )
{
  using namespace types2D;

  const double VERT_VAL = 0.5;
  const double DIAG_VAL = VERT_VAL * std::sqrt(2.)/2.;
  const std::string feCollName = "Quadratic";
  const double jitterFactor = .15;

  SLIC_INFO("Construction and querying PointInCell structure"
      << " over jittered quadratic quad mesh with " << NREFINE << " refinement levels.");


  // Create a simple high order mfem mesh over a single quad
  //  -- a diamond with verts at +-(VAL,0) and +-(0, VAL) and nodes along diagonals (+- DIAG_VAL}
  std::stringstream sstr;
  sstr << meshPrefix
       << fmt::format(highOrderNodes, VERT_VAL, DIAG_VAL, feCollName);

  mfem::Mesh mesh(sstr);

  // Test single element mesh
  {
    mesh.UniformRefinement();
    jitterNodalValues<DIM>(&mesh, jitterFactor);

    SCOPED_TRACE("point_in_cell_curved_jittered_single_elt");
    std::string meshTypeStr = fmt::format("{}_jittered_{}", feCollName, 1);
    testRandomPointsOnMesh(&mesh, ExpectedValueCurved(VERT_VAL), meshTypeStr);
    testIsoGridPointsOnMesh(&mesh, meshTypeStr);
  }

  // Test refined mesh
  {
    for(int i=0; i< NREFINE -1; ++i)
    {
      mesh.UniformRefinement();
      jitterNodalValues<DIM>(&mesh, jitterFactor);
    }

    SCOPED_TRACE("point_in_cell_curved_jittered_refined");
    std::string meshTypeStr = fmt::format("{}_jittered_{}", feCollName, NREFINE);
    testRandomPointsOnMesh(&mesh, ExpectedValueCurved(VERT_VAL), meshTypeStr);
    testIsoGridPointsOnMesh(&mesh, meshTypeStr);
  }
}


TEST( quest_point_in_cell, pic_curved_quad_jittered_pos )
{
  using namespace types2D;

  const double VERT_VAL = 0.5;
  const double DIAG_VAL = VERT_VAL * (std::sqrt(2.) - 0.5);
  const std::string feCollName = "QuadraticPos";
  const double jitterFactor = .15;

  SLIC_INFO("Construction and querying PointInCell structure"
      << " over jittered quadratic quad mesh with " << NREFINE << " refinement levels."
      << " Mesh encoded in positive (Bernstein) basis");

  // Create a simple high order mfem mesh over a single quad
  //  -- a diamond with verts at +-(VAL,0) and +-(0, VAL) and nodes along diagonals (+- DIAG_VAL}
  std::stringstream sstr;
  sstr << meshPrefix
       << fmt::format(highOrderNodes, VERT_VAL, DIAG_VAL, feCollName);

  mfem::Mesh mesh(sstr);

  // Test single element mesh
  {
    mesh.UniformRefinement();
    jitterNodalValues<DIM>(&mesh, jitterFactor);

    SCOPED_TRACE("point_in_cell_curved_jittered_positive_single_elet");
    std::string meshTypeStr = fmt::format("{}_jittered_{}", feCollName, 1);
    testRandomPointsOnMesh(&mesh, ExpectedValueCurved(VERT_VAL), meshTypeStr);
    testIsoGridPointsOnMesh(&mesh, meshTypeStr);
  }

  // Test refined mesh
  {
    for(int i=0; i< NREFINE -1; ++i)
    {
      mesh.UniformRefinement();
      jitterNodalValues<DIM>(&mesh, jitterFactor);
    }

    SCOPED_TRACE("point_in_cell_curved_jittered_positive_refined");
    std::string meshTypeStr = fmt::format("{}_jittered_{}", feCollName, NREFINE);
    testRandomPointsOnMesh(&mesh, ExpectedValueCurved(VERT_VAL), meshTypeStr);
    testIsoGridPointsOnMesh(&mesh, meshTypeStr);
  }
}


TEST( quest_point_in_cell, pic_curved_quad_c_shaped )
{
  using namespace types2D;

  const std::string feCollName = "Quadratic";

  SLIC_INFO("Construction and querying PointInCell structure"
      << " over C-shaped biquadratic element");


  // Here we are testing a very curved C-shaped mesh
  // We are comparing a single element mesh to
  // a refined mesh (one level of refinement), and checking
  // that the PointInCell query gives the same result for both
  // (i.e. both inside or both outside)

  // Note: Mesh1 only has one element, with index 0
  std::stringstream sstr1;
  sstr1 << meshPrefix
       << fmt::format(c_shapedNodes, feCollName);
  mfem::Mesh mesh1(sstr1);

  std::stringstream sstr2;
  sstr2 << meshPrefix
       << fmt::format(c_shapedNodes, feCollName);
  mfem::Mesh mesh2(sstr2);
  for(int i=0; i< 1; ++i)
  {
    mesh2.UniformRefinement();
  }

  // output refined mesh to logger in mfem and vtk formats
  std::string filename = "quest_point_in_cell_c_shaped_quad";
  {
    mfem::VisItDataCollection dataCol(filename + "001_mesh", &mesh1);
    dataCol.Save();
  }
  {
    mfem::VisItDataCollection dataCol(filename + "002_mesh", &mesh2);
    dataCol.Save();
  }

  SLIC_INFO("Characteristics for the single element quadratic quad mesh:");
  mesh1.PrintCharacteristics();
  SLIC_INFO("\n-- Elt size sqrt det: " << mesh1.GetElementSize(0, 0) );
  SLIC_INFO("\n-- Elt size h_min: " << mesh1.GetElementSize(0, 1) );
  SLIC_INFO("\n-- Elt size h_max: " << mesh1.GetElementSize(0, 2) );


  // Add mesh to the grid
  axom::utilities::Timer constructTimer(true);
  PointInCellType spatialIndex1(&mesh1, GridCell(25));
  SLIC_INFO(fmt::format("Constructing index over curved quad mesh1 with {} elems took {} s",
      mesh1.GetNE(),
      constructTimer.elapsed() ) );

  axom::utilities::Timer constructTimer2(true);
  PointInCellType spatialIndex2(&mesh2, GridCell(25));
  SLIC_INFO(fmt::format("Constructing index over curved quad mesh2 with {} elems took {} s",
      mesh2.GetNE(),
      constructTimer2.elapsed() ) );

  // For each random point, check that both queries agree
  SLIC_INFO("Querying random points in domain of C-shaped mesh.");
  int numUntransformed = 0;
  axom::utilities::Timer queryTimer(true);
  const int num_pts = NUM_TEST_PTS;
  for( int i=0; i< num_pts; ++i)
  {
    // Create a random point in the domain
    SpacePt queryPoint = axom::quest::utilities::randomSpacePt<DIM>(0, 8.5);
    queryPoint[0] *= -1;

    // Try to find point in mesh1
    SpacePt isoPar1;
    int idx1 = spatialIndex1.locatePoint( queryPoint, &isoPar1);
    bool isInMesh1 = (idx1 != PointInCellType::NO_CELL);

    // Try to find point in mesh2
    SpacePt isoPar2;
    int idx2 = spatialIndex2.locatePoint( queryPoint, &isoPar2);
    bool isInMesh2 = (idx2 != PointInCellType::NO_CELL);

    // Results should be the same
    EXPECT_EQ(isInMesh1, isInMesh2)
        << "Point " << queryPoint
        << " was " << (isInMesh1 ? "in" : "not in") << " mesh 1"
        << " and was " << (isInMesh2 ? "in" : "not in") << " mesh 2"
        << " They should agree.";

    // Try to transform back to space using spatialIndex1
    // Transformed point should match query point
    if(isInMesh1)
    {
      ++numUntransformed;
      SpacePt untransformPt = spatialIndex1.findInSpace(idx1, isoPar1);

      if(! isInMesh2)
      {
        SLIC_INFO("UntransformedPt for mesh 1 was " << untransformPt);
      }

      for(int i=0; i< DIM; ++i)
      {
        EXPECT_NEAR(queryPoint[i],untransformPt[i], EPS);
      }

    }

    // Try to transform back to space using spatialIndex2
    // Transformed point should match query point
    if(isInMesh2)
    {
      ++numUntransformed;
      SpacePt untransformPt = spatialIndex2.findInSpace(idx2, isoPar2);

      if(! isInMesh1)
      {
        SLIC_INFO("UntransformedPt for mesh 2 was " << untransformPt);
      }

      for(int i=0; i< DIM; ++i)
      {
        EXPECT_NEAR(queryPoint[i],untransformPt[i], EPS);
      }

    }
  }

  SLIC_INFO(fmt::format("Querying {} random pts on two C-shaped quadratic quad meshes took {} s -- rate: {} q/s",
      num_pts * 2,
      queryTimer.elapsed(),
      num_pts * 2 / queryTimer.elapsed()  )
    << "\n\t (includes " << numUntransformed << " transfomations back into space)"
  );

    // Test that a fixed set of isoparametric coords on each cell maps to the correct place.
    std::vector<SpacePt> pts = generate2DIsoParTestPoints(10);
    axom::utilities::Timer queryTimer2(true);
    SpacePt foundIsoPar1, foundIsoPar2;
    for(int eltId=0; eltId< mesh1.GetNE(); ++eltId)
    {
      for( std::vector<SpacePt>::iterator it= pts.begin(); it != pts.end(); ++it)
      {
        SpacePt& isoparCenter = *it;
        
        SpacePt spacePt = spatialIndex1.findInSpace(eltId, isoparCenter);
        
        int foundCellId1 = spatialIndex1.locatePoint(spacePt, &foundIsoPar1);
        bool isInMesh1 = (foundCellId1 != PointInCellType::NO_CELL);
        EXPECT_TRUE(isInMesh1) 
            << "Failed for mesh1 on mesh2's cell " << eltId << " with isopar " << isoparCenter << " and spacePt " << spacePt;

        int foundCellId2 = spatialIndex2.locatePoint(spacePt, &foundIsoPar2);
        bool isInMesh2 = (foundCellId2 != PointInCellType::NO_CELL);
        EXPECT_TRUE(isInMesh2) 
            << "Failed for mesh2 on mesh2's cell " << eltId << " with isopar " << isoparCenter << " and spacePt " << spacePt;
        

        // Check that the cell ids agree
        EXPECT_EQ( eltId, foundCellId1)
            << fmt::format("For element {} -- computed midpoint {} @ {} -- found isoPar is {}",
                eltId, spacePt, isoparCenter, foundIsoPar1);

        // Check that the isoparametric coordinates agree
        for(int i=0; i< DIM; ++i)
        {
          EXPECT_NEAR(isoparCenter[i],foundIsoPar1[i], EPS);
        }

      }
    }
    SLIC_INFO(fmt::format("Verifying {} pts on curved quad jittered mesh took {} s -- rate: {} q/s",
        pts.size() * mesh2.GetNE() * 2,
        queryTimer2.elapsed(),
        pts.size() * mesh2.GetNE() * 2 / queryTimer2.elapsed()  ) );

    //spatialIndex1.printDebugMesh( filename + "001.vtk");
    //spatialIndex2.printDebugMesh( filename + "002.vtk");


}


TEST(quest_point_in_cell, pic_curved_quad_c_shaped_output_mesh)
{
    using namespace types2D;
    
    const std::string feCollName = "Quadratic";
    
    SLIC_INFO("Generating diagnostic mesh for"
        << " C-shaped biquadratic element");


    // Here we are testing a very curved C-shaped mesh
    // We are comparing a single element mesh to
    // a refined mesh, and checking that the PointInCell query
    // gives the same result for both (i.e. both inside or both outside)
    
    std::stringstream sstr1;
    sstr1 << meshPrefix
         << fmt::format(c_shapedNodes, feCollName);
    mfem::Mesh mesh1(sstr1);

    const int res = 25;

    PointInCellType spatialIndex1(&mesh1, GridCell(res));
    
    int ext[4] = { 0,res,0,res};
    axom::mint::CurvilinearMesh cmesh(2, ext);

    {
      axom::primal::Point<int,3> ext_size;
      cmesh.getExtentSize( ext_size.data() );
      SLIC_INFO( "Extents of curvilinear mesh: " << ext_size    );
    }
  
    // Set the positions of the nodes of the diagnostic mesh
    const double denom = res;
    for(int i=0; i <= res; ++i)
    {
        for(int j=0; j<= res; ++j)
        {
            SpacePt pt = spatialIndex1.findInSpace(0, SpacePt::make_point(i/denom,j/denom) );
            cmesh.setNode(i,j,pt[0],pt[1]);
        }
    }

    // Add a scalar field on the cells -- value is 1 when isoparametric transform succeeds, 0 otherwise
    //if(false)
    {
        std::string name = "query_status";
        axom::mint::FieldData* CD = cmesh.getCellFieldData();     
        int numSuccesses = 0;
        const int numCells =   cmesh.getMeshNumberOfCells();
        SLIC_INFO("Mesh has " << numCells << " cells.");
        CD->addField( new axom::mint::FieldVariable< int >(name, numCells ) );
        int* fld = CD->getField( name )->getIntPtr();

        for(int i=0; i < res; ++i)
        {
            for(int j=0; j< res; ++j)
            {
                // Forward map
                double midX = (2.*i+1)/(2.*denom);
                double midY = (2.*j+1)/(2.*denom);
                
                SpacePt pt = spatialIndex1.findInSpace(0, SpacePt::make_point(midX, midY) );

                // Reverse map
                SpacePt isoPt;
                bool found = spatialIndex1.getIsoparametricCoords(0, pt, &isoPt);

                int idx = cmesh.getCellLinearIndex(i,j);
                fld[idx] = found ? 1 : -1;

                if(found) { ++numSuccesses; }
            }
        }

        SLIC_INFO(fmt::format("Found {} of {} points ({}%)",
            numSuccesses, numCells, (100. * numSuccesses)/numCells ) );

    }      

    std::stringstream filenameStr;
    filenameStr << "quest_point_in_cell_c_shaped_quad_001_mint_" << res << ".vtk";
    SLIC_INFO("About to write file " << filenameStr.str());
    axom::mint::write_vtk(&cmesh, filenameStr.str() );

    std::string filename = "quest_point_in_cell_c_shaped_quad";
    spatialIndex1.printDebugMesh( filename + "_001.vtk");

}

TEST(quest_point_in_cell, printIsoparams)
{
  const int geom = mfem::Geometry::SQUARE;
  const int dim = 2;

  for(int order = 0; order < 3; ++order)
  {
    const mfem::IntegrationRule* ir = &(mfem::RefinedIntRules.Get(geom, order ) );

    const int npts = ir->GetNPoints();
    for(int i=0; i < npts; ++i)
    {
        double x[3];
        ir->IntPoint(i).Get(x, dim);
        SLIC_INFO( "integration point  " << i
            << "\n\t x: " << x[0]
            << "\n\t y: " << x[1]
            << "\n\t order " << order
            );
    }
  }
}




int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,
  axom::slic::setLoggingMsgLevel( axom::slic::message::Info );

  std::srand( SRAND_SEED );

  result = RUN_ALL_TESTS();
  return result;
}
