/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */



#include "gtest/gtest.h"

#include "axom_utils/Timer.hpp"


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
      pts.push_back( SpacePt2D::make_point( static_cast<double>(i)/res, 
                                            static_cast<double>(j)/res) );
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
void testMesh(mfem::Mesh* mesh, ExpectedValueFunctor exp, const std::string& meshTypeStr)
{
  using namespace types2D;

  // output mesh in mfem and vtk format
  // output refined mesh to logger in mfem and vtk formats
  bool isCurved = mesh->GetNodes() != AXOM_NULLPTR;
  std::string filename = fmt::format("quest_point_in_cell_{}_quad",meshTypeStr);
  {
    mfem::VisItDataCollection dataCol(filename, mesh);
    if(mesh->GetNodes())
      dataCol.RegisterField("nodes", mesh->GetNodes());
    dataCol.Save();
  }
//  {
//    std::ofstream mfem_stream(filename + ".vtk");
//    mesh->PrintVTK(mfem_stream);
//  }

  // Add mesh to the grid
  axom::utilities::Timer constructTimer(true);
  PointInCellType spatialIndex(mesh, GridCell(25));
  SLIC_INFO(fmt::format("Constructing index over {} quad mesh with {} elems took {} s",
      meshTypeStr,
      mesh->GetNE(),
      constructTimer.elapsed() ) );

  std::vector<SpacePt> pts = generate2DTestPoints( exp.radius() );
  SpacePt isoPar;
  int numTested = 0;

  axom::utilities::Timer queryTimer(true);
  for( std::vector<SpacePt>::iterator it= pts.begin(); it != pts.end(); ++it)
  {
    SpacePt& queryPoint = *it;

    int idx = spatialIndex.locatePoint( queryPoint, &isoPar );
    bool isInMesh = (idx != PointInCellType::NO_CELL);

    if( exp.canTestPoint(queryPoint) )
    {
      numTested++;

      bool expectedInMesh = exp.expectedInMesh(queryPoint);
      EXPECT_EQ( expectedInMesh, isInMesh) << "Point " << *it;
    }

    if(isInMesh)
    {
      SpacePt untransformPt = spatialIndex.findInSpace(idx, isoPar);

      for(int i=0; i< DIM; ++i)
      {
        EXPECT_NEAR(queryPoint[i],untransformPt[i], EPS);
      }
    }
  }


  SLIC_INFO(fmt::format("Querying {} pts on {} quad mesh took {} s -- rate: {} q/s",
      pts.size(),
      meshTypeStr,
      queryTimer.elapsed(),
      pts.size() / queryTimer.elapsed()  ) );

  if(isCurved)
  {
    SLIC_INFO(fmt::format("On {} mesh, verified plausibility in {} of {} cases ({}%)",
        meshTypeStr, numTested, pts.size(), static_cast<double>(numTested)/ pts.size() ));
  }


  // Test that a fixed set of isoparametric coords on each cell maps to the correct place.
  pts = generate2DIsoParTestPoints(10);
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
  int expectedNE = std::pow<int>(ELT_MULT_FAC, NREFINE);
  EXPECT_EQ( expectedNE, mesh.GetNE() );
  EXPECT_EQ(AXOM_NULLPTR, mesh.GetNodes());

  SCOPED_TRACE("point_in_cell_flat");
  std::string meshTypeStr = fmt::format("flat_{}", NREFINE);;
  testMesh(&mesh, ExpectedValueFlat(VAL), meshTypeStr);

}

TEST( quest_point_in_cell, pic_curved_quad )
{
  using namespace types2D;

  const double VERT_VAL = 0.5;
  const double DIAG_VAL = VERT_VAL * std::sqrt(2.)/2.;
  const std::string feCollName = "Quadratic";

  // Create a simple high order mfem mesh over a single quad
  //  -- a diamond with verts at +-(VAL,0) and +-(0, VAL) and nodes along diagonals (+- VAL2, +-VAL2}
  std::stringstream sstr;
  sstr << meshPrefix
       << fmt::format(highOrderNodes, VERT_VAL, DIAG_VAL, feCollName);

  mfem::Mesh mesh(sstr);

  // Test single element mesh
  {
    SCOPED_TRACE("point_in_cell_curved");
    std::string meshTypeStr = fmt::format("{}_{}", feCollName, 1);
    testMesh(&mesh, ExpectedValueCurved(VERT_VAL), meshTypeStr);
  }

  for(int i=0; i< NREFINE; ++i)
    mesh.UniformRefinement();

  int expectedNE = std::pow<int>(ELT_MULT_FAC, NREFINE);
  EXPECT_EQ( expectedNE, mesh.GetNE() );
  EXPECT_NE(AXOM_NULLPTR, mesh.GetNodes());
  EXPECT_EQ(feCollName, mesh.GetNodalFESpace()->FEColl()->Name());

  // Test refined mesh
  {
    SCOPED_TRACE("point_in_cell_curved");
    std::string meshTypeStr = fmt::format("{}_{}", feCollName, NREFINE);
    testMesh(&mesh, ExpectedValueCurved(VERT_VAL), meshTypeStr);
  }
}

TEST( quest_point_in_cell, pic_curved_quad_jittered )
{
  using namespace types2D;

  const double VERT_VAL = 0.5;
  const double DIAG_VAL = VERT_VAL * std::sqrt(2.)/2.;
  const std::string feCollName = "Quadratic";
  const double jitterFactor = .15;

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

    SCOPED_TRACE("point_in_cell_curved_jittered_1");
    std::string meshTypeStr = fmt::format("{}_jittered_{}", feCollName, 1);
    testMesh(&mesh, ExpectedValueCurved(VERT_VAL), meshTypeStr);
  }

  for(int i=0; i< NREFINE -1; ++i)
  {
    mesh.UniformRefinement();
    jitterNodalValues<DIM>(&mesh, jitterFactor);
  }

  SCOPED_TRACE("point_in_cell_curved_jittered");
  std::string meshTypeStr = fmt::format("{}_jittered_{}", feCollName, NREFINE);
  testMesh(&mesh, ExpectedValueCurved(VERT_VAL), meshTypeStr);

}


TEST( quest_point_in_cell, pic_curved_quad_jittered_pos )
{
  using namespace types2D;

  const double VERT_VAL = 0.5;
  const double DIAG_VAL = VERT_VAL * (std::sqrt(2.) - 0.5);
  const std::string feCollName = "QuadraticPos";
  const double jitterFactor = .25;

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

    SCOPED_TRACE("point_in_cell_curved_jittered_positive_1");
    std::string meshTypeStr = fmt::format("{}_jittered_{}", feCollName, 1);
    testMesh(&mesh, ExpectedValueCurved(VERT_VAL), meshTypeStr);
  }

  for(int i=0; i< NREFINE -1; ++i)
  {
    mesh.UniformRefinement();
    jitterNodalValues<DIM>(&mesh, jitterFactor);
  }

  SCOPED_TRACE("point_in_cell_curved_jittered_positive");
  std::string meshTypeStr = fmt::format("{}_jittered_{}", feCollName, NREFINE);
  testMesh(&mesh, ExpectedValueCurved(VERT_VAL), meshTypeStr);

}


TEST( quest_point_in_cell, pic_curved_quad_c_shaped )
{
  using namespace types2D;

  const std::string feCollName = "Quadratic";

  // Here we are testing a very curved C-shaped mesh
  // We are comparing a single element mesh to
  // a refined mesh, and checking that the PointInCell query
  // gives the same result for both (i.e. both inside or both outside)

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
    mfem::VisItDataCollection dataCol(filename + "001", &mesh1);
    dataCol.Save();
  }
  {
    mfem::VisItDataCollection dataCol(filename + "002", &mesh2);
    dataCol.Save();
  }

  SLIC_INFO("Characteristics for the single element quadratic quad mesh:");
  mesh1.PrintCharacteristics();
  SLIC_INFO("\n-- Elt size sqrt det: " << mesh1.GetElementSize(0, 0) );
  SLIC_INFO("\n-- Elt size h_min: " << mesh1.GetElementSize(0, 1) );
  SLIC_INFO("\n-- Elt size h_max: " << mesh1.GetElementSize(0, 2) );

  if(true )
    return;


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

  // Test that both queries agree
  axom::utilities::Timer queryTimer(true);
  const int num_pts = 100;  // NUM_TEST_PTS;
  for( int i=0; i< num_pts; ++i)
  {
    SpacePt queryPoint = axom::quest::utilities::randomSpacePt<DIM>(-10, 10);

    SpacePt isoPar1;
    int idx1 = spatialIndex1.locatePoint( queryPoint, &isoPar1);
    bool isInMesh1 = (idx1 != PointInCellType::NO_CELL);

    SpacePt isoPar2;
    int idx2 = spatialIndex2.locatePoint( queryPoint, &isoPar2);
    bool isInMesh2 = (idx2 != PointInCellType::NO_CELL);

    EXPECT_EQ(isInMesh1, isInMesh2)
        << "Point " << queryPoint
        << " was " << (isInMesh1 ? "in" : "not in") << " mesh 1"
        << " and was " << (isInMesh2 ? "in" : "not in") << " mesh 2"
        << " They should agree.";

    if(isInMesh1)
    {
      SpacePt untransformPt = spatialIndex1.findInSpace(idx1, isoPar1);

      if(! isInMesh2)
      {
        SLIC_INFO("UntransformedPt for mesh 1 was " << untransformPt);
      }
//      for(int i=0; i< DIM; ++i)
//      {
//        EXPECT_NEAR(queryPoint[i],untransformPt[i], EPS);
//      }

    }
    if(isInMesh2)
    {
      SpacePt untransformPt = spatialIndex2.findInSpace(idx2, isoPar2);

      if(! isInMesh1)
      {
        SLIC_INFO("UntransformedPt for mesh 2 was " << untransformPt);
      }

//      for(int i=0; i< DIM; ++i)
//      {
//        EXPECT_NEAR(queryPoint[i],untransformPt[i], EPS);
//      }

    }
  }

  SLIC_INFO(fmt::format("Querying {} pts on C-shaped quadratic quad mesh took {} s -- rate: {} q/s",
      num_pts * 2,
      queryTimer.elapsed(),
      num_pts * 2 / queryTimer.elapsed()  ) );



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
        //EXPECT_EQ( eltId, foundCellId)
        //    << fmt::format("For element {} -- computed midpoint {} @ {} -- found isoPar is {}",
        //        eltId, spacePt, isoparCenter, foundIsoPar);

        // Check that the isoparametric coordinates agree
        //for(int i=0; i< DIM; ++i)
        //{
        //  EXPECT_NEAR(isoparCenter[i],foundIsoPar[i], EPS);
        //}

      }
    }
    SLIC_INFO(fmt::format("Verifying {} pts on curved quad jittered mesh took {} s -- rate: {} q/s",
        pts.size() * mesh2.GetNE() * 2,
        queryTimer2.elapsed(),
        pts.size() * mesh2.GetNE() * 2 / queryTimer2.elapsed()  ) );



}




int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,
  axom::slic::setLoggingMsgLevel( axom::slic::message::Info );

  // finalized when exiting main scope

  std::srand( SRAND_SEED );

  result = RUN_ALL_TESTS();

  return result;
}
