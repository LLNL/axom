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

  const double EPS = 1e-8;

#ifdef AXOM_DEBUG
  const int NREFINE = 5;
  const int NUM_TEST_PTS = 10000;
#else
  const int NREFINE = 5;
  const int NUM_TEST_PTS = 100000;
#endif

}

enum MeshType {
  FLAT_MESH,
  QUADRATIC_MESH,
  QUADRATIC_POS_MESH,
  C_SHAPED_MESH
};


class PointInCell2DTest : public ::testing::Test
{
public:
  static const int DIM = 2;
  static const int ELT_MULT_FAC = 4;

  typedef axom::primal::BoundingBox<double, DIM> BBox;
  typedef axom::primal::Point<double,DIM> SpacePt;
  typedef axom::primal::Vector<double,DIM> SpaceVec;
  typedef axom::primal::Point<int,DIM> GridCell;
  typedef axom::quest::PointInCell<DIM> PointInCellType;

protected:
  virtual void SetUp()
  {
    m_mesh = AXOM_NULLPTR;

    /// Setup mesh strings

    // Prefix string for a single element 2D mfem quad mesh
    m_meshPrefixStr =
        "MFEM mesh v1.0"  "\n\n"
        "dimension"         "\n"
        "2"               "\n\n"
        "elements"          "\n"
        "1"                 "\n"
        "1 3 0 1 2 3"     "\n\n"
        "boundary"          "\n"
        "0"               "\n\n";

    // Vertex positions for a single element quad mesh
    //   -- a diamond with verts at +-(VAL0,0) and +-(0, VAL0)
    // This is a circle with radius VAL in the L1 norm
    // Requires one value to be passed in for string interpolation
    m_lowOrderVertsStr =
        "vertices"       "\n"
        "4"              "\n"
        "2"              "\n"
        "0     -{0}"     "\n"
        "{0}      0"     "\n"
        "0      {0}"     "\n"
        "-{0}     0"     "\n";

    // Nodal grid function for a single element quadratic quad mesh
    //  -- a diamond with verts at +-(VAL0,0) and +-(0, VAL0)
    //  -- and nodes along diagonals (+- VAL1, +-VAL1}
    // Requires two values to be passed in for string interpolation
    m_highOrderNodesStr =
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

    // Nodal grid function for a C-shaped quadratic quadrilateral
    m_CShapedNodesStr =
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
        "-6 4"                          "\n"
        "0 7"                           "\n"
        "-8 4"                          "\n"
        "-7 4"                          "\n";
  }

  virtual void TearDown()
  {
    if(m_mesh != AXOM_NULLPTR)
    {
      delete m_mesh;
      m_mesh = AXOM_NULLPTR;
    }
  }

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

  /*!
   * Jitter all degrees of freedom (dofs) of the mesh's nodal grid function
   * Implementation borrowed from mfem's mesh-explorer mini-app
   */
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

public:
  /*!
   * Sets up a test mesh of the desired type with the specified parameters
   *
   * \param meshType Enum of type MeshType
   * \param numRefine The number of times to uniformly refine the mesh, default: 0
   * \param vertVal Parameter value for the mesh, default: 0.5
   * \param jitterFactor Factor by which to jitter the nodes, default: 0
   */
  void setupTestMesh(MeshType meshType, int numRefine=0, double vertVal = 0.5, double jitterFactor = 0.)
  {
    // Compose mesh input string
    std::stringstream sstr;
    sstr << m_meshPrefixStr;

    std::string feCollName;

    std::stringstream meshDescSstr; // used to compose the mesh descriptor string

    switch(meshType)
    {
    case FLAT_MESH:
      feCollName = "Linear";
      meshDescSstr << "flat";
      sstr << fmt::format(m_lowOrderVertsStr, vertVal);
      break;
    case QUADRATIC_MESH:
      feCollName = "Quadratic";
      meshDescSstr << "curved";

      sstr << fmt::format(m_highOrderNodesStr, vertVal, vertVal * std::sqrt(2.)/2., feCollName);
      break;
    case QUADRATIC_POS_MESH:
      feCollName = "QuadraticPos";
      meshDescSstr << "curved_pos";

      sstr << fmt::format(m_highOrderNodesStr, vertVal, vertVal * (std::sqrt(2.) - 0.5), feCollName);
      break;
    case C_SHAPED_MESH:
      feCollName = "Quadratic";
      meshDescSstr << "c_shaped";

      sstr << fmt::format(m_CShapedNodesStr, feCollName);
      break;
    default:
      FAIL() << "Did not provide a valid MeshType.";
      break;
    }

    // compose the mesh descriptor string
    {
      if(numRefine > 0)
        meshDescSstr << "_refined_" << numRefine;
      else
        meshDescSstr << "_single";

      if(jitterFactor > 0)
      {
        meshDescSstr << "_jittered";
      }

      m_meshDescriptorStr = meshDescSstr.str();

      SLIC_INFO(fmt::format("Generating {} mfem quad mesh", feCollName)
          << (jitterFactor > 0 ? fmt::format(" with jitter factor {}", jitterFactor) : "" )
          << (numRefine > 0 ? fmt::format(" refined to level {}.", numRefine) : "")
          << "\nDescriptor string: " << m_meshDescriptorStr);
    }

    m_mesh = new mfem::Mesh(sstr);
    EXPECT_NE(AXOM_NULLPTR, m_mesh);

    // Refine (and possibly jitter) the mesh several times
    for(int i=0; i< numRefine; ++i)
    {
      m_mesh->UniformRefinement();

      if(jitterFactor > 0)
      {
        jitterNodalValues(m_mesh, jitterFactor);
      }
    }

    // Sanity checks on number of elements
    int expectedNE = expectedNumElts(ELT_MULT_FAC, numRefine);
    EXPECT_EQ( expectedNE, m_mesh->GetNE() );

    // Sanity checks on whether mesh is low or high order
    if(meshType == FLAT_MESH)
    {
      EXPECT_EQ(AXOM_NULLPTR, m_mesh->GetNodes());
    }
    else
    {
      EXPECT_NE(AXOM_NULLPTR, m_mesh->GetNodes());
      EXPECT_EQ(feCollName, m_mesh->GetNodalFESpace()->FEColl()->Name());
    }

  }

  template<typename ExpectedValueFunctor>
  void testRandomPointsOnMesh(ExpectedValueFunctor exp, const std::string& meshTypeStr)
  {
    // output mesh in mfem and vtk format
    bool isCurved = m_mesh->GetNodes() != AXOM_NULLPTR;
    std::string filename = fmt::format("quest_point_in_cell_{}_quad",meshTypeStr);

    if(outputMeshMFEM)
    {
      mfem::VisItDataCollection dataCol(filename, m_mesh);
      if(m_mesh->GetNodes())
        dataCol.RegisterField("nodes", m_mesh->GetNodes());
      dataCol.Save();
    }
    if(outputMeshVTK)
    {
      std::ofstream mfem_stream(filename + ".vtk");
      m_mesh->PrintVTK(mfem_stream);
    }

    // Add mesh to the grid
    axom::utilities::Timer constructTimer(true);
    PointInCellType spatialIndex(m_mesh, GridCell(25));
    SLIC_INFO(fmt::format("Constructing index over {} quad mesh with {} elems took {} s",
        meshTypeStr,
        m_mesh->GetNE(),
        constructTimer.elapsed() ) );

    std::vector<SpacePt> pts = generateRandomTestPoints( exp.radius() );
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
          EXPECT_NEAR(queryPoint[i],untransformPt[i], ::EPS);
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

  void testIsoGridPointsOnMesh(const std::string& meshTypeStr)
  {
    std::string filename = fmt::format("quest_point_in_cell_{}_quad",meshTypeStr);

    // Add mesh to the grid
    axom::utilities::Timer constructTimer(true);
    PointInCellType spatialIndex(m_mesh, GridCell(25));
    SLIC_INFO(fmt::format("Constructing index over {} quad mesh with {} elems took {} s",
        meshTypeStr,
        m_mesh->GetNE(),
        constructTimer.elapsed() ) );


    // Test that a fixed set of isoparametric coords on each cell maps to the correct place.
    std::vector<SpacePt> pts = generateIsoParTestPoints(10);
    axom::utilities::Timer queryTimer2(true);
    SpacePt foundIsoPar;
    for(int eltId=0; eltId< m_mesh->GetNE(); ++eltId)
    {
      for( std::vector<SpacePt>::iterator it= pts.begin(); it != pts.end(); ++it)
      {
        SpacePt& isoparCenter = *it;

        // Check if isoparCenter is on element boundary
        bool isBdry = false;
        for(int i=0; i< DIM; ++i)
        {
          if( axom::utilities::isNearlyEqual(isoparCenter[i], 0.)
            || axom::utilities::isNearlyEqual(isoparCenter[i], 1.))
          {
            isBdry = true;
          }
        }

        SpacePt spacePt = spatialIndex.findInSpace(eltId, isoparCenter);

        int foundCellId = spatialIndex.locatePoint(spacePt, &foundIsoPar);

        // Check that we found a cell
        EXPECT_NE( PointInCellType::NO_CELL, foundCellId)
          << "element: " << eltId
          << " -- isopar: " << isoparCenter
          << " -- spacePt: " << spacePt
          << " -- isBdry: " << (isBdry ? "yes" : "no");

        if( !isBdry)
        {
          EXPECT_EQ( eltId, foundCellId)
              << fmt::format("For element {} -- computed space point {} from isoPar {} -- found isoPar is {}",
                  eltId, spacePt, isoparCenter, foundIsoPar);
        }

        // If we found the same cell, check that found isoparametric coordinates agree with original
        if(eltId == foundCellId)
        {
          for(int i=0; i< DIM; ++i)
          {
            EXPECT_NEAR(isoparCenter[i],foundIsoPar[i], ::EPS);
          }
        }

        // Convert point back to space, and check that it matches our original point
        if(foundCellId != PointInCellType::NO_CELL)
        {
          SpacePt transformedPt = spatialIndex.findInSpace(foundCellId, foundIsoPar);

          for(int i=0; i< DIM; ++i)
          {
            EXPECT_NEAR(spacePt[i],transformedPt[i], ::EPS);
          }
        }

      }
    }

    SLIC_INFO(fmt::format("Verifying {} pts on {} quad mesh took {} s -- rate: {} q/s",
        pts.size() * m_mesh->GetNE(),
        meshTypeStr,
        queryTimer2.elapsed(),
        pts.size() * m_mesh->GetNE() / queryTimer2.elapsed()  ) );

  }

  std::vector< SpacePt > generateRandomTestPoints(double val)
  {
    const int SZ = ::NUM_TEST_PTS;

    std::vector<SpacePt> pts;
    pts.reserve( SZ + 6 );

    pts.push_back( SpacePt::zero() );
    pts.push_back( SpacePt::ones() );
    pts.push_back( SpacePt::make_point(.1, .1) );
    pts.push_back( SpacePt::make_point(.4, .4) );
    pts.push_back( SpacePt::make_point(val, 0) );
    pts.push_back( SpacePt::make_point(val  +0.005, 0) );


    for(int i=0; i< SZ; ++i)
    {
      pts.push_back( axom::quest::utilities::randomSpacePt<DIM>(-1.25 * val, 1.25 * val));
    }

    return pts;
  }

  std::vector< SpacePt > generateIsoParTestPoints(int res)
  {
    std::vector<SpacePt> pts;

    for(int i=0; i <= res; ++i)
      for(int j=0; j <= res; ++j)
      {
        // Get the corresponding isoparametric value
        SpacePt pt = SpacePt::make_point( static_cast<double>(i)/res,
                                          static_cast<double>(j)/res);

//        // Add a small random offset
//        SpacePt off = axom::quest::utilities::randomSpacePt<DIM>( -::EPS, +::EPS);
//
//        // Clamp to ensure isoparametric pt is in element
//        for(int i=0; i< SpacePt::DIMENSION; ++i)
//          pt[i] = axom::utilities::clampVal(pt[i]+off[i], 0., 1.);


        pts.push_back( pt );
      }

    return pts;
  }

  mfem::Mesh* getMesh() { return m_mesh; }

  const std::string& getMeshDescriptor() const { return m_meshDescriptorStr; }

private:
  std::string m_meshPrefixStr;
  std::string m_lowOrderVertsStr;
  std::string m_highOrderNodesStr;
  std::string m_CShapedNodesStr;

  std::string m_meshDescriptorStr;

  mfem::Mesh* m_mesh;
};


template<int DIM>
struct ExpectedValueFlat
{
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
template<int DIM>
struct ExpectedValueCurved
{
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

TEST_F( PointInCell2DTest, pic_flat_single_quad )
{
  const double  vertVal = 0.5;
  const int numRefine = 0;

  this->setupTestMesh(FLAT_MESH, numRefine, vertVal);

  std::string meshTypeStr = this->getMeshDescriptor();
  SCOPED_TRACE(fmt::format("point_in_cell_{}",meshTypeStr));

  // Add a bilinear gridfunction
  mfem::Mesh& mesh = *this->getMesh();
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

  this->testRandomPointsOnMesh( ExpectedValueFlat<DIM>(vertVal), meshTypeStr);
  this->testIsoGridPointsOnMesh(meshTypeStr);
}

TEST_F( PointInCell2DTest, pic_flat_refined_quad )
{
  const double  vertVal = 0.5;
  const int numRefine = ::NREFINE;
  const int DIM = PointInCell2DTest::DIM;

  this->setupTestMesh(FLAT_MESH, numRefine, vertVal);

  std::string meshTypeStr = this->getMeshDescriptor();
  SCOPED_TRACE(fmt::format("point_in_cell_{}",meshTypeStr));

  this->testRandomPointsOnMesh( ExpectedValueFlat<DIM>(vertVal), meshTypeStr);
  this->testIsoGridPointsOnMesh(meshTypeStr);
}


TEST_F( PointInCell2DTest, pic_curved_single_quad )
{
  const double  vertVal = 0.5;
  const int numRefine = 0;
  const int DIM = PointInCell2DTest::DIM;

  this->setupTestMesh(QUADRATIC_MESH, numRefine, vertVal);

  std::string meshTypeStr = this->getMeshDescriptor();
  SCOPED_TRACE(fmt::format("point_in_cell_{}",meshTypeStr));

  this->testRandomPointsOnMesh(ExpectedValueCurved<DIM>(vertVal), meshTypeStr);
  this->testIsoGridPointsOnMesh(meshTypeStr);
}

TEST_F( PointInCell2DTest, pic_curved_refined_quad )
{
  const double  vertVal = 0.5;
  const int numRefine = ::NREFINE;
  const int DIM = PointInCell2DTest::DIM;

  this->setupTestMesh(QUADRATIC_MESH, numRefine, vertVal);

  std::string meshTypeStr = this->getMeshDescriptor();
  SCOPED_TRACE(fmt::format("point_in_cell_{}",meshTypeStr));

  this->testRandomPointsOnMesh(ExpectedValueCurved<DIM>(vertVal), meshTypeStr);
  this->testIsoGridPointsOnMesh(meshTypeStr);
}

TEST_F( PointInCell2DTest, pic_curved_single_quad_jittered )
{
  const double vertVal = 0.5;
  const double jitterFactor = .15;
  const int numRefine = 1;
  const int DIM = PointInCell2DTest::DIM;

  this->setupTestMesh(QUADRATIC_MESH, numRefine, vertVal, jitterFactor);

  std::string meshTypeStr = this->getMeshDescriptor();
  SCOPED_TRACE(fmt::format("point_in_cell_{}",meshTypeStr));

  this->testRandomPointsOnMesh(ExpectedValueCurved<DIM>(vertVal), meshTypeStr);
  this->testIsoGridPointsOnMesh(meshTypeStr);
}

TEST_F( PointInCell2DTest, pic_curved_refined_quad_jittered )
{
  const double vertVal = 0.5;
  const double jitterFactor = .15;
  const int numRefine = ::NREFINE;
  const int DIM = PointInCell2DTest::DIM;

  this->setupTestMesh(QUADRATIC_MESH, numRefine, vertVal, jitterFactor);

  std::string meshTypeStr = this->getMeshDescriptor();
  SCOPED_TRACE(fmt::format("point_in_cell_{}",meshTypeStr));

  this->testRandomPointsOnMesh(ExpectedValueCurved<DIM>(vertVal), meshTypeStr);
  this->testIsoGridPointsOnMesh(meshTypeStr);
}

TEST_F( PointInCell2DTest, pic_curved_single_quad_jittered_positive )
{
  const double vertVal = 0.5;
  const double jitterFactor = .15;
  const int numRefine = 1;
  const int DIM = PointInCell2DTest::DIM;

  this->setupTestMesh(QUADRATIC_POS_MESH, numRefine, vertVal, jitterFactor);

  std::string meshTypeStr = this->getMeshDescriptor();
  SCOPED_TRACE(fmt::format("point_in_cell_{}",meshTypeStr));

  this->testRandomPointsOnMesh(ExpectedValueCurved<DIM>(vertVal), meshTypeStr);
  this->testIsoGridPointsOnMesh(meshTypeStr);
}

TEST_F( PointInCell2DTest, pic_curved_refined_quad_jittered_positive )
{
  const double vertVal = 0.5;
  const double jitterFactor = .15;
  const int numRefine = ::NREFINE;
  const int DIM = PointInCell2DTest::DIM;

  this->setupTestMesh(QUADRATIC_POS_MESH, numRefine, vertVal, jitterFactor);

  std::string meshTypeStr = this->getMeshDescriptor();
  SCOPED_TRACE(fmt::format("point_in_cell_{}",meshTypeStr));

  this->testRandomPointsOnMesh(ExpectedValueCurved<DIM>(vertVal), meshTypeStr);
  this->testIsoGridPointsOnMesh(meshTypeStr);
}

TEST_F( PointInCell2DTest, pic_curved_quad_c_shaped )
{
  // Here we are testing a very curved C-shaped mesh
  // We are comparing a single element mesh to
  // a refined mesh (one level of refinement), and checking
  // that the PointInCell query gives the same result for both
  // (i.e. both inside or both outside)

  SLIC_INFO("Constructing and querying PointInCell structure"
      << " over C-shaped biquadratic element");


  this->setupTestMesh(C_SHAPED_MESH);

  // mesh1 contains a single C-shaped quadratic quad element
  mfem::Mesh& mesh1 = *this->getMesh();
  EXPECT_EQ(1, mesh1.GetNE() );

  // mesh2 is the same mesh with one additional refinement level
  mfem::Mesh mesh2( mesh1, true);
  mesh2.UniformRefinement();

  // output refined mesh to logger in mfem and vtk formats
  std::string filename = "quest_point_in_cell_c_shaped_quad";
  {
    mfem::VisItDataCollection dataCol(filename + "001_mesh", &mesh1);
    dataCol.Save();

    SLIC_INFO("Characteristics for the single element quadratic quad mesh:");
    mesh1.PrintCharacteristics();
    SLIC_INFO("\n-- Elt size sqrt det: " << mesh1.GetElementSize(0, 0) );
    SLIC_INFO("\n-- Elt size h_min: " << mesh1.GetElementSize(0, 1) );
    SLIC_INFO("\n-- Elt size h_max: " << mesh1.GetElementSize(0, 2) );
  }
  {
    mfem::VisItDataCollection dataCol(filename + "002_mesh", &mesh2);
    dataCol.Save();
  }

  // Create PointInCell structures over mesh1 and mesh2
  axom::utilities::Timer constructTimer(true);
  PointInCellType spatialIndex1(&mesh1, GridCell(10));
  SLIC_INFO(fmt::format("Constructing index over curved quad mesh1 with {} elems took {} s",
      mesh1.GetNE(),
      constructTimer.elapsed() ) );

  axom::utilities::Timer constructTimer2(true);
  PointInCellType spatialIndex2(&mesh2, GridCell(10));
  SLIC_INFO(fmt::format("Constructing index over curved quad mesh2 with {} elems took {} s",
      mesh2.GetNE(),
      constructTimer2.elapsed() ) );

  // Create random points; for each point check that both queries agree
  SLIC_INFO("Querying random points in domain of C-shaped mesh.");
  int numUntransformed = 0;
  axom::utilities::Timer queryTimer(true);
  const int num_pts = ::NUM_TEST_PTS;
  for( int i=0; i< num_pts; ++i)
  {
    // Create a random point in the domain of mesh1
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
    if(isInMesh1)
    {
      ++numUntransformed;
      SpacePt untransformPt = spatialIndex1.findInSpace(idx1, isoPar1);

      if(! isInMesh2)
      {
        SLIC_INFO("untransformed point for mesh 1 was " << untransformPt);
      }

      // Transformed point should match query point
      for(int i=0; i< DIM; ++i)
      {
        EXPECT_NEAR(queryPoint[i],untransformPt[i], ::EPS);
      }

    }

    // Try to transform back to space using spatialIndex2
    if(isInMesh2)
    {
      ++numUntransformed;
      SpacePt untransformPt = spatialIndex2.findInSpace(idx2, isoPar2);

      if(! isInMesh1)
      {
        SLIC_INFO("untransformed point for mesh 2 was " << untransformPt);
      }

      // Transformed point should match query point
      for(int i=0; i< DIM; ++i)
      {
        EXPECT_NEAR(queryPoint[i],untransformPt[i], ::EPS);
      }

    }
  }

  SLIC_INFO(fmt::format("Querying {} random pts on two C-shaped quadratic quad meshes took {} s -- rate: {} q/s",
      num_pts * 2,
      queryTimer.elapsed(),
      num_pts * 2 / queryTimer.elapsed()  )
    << "\n\t (includes " << numUntransformed << " transformations back into space)"
  );


  /// Test that a fixed set of isoparametric coords on each cell maps to the correct place.
  std::vector<SpacePt> pts = this->generateIsoParTestPoints(10);
  axom::utilities::Timer queryTimer2(true);
  SpacePt foundIsoPar1, foundIsoPar2;
  const int eltId=0; // Recall: mesh1 only has a single element

  for( std::vector<SpacePt>::iterator it= pts.begin(); it != pts.end(); ++it)
  {
    // Iteration point is our isoparametric coordinate
    SpacePt& isoparCenter = *it;

    // Find the corresponding point in space
    SpacePt spacePt = spatialIndex1.findInSpace(eltId, isoparCenter);

    // Check that we can find this point in mesh1
    int foundCellId1 = spatialIndex1.locatePoint(spacePt, &foundIsoPar1);
    EXPECT_NE(PointInCellType::NO_CELL, foundCellId1)
        << fmt::format("Pt {} was transformed from mesh1 using isoparametric coordinates {}.",
            spacePt, isoparCenter) << " Failed to reverse the transformation.";

    // Check that the isoparametric coordinates agree
    for(int i=0; i< DIM; ++i)
    {
      EXPECT_NEAR(isoparCenter[i],foundIsoPar1[i], EPS);
    }

    // Check that we can find this point in mesh2
    int foundCellId2 = spatialIndex2.locatePoint(spacePt, &foundIsoPar2);
    EXPECT_NE(PointInCellType::NO_CELL, foundCellId2)
        << fmt::format("Pt {} has isoparametric coordinates {} in mesh1, should be in mesh2",
            spacePt, isoparCenter);

  }

  SLIC_INFO(fmt::format("Verifying {} pts on curved quad jittered mesh took {} s -- rate: {} q/s",
      pts.size() * mesh2.GetNE() * 2,
      queryTimer2.elapsed(),
      pts.size() * mesh2.GetNE() * 2 / queryTimer2.elapsed()  ) );

    //spatialIndex1.printDebugMesh( filename + "001.vtk");
    //spatialIndex2.printDebugMesh( filename + "002.vtk");
}

TEST_F(PointInCell2DTest, pic_curved_quad_c_shaped_output_mesh)
{
    SLIC_INFO("Generating diagnostic mesh for"
        << " C-shaped biquadratic element");

    this->setupTestMesh(C_SHAPED_MESH);

    // mesh1 contains a single C-shaped quadratic quad element
    mfem::Mesh& mesh1 = *this->getMesh();
    
    const int res = 25;
    PointInCellType spatialIndex1(&mesh1, GridCell(res));
    
    // Setup linear mint mesh to approximate our mesh
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
            SpacePt isoparPt = SpacePt::make_point(i/denom,j/denom);
            SpacePt spacePt = spatialIndex1.findInSpace(0, isoparPt);
            cmesh.setNode(i,j,spacePt[0],spacePt[1]);
        }
    }

    // Add a scalar field on the cells
    // -- value is 1 when isoparametric transform succeeds, 0 otherwise
    {
      int numSuccesses = 0;
      const int numCells =  cmesh.getMeshNumberOfCells();
      SLIC_INFO("Mesh has " << numCells << " cells.");

      std::string name = "query_status";
      axom::mint::FieldData* CD = cmesh.getCellFieldData();
      CD->addField( new axom::mint::FieldVariable< int >(name, numCells ) );
      int* fld = CD->getField( name )->getIntPtr();

      for(int i=0; i < res; ++i)
      {
          for(int j=0; j< res; ++j)
          {
              // Forward map
              double midX = (2.*i+1)/(2.*denom);
              double midY = (2.*j+1)/(2.*denom);
              SpacePt origIsoPt = SpacePt::make_point(midX, midY);
              SpacePt pt = spatialIndex1.findInSpace(0,  origIsoPt);

              // Reverse map
              SpacePt isoPt;
              bool found = spatialIndex1.getIsoparametricCoords(0, pt, &isoPt);

              // Check that we were able to reverse the xform
              EXPECT_TRUE(found);
              for(int i=0; i< DIM; ++i)
              {
                EXPECT_NEAR(origIsoPt[i],isoPt[i], ::EPS);
              }

              int idx = cmesh.getCellLinearIndex(i,j);
              fld[idx] = found ? 1 : -1;

              if(found) { ++numSuccesses; }
          }
      }

        SLIC_INFO(fmt::format("Found {} of {} points ({}%)",
            numSuccesses, numCells, (100. * numSuccesses)/numCells ) );

    }      

    // Dump the mint mesh
    std::stringstream filenameStr;
    filenameStr << "quest_point_in_cell_c_shaped_quad_001_mint_" << res << ".vtk";
    SLIC_INFO("About to write file " << filenameStr.str());
    axom::mint::write_vtk(&cmesh, filenameStr.str() );

    // Dump the PIC debug mesh containing the Newton-Raphson paths
    // Note: This is debug code and will be removed
    std::string filename = "quest_point_in_cell_c_shaped_quad";
    spatialIndex1.printDebugMesh( filename + "_001.vtk");
}

TEST(quest_point_in_cell, printIsoparams)
{
  const int geom = mfem::Geometry::SQUARE;
  const int dim = 2;

  //for(int order=0; order < 5; ++order)
  const int order = 2;
  {
    const mfem::IntegrationRule* ir = &(mfem::RefinedIntRules.Get(geom, order ) );

    const int npts = ir->GetNPoints();
    SLIC_INFO("There are " << npts <<" refinement integration points of order " << order);

    double x[3];
    for(int i=0; i < npts; ++i)
    {
        ir->IntPoint(i).Get(x, dim);
        SLIC_INFO( "integration point  " << i
            << "\t x: " << x[0]
            << "\t y: " << x[1]
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
