// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MIR_MESH_TEST_H_
#define MIR_MESH_TEST_H_

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/mir.hpp"


namespace mir = axom::mir;

class MirMeshTest : public ::testing::Test
{
public:

   template<typename T>
   using Vec = std::vector<T>;

   using IndexVec = Vec<mir::PosType>;
   using VolFracVec = Vec<axom::float64>;
   using VolumeFractions = Vec<VolFracVec>;

   enum { GREEN = 0, BLUE = 1 };

public:

   mir::MIRMesh& getMesh() { return m_mesh; }
   const mir::MIRMesh& getMesh() const { return m_mesh; }

protected:
  void SetUp() override
  {
     m_verts= mir::VertSet(16);
     m_elems= mir::ElemSet(9);
     m_nMats = 2;

     // Create the mesh connectivity information
     setupTopoData();
     setupMapData();
     setupVolumeFractions();

     m_mesh.initializeMesh(this->m_verts, this->m_elems, this->m_nMats,
           this->m_topoData, this->m_mapData, this->m_volFracs);
  }

  // Set up the mesh's topological connectivity
  void setupTopoData()
  {
     m_topoData.m_evInds = {
         0,4,5,1,     // elem 0, card 4, start 0
         1,5,6,2,     // elem 1, card 4, start 4
         2,6,7,3,     // elem 2, card 4, start 8
         4,8,9,5,     // elem 3, card 4, start 12
         5,9,10,6,    // elem 4, card 4, start 16
         6,10,11,7,   // elem 5, card 4, start 20
         8,12,13,9,   // elem 6, card 4, start 24
         9,13,14,10,  // elem 7, card 4, start 28
         10,14,15,11  // elem 8, card 4, start 32, end 36
       };

     m_topoData.m_evBegins = {
         0,4,8,12,16,20,24,28,32,36
       };

     m_topoData.m_veInds = {
         0,          // vert  0, card 1, start 0
         0,1,        // vert  1, card 2, start 1
         1,2,        // vert  2, card 2, start 3
         2,          // vert  3, card 1, start 5
         0,3,        // vert  4, card 2, start 6
         0,1,3,4,    // vert  5, card 4, start 8
         1,2,4,5,    // vert  6, card 4, start 12
         2,5,        // vert  7, card 2, start 16
         3,6,        // vert  8, card 2, start 18
         3,4,6,7,    // vert  9, card 4, start 20
         4,5,7,8,    // vert  10, card 4, start 24
         5,8,        // vert  11, card 2, start 28
         6,          // vert  12, card 1, start 30
         6,7,        // vert  13, card 2, start 31
         7,8,        // vert  14, card 2, start 33
         8,          // vert  15, card 1, start 35, end 36
       };

     m_topoData.m_veBegins = {
         0,1,3,5,6,8,12,16,18,20,24,28,30,31,33,35,36
       };
  }

  // Set up the mesh's map data
  void setupMapData()
  {
     m_mapData.m_vertexPositions = {
       mir::Point2( 0.0, 3.0 ),
       mir::Point2( 1.0, 3.0 ),
       mir::Point2( 2.0, 3.0 ),
       mir::Point2( 3.0, 3.0 ),

       mir::Point2( 0.0, 2.0 ),
       mir::Point2( 1.0, 2.0 ),
       mir::Point2( 2.0, 2.0 ),
       mir::Point2( 3.0, 2.0 ),

       mir::Point2( 0.0, 1.0 ),
       mir::Point2( 1.0, 1.0 ),
       mir::Point2( 2.0, 1.0 ),
       mir::Point2( 3.0, 1.0 ),

       mir::Point2( 0.0, 0.0 ),
       mir::Point2( 1.0, 0.0 ),
       mir::Point2( 2.0, 0.0 ),
       mir::Point2( 3.0, 0.0 )
     };

     m_mapData.m_elementDominantMaterials = Vec<int>(m_elems.size(), mir::NULL_MAT);
     m_mapData.m_elementParents = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
     m_mapData.m_shapeTypes = Vec<mir::Shape>(m_elems.size(), mir::Shape::Quad);
  }

  // Set up the mesh's volume fraction data
  void setupVolumeFractions()
  {
     m_volFracs.resize(m_nMats);
     m_volFracs[GREEN] = {1.0, 1.0, 1.0, 1.0, 0.5, 0.2, 0.2, 0.0, 0.0};
     m_volFracs[BLUE]  = {0.0, 0.0, 0.0, 0.0, 0.5, 0.8, 0.8, 1.0, 1.0};
  }

protected:
  mir::VertSet m_verts;
  mir::ElemSet m_elems;
  int m_nMats;

  mir::CellTopologyData m_topoData;
  mir::CellMapData m_mapData;
  mir::CellData m_cellData;
  VolumeFractions m_volFracs;

  mir::MIRMesh m_mesh;
};


TEST_F(MirMeshTest,default_ctor)
{
   mir::MIRMesh mesh;
   EXPECT_TRUE(mesh.isValid(true));
}


TEST_F(MirMeshTest,initialize)
{
   mir::MIRMesh mesh;
   mesh.initializeMesh(
         this->m_verts,
         this->m_elems,
         this->m_nMats,
         this->m_topoData,
         this->m_mapData,
         this->m_volFracs);

   EXPECT_TRUE(mesh.isValid(true));

   mesh.print();
}

TEST_F(MirMeshTest,copy_ctor)
{
   // test copy constructor
   {
     mir::MIRMesh& mesh = this->getMesh();
     EXPECT_TRUE(mesh.isValid(true));

     mir::MIRMesh mirCopy(mesh);
     EXPECT_TRUE( mirCopy.isValid(true) );
   }

   // test const copy constructor
   {
     const mir::MIRMesh& cmesh = this->getMesh();
     EXPECT_TRUE(cmesh.isValid(true));

     // test copy constructor
     const mir::MIRMesh mirCopy(cmesh);
     EXPECT_TRUE( mirCopy.isValid(true) );
   }
}

TEST_F(MirMeshTest,copy_assign)
{
   // test copy assignment from a non-const MIRMesh
   {
     mir::MIRMesh& mesh = this->getMesh();
     EXPECT_TRUE(mesh.isValid(true));

     mir::MIRMesh mirCopy;
     mirCopy = mesh;
     EXPECT_TRUE( mirCopy.isValid(true) );
   }

   // test copy assignment from a const MIRMesh
   {
     const mir::MIRMesh& cmesh = this->getMesh();
     EXPECT_TRUE(cmesh.isValid(true));

     mir::MIRMesh mirCopy;
     mirCopy = cmesh;
     EXPECT_TRUE( mirCopy.isValid(true) );
   }
}


//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::UnitTestLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();
  return result;
}


#endif //  MIR_MESH_TEST_H_
