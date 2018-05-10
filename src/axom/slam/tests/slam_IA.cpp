/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * \file slam_IA.cpp
 *
 * \brief ...
 */

#include "gtest/gtest.h"

#include "axom/config.hpp"    // for AXOM_USE_BOOST
#include "slic/slic.hpp"        // for SLIC_INFO

#include "slam/IA.hpp"
#include "slam/Utilities.hpp"


typedef axom::slam::util::Point3<double>  PointType;


struct BasicTriMeshData
{
  //This is a regular cube composed of triangles
  std::vector<double> points;
  std::vector<int> elem;

  static const double point_arr[];
  static const int tri_arr [];
  static const int vert_to_el_num []; // for verification purpose
  static const int el_nbr_rel[];

  BasicTriMeshData()
    : points (point_arr, point_arr + 24 )
    , elem(tri_arr, tri_arr + 36 )
  { }
};
const double BasicTriMeshData::point_arr [] = {
  0.0, 0.0, 0.0,
  0.0, 0.0, 1.0,
  0.0, 1.0, 0.0,
  0.0, 1.0, 1.0,
  1.0, 0.0, 0.0,
  1.0, 0.0, 1.0,
  1.0, 1.0, 0.0,
  1.0, 1.0, 1.0
};
const int BasicTriMeshData::tri_arr [] = {
  0, 6, 4,
  0, 2, 6,
  0, 3, 2,
  2, 7, 6,
  2, 3, 7,
  4, 6, 7,
  4, 7, 5,
  0, 4, 5,
  0, 1, 3,
  0, 5, 1,
  1, 5, 7,
  1, 7, 3
};
const int BasicTriMeshData::vert_to_el_num [] = {
  6, 4, 4, 4, 4, 4, 4, 6
};
const int BasicTriMeshData::el_nbr_rel [] = {
  1,5,7,
  2,3,0,
  8,4,1,
  4,5,1,
  2,11,3,
  0,3,6,
  5,10,7,
  0,6,9,
  9,11,2,
  7,10,8,
  9,6,11,
  10,4,8
};


struct BasicTetMeshData
{
  //cube divided into tets

  std::vector<double> points;
  std::vector<int> elem;

  static const double point_arr[];
  static const int tet_arr[];
  static const int vert_to_el_num [];
  static const int el_nbr_rel[];

  BasicTetMeshData()
    : points (point_arr, point_arr + 24 )
    , elem(tet_arr, tet_arr + 24 )
  { }

};
const double BasicTetMeshData::point_arr[] = {
  -1.0, -1.0, -1.0,
  -1.0, -1.0,  1.0,
  -1.0,  1.0, -1.0,
  -1.0,  1.0,  1.0,
  1.0, -1.0, -1.0,
  1.0, -1.0,  1.0,
  1.0,  1.0, -1.0,
  1.0,  1.0,  1.0
};
const int BasicTetMeshData::tet_arr [] = {
  3,2,4,0,
  3,1,4,0,
  3,6,2,4,
  3,6,7,4,
  3,5,1,4,
  3,5,7,4
};
const int BasicTetMeshData::vert_to_el_num [] = {
  2, 2, 2, 6, 6, 2, 2, 2
};
const int BasicTetMeshData::el_nbr_rel [] = {
  2,-1,1,-1,
  4,-1,0,-1,
  -1,-1,0,3,
  -1,-1,5,2,
  -1,-1,1,5,
  -1,-1,3,4
};

TEST(gtest_slam_IA, empty_mesh)
{
  SLIC_INFO("Testing creating an empty IA mesh...");

  typedef axom::slam::IAMesh<2,3, PointType>   IAMeshType;
  typedef IAMeshType::IndexListType IndexListType;

  IAMeshType ia_mesh; //empty mesh

  EXPECT_TRUE( ia_mesh.isValid(true) );
  EXPECT_EQ(ia_mesh.getNumberOfElements(), 0);
  EXPECT_EQ(ia_mesh.getNumberOfVertices(), 0);

  // How to check for SLIC_WARNINGS?
  IndexListType r;
  r = ia_mesh.getVerticesInElement(0);
  EXPECT_EQ((int)r.size(), 0);
  r = ia_mesh.getElementsWithVertex(0);
  EXPECT_EQ((int)r.size(), 0);
  r = ia_mesh.getElementNeighbors(0);
  EXPECT_EQ((int)r.size(), 0);

  //expect assert failure
  //EXPECT_EQ( PointType() , ia_mesh.getVertexPoint(0) );

  ia_mesh.removeElement(0);
  ia_mesh.removeVertex(0);

}


TEST(gtest_slam_IA, basic_tri_mesh)
{
  SLIC_INFO("Testing constructing basic triangle mesh...");
  const int vert_per_elem = 3;
  const int coord_per_vert = 3;
  typedef axom::slam::
    IAMesh<vert_per_elem-1, coord_per_vert,PointType>   IAMeshType;
  typedef IAMeshType::IndexListType IndexListType;

  BasicTriMeshData basic_mesh_data;
  IAMeshType ia_mesh(basic_mesh_data.points, basic_mesh_data.elem);

  EXPECT_TRUE( ia_mesh.isValid(true) );

  EXPECT_EQ(
    (int)basic_mesh_data.elem.size() / vert_per_elem,
    ia_mesh.getNumberOfElements() );

  EXPECT_EQ(
    (int)basic_mesh_data.points.size() / coord_per_vert,
    ia_mesh.getNumberOfVertices() );

  EXPECT_EQ( coord_per_vert, (int)ia_mesh.COORDS_PER_VERT);
  EXPECT_EQ( vert_per_elem, (int)ia_mesh.VERTS_PER_ELEM);

  // test mesh.getVerticesInElement( element_idx )
  for(int el_i = 0 ; el_i < ia_mesh.getNumberOfElements() ; el_i++ )
  {
    IndexListType r = ia_mesh.getVerticesInElement(el_i);
    EXPECT_EQ((int)r.size(), vert_per_elem);
    for(int i = 0 ; i<(int)r.size() ; i++)
    {
      EXPECT_EQ(r[i], basic_mesh_data.elem[el_i*vert_per_elem+i]);
    }
  }

  // test mesh.getElementsWithVertex( element_idx )
  for(int vert_i = 0 ; vert_i < ia_mesh.getNumberOfVertices() ; vert_i++)
  {
    IndexListType r = ia_mesh.getElementsWithVertex(vert_i);
    EXPECT_EQ( (int)r.size(), basic_mesh_data.vert_to_el_num[vert_i]);
    for(int i = 0 ; i<(int)r.size() ; i++)
    {
      int el = r[i];
      bool bContains = false;
      for(int j=0 ; j<vert_per_elem ; j++)
      {
        if(basic_mesh_data.elem[el*vert_per_elem+j] == vert_i)
        {
          bContains = true;
          break;
        }
      }
      EXPECT_TRUE(bContains);
    }
  }
  EXPECT_TRUE(ia_mesh.isManifold(true));

  SLIC_INFO("Done");
}


TEST(gtest_slam_IA, dynamically_build_tri_mesh)
{
  SLIC_INFO("Testing dynamically modifying a triangle mesh...");

  BasicTriMeshData basic_mesh_data;
  const int vert_per_elem = 3;
  const int coord_per_vert = 3;

  typedef axom::slam::IAMesh<vert_per_elem-1, coord_per_vert,
                             PointType>   IAMeshType;
  IAMeshType ia_mesh; //empty mesh

  EXPECT_TRUE( ia_mesh.isValid(true) );

  //build the mesh from nothing

  //adding the vertices
  for(int vert_i = 0 ; vert_i < (int)basic_mesh_data.points.size() ;
      vert_i+=vert_per_elem)
  {
    PointType pt( &basic_mesh_data.points[vert_i]);
    ia_mesh.addVertex( pt );
    EXPECT_TRUE( ia_mesh.isValid(true) );
  }

  //adding the elements
  for(int elem_i = 0 ; elem_i < (int)basic_mesh_data.elem.size() ;
      elem_i+=vert_per_elem)
  {
    ia_mesh.addElement(
      basic_mesh_data.elem[elem_i],
      basic_mesh_data.elem[elem_i + 1],
      basic_mesh_data.elem[elem_i + 2]);

    EXPECT_TRUE( ia_mesh.isValid(true) );
  }

  //check the ev_rel entries are correct
  for(int i=0 ; i< ia_mesh.ev_rel.size() ; i++)
  {
    EXPECT_EQ(ia_mesh.ev_rel[i].size(), vert_per_elem);
    for(int j=0 ; j<ia_mesh.ev_rel[i].size() ; j++)
    {
      EXPECT_EQ( ia_mesh.ev_rel[i][j], basic_mesh_data.elem[i*vert_per_elem+j]);
    }
  }

  //check the ve_rel entries are correct
  for(int v_i=0 ; v_i< ia_mesh.ve_rel.size() ; v_i++)
  {
    EXPECT_EQ(ia_mesh.ve_rel[v_i].size(), 1);
    int e_i = ia_mesh.ve_rel[v_i][0];
    bool bContains = false;
    for(int j=0 ; j<vert_per_elem ; j++)
    {
      bContains |= basic_mesh_data.elem[e_i*vert_per_elem+j];
    }
    EXPECT_TRUE(bContains);
  }

  //check the ee_rel entries are correct
  for(int i=0 ; i< ia_mesh.ee_rel.size() ; i++)
  {
    EXPECT_EQ(ia_mesh.ee_rel[i].size(), vert_per_elem);
    for(int j=0 ; j<ia_mesh.ee_rel[i].size() ; j++)
    {
      EXPECT_EQ( ia_mesh.ee_rel[i][j],
                 basic_mesh_data.el_nbr_rel[i*vert_per_elem+j] );
    }
  }

  //removing elements bigger than 7
  ia_mesh.removeElement(8);
  ia_mesh.removeElement(9);
  ia_mesh.removeElement(10);
  ia_mesh.removeElement(11);

  //removing vertex 1, which the removed elements contain
  ia_mesh.removeVertex(1);

  ia_mesh.isValid(true);

  EXPECT_EQ( ia_mesh.getNumberOfElements(),
             (int)basic_mesh_data.elem.size() / vert_per_elem - 4 );

  EXPECT_EQ( ia_mesh.getNumberOfVertices(),
             (int)basic_mesh_data.points.size() / coord_per_vert - 1 );

  //check the validity of the set entries
  for(int i=0 ; i< ia_mesh.getNumberOfElements() ; i++)
  {
    EXPECT_EQ( ia_mesh.isValidElementEntry(i), i<8 );
  }

  //check the ee_rel entries are correct after removing the elements
  for(int el_i=0 ; el_i< ia_mesh.ee_rel.size() ; el_i++)
  {
    EXPECT_EQ(ia_mesh.ee_rel[el_i].size(), vert_per_elem);
    for(int j=0 ; j<ia_mesh.ee_rel[el_i].size() ; j++)
    {
      int orig_nbr = basic_mesh_data.el_nbr_rel[el_i*vert_per_elem+j];
      if( orig_nbr > 7 || el_i > 7 )
      {
        EXPECT_EQ( ia_mesh.ee_rel[el_i][j],
                   (int)IAMeshType::ElementToElementRelation::INVALID_INDEX );
      }
      else
      {
        EXPECT_EQ( ia_mesh.ee_rel[el_i][j], orig_nbr );
      }
    }
  }

  //removing vertex 0
  ia_mesh.removeVertex(0);

  ia_mesh.isValid(true);

  //check that the elements with the removed vertices 0 are also removed.
  for(int i = 0 ; i < (int) basic_mesh_data.elem.size()/vert_per_elem ; i++)
  {
    bool bDeleted = false;
    for(int j=0 ; j<vert_per_elem ; j++)
    {
      bDeleted |= basic_mesh_data.elem[i*vert_per_elem+j] == 0;
    }
    EXPECT_EQ( ia_mesh.isValidElementEntry(i), !bDeleted && i<=7 );
  }
}


TEST(gtest_slam_IA, basic_tet_mesh)
{
  SLIC_INFO("Testing constructing basic tetrahedral mesh...");

  const int vert_per_elem = 4;
  const int coord_per_vert = 3;
  typedef axom::slam::
    IAMesh<vert_per_elem-1, coord_per_vert,PointType>   IAMeshType;
  typedef IAMeshType::IndexListType IndexListType;

  BasicTetMeshData basic_mesh_data;
  IAMeshType ia_mesh(basic_mesh_data.points, basic_mesh_data.elem);

  EXPECT_EQ((int)ia_mesh.VERTS_PER_ELEM, vert_per_elem);
  EXPECT_EQ((int)ia_mesh.COORDS_PER_VERT, coord_per_vert);

  EXPECT_TRUE( ia_mesh.isValid(true) );

  EXPECT_EQ(
    (int)basic_mesh_data.elem.size() / vert_per_elem,
    ia_mesh.getNumberOfElements());
  EXPECT_EQ(
    (int)basic_mesh_data.points.size() / coord_per_vert,
    ia_mesh.getNumberOfVertices());

  for(int el_i = 0 ; el_i < ia_mesh.getNumberOfElements() ; el_i++ )
  {
    IndexListType r = ia_mesh.getVerticesInElement(el_i);
    EXPECT_EQ((int)r.size(), vert_per_elem);
    for(int i = 0 ; i<(int)r.size() ; i++)
    {
      EXPECT_EQ(r[i], basic_mesh_data.elem[el_i*vert_per_elem+i]);
    }
  }

  for(int vert_i = 0 ; vert_i < ia_mesh.getNumberOfVertices() ; vert_i++)
  {
    IndexListType r = ia_mesh.getElementsWithVertex(vert_i);
    EXPECT_EQ( (int)r.size(), basic_mesh_data.vert_to_el_num[vert_i]);
    for(int i = 0 ; i<(int)r.size() ; i++)
    {
      int el = r[i];
      bool bContains = false;
      for(int j=0 ; j<vert_per_elem ; j++)
      {
        if(basic_mesh_data.elem[el*vert_per_elem+j] == vert_i)
        {
          bContains = true;
          break;
        }
      }
      EXPECT_TRUE(bContains);
    }
  }

  ia_mesh.print_all();

  EXPECT_TRUE(ia_mesh.isManifold(true));

  SLIC_INFO("Done");
}




TEST(gtest_slam_IA, dynamically_build_tet_mesh)
{
  SLIC_INFO("Testing dynamically modifying a tetrahedral mesh...");

  BasicTetMeshData basic_mesh_data;
  const int vert_per_elem = 4;
  const int coord_per_vert = 3;

  typedef axom::slam::
    IAMesh<vert_per_elem-1, coord_per_vert,PointType>   IAMeshType;
  IAMeshType ia_mesh; //empty mesh

  EXPECT_TRUE( ia_mesh.isValid(true) );

  //build the mesh from nothing

  //adding the vertices
  for(int vert_i = 0 ; vert_i < (int)basic_mesh_data.points.size() ;
      vert_i+=coord_per_vert)
  {
    PointType pt( &basic_mesh_data.points[vert_i]);
    ia_mesh.addVertex( pt );
    EXPECT_TRUE( ia_mesh.isValid(true) );
  }

  //adding the elements
  for(int elem_i = 0 ; elem_i < (int)basic_mesh_data.elem.size() ;
      elem_i+=vert_per_elem)
  {
    ia_mesh.addElement(
      basic_mesh_data.elem[elem_i],
      basic_mesh_data.elem[elem_i + 1],
      basic_mesh_data.elem[elem_i + 2],
      basic_mesh_data.elem[elem_i + 3]);

    EXPECT_TRUE( ia_mesh.isValid(true) );
  }

  //check the ev_rel entries are correct
  for(int i=0 ; i< ia_mesh.ev_rel.size() ; i++)
  {
    EXPECT_EQ(ia_mesh.ev_rel[i].size(), vert_per_elem);
    for(int j=0 ; j<ia_mesh.ev_rel[i].size() ; j++)
    {
      EXPECT_EQ( ia_mesh.ev_rel[i][j], basic_mesh_data.elem[i*vert_per_elem+j]);
    }
  }

  //check the ve_rel entries are correct
  for(int v_i=0 ; v_i< ia_mesh.ve_rel.size() ; v_i++)
  {
    EXPECT_EQ(ia_mesh.ve_rel[v_i].size(), 1);
    int e_i = ia_mesh.ve_rel[v_i][0];
    bool bContains = false;
    for(int j=0 ; j<vert_per_elem ; j++)
    {
      bContains |= basic_mesh_data.elem[e_i*vert_per_elem+j];
    }
    EXPECT_TRUE(bContains);
  }

  //check the ee_rel entries are correct
  for(int i=0 ; i< ia_mesh.ee_rel.size() ; i++)
  {
    EXPECT_EQ(ia_mesh.ee_rel[i].size(), vert_per_elem);
    for(int j=0 ; j<ia_mesh.ee_rel[i].size() ; j++)
    {
      EXPECT_EQ( ia_mesh.ee_rel[i][j],
                 basic_mesh_data.el_nbr_rel[i*vert_per_elem+j] );
    }
  }

  //removing elements bigger than 3
  ia_mesh.removeElement(4);
  ia_mesh.removeElement(5);

  //removing vertex 1, which only the removed elements contain
  ia_mesh.removeVertex(5);

  ia_mesh.isValid(true);

  EXPECT_EQ( ia_mesh.getNumberOfElements(),
             (int)basic_mesh_data.elem.size() / vert_per_elem - 2 );

  EXPECT_EQ( ia_mesh.getNumberOfVertices(),
             (int)basic_mesh_data.points.size() / coord_per_vert - 1 );

  //check the validity of the set entries
  for(int i=0 ; i< ia_mesh.getNumberOfElements() ; i++)
  {
    EXPECT_EQ( ia_mesh.isValidElementEntry(i), i < 4 );
  }

  //check the ee_rel entries are correct after removing the elements
  for(int el_i=0 ; el_i< ia_mesh.ee_rel.size() ; el_i++)
  {
    EXPECT_EQ(ia_mesh.ee_rel[el_i].size(), vert_per_elem);
    for(int j=0 ; j<ia_mesh.ee_rel[el_i].size() ; j++)
    {
      int orig_nbr = basic_mesh_data.el_nbr_rel[el_i*vert_per_elem+j];
      if( orig_nbr < 4 && el_i < 4 )
      {
        EXPECT_EQ( ia_mesh.ee_rel[el_i][j], orig_nbr );
      }
      else
      {
        EXPECT_EQ( ia_mesh.ee_rel[el_i][j],
                   (int)IAMeshType::ElementToElementRelation::INVALID_INDEX );
      }
    }
  }

  //removing vertex 6
  ia_mesh.removeVertex(6);

  //check that the elements with the removed vertices are also removed.
  for(int i = 0 ; i < (int) basic_mesh_data.elem.size()/vert_per_elem ; i++)
  {
    bool bDeleted = false;
    for(int j=0 ; j<vert_per_elem ; j++)
    {
      bDeleted |= basic_mesh_data.elem[i*vert_per_elem+j] == 6;
    }
    EXPECT_EQ( ia_mesh.isValidElementEntry(i), !bDeleted && i<4 );
  }

  SLIC_INFO("Done");
}



//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  // create & initialize test logger. finalized when exiting main scope
  UnitTestLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
