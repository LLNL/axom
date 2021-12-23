// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file slam_IA.cpp
 *
 * \brief This file contains tests for slam's IA mesh
 */

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/slam/mesh_struct/IA.hpp"
#include "axom/slam/Utilities.hpp"

namespace slam = axom::slam;

namespace
{
using PointType = slam::util::Point3<double>;
using IndexType = axom::IndexType;

struct BasicTriMeshData
{
  //This is the surface of a regular cube composed of 12 triangles
  std::vector<double> points;
  std::vector<IndexType> elem;

  static const double point_arr[];
  static const IndexType tri_arr[];
  static const IndexType vert_to_el_num[];  // for verification purpose
  static const IndexType el_nbr_rel[];

  BasicTriMeshData()
    : points(point_arr, point_arr + 24)
    , elem(tri_arr, tri_arr + 36)
  { }

  static int numVertices() { return 8; }
  static int numTriangles() { return 12; }
};

// clang-format off
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
const IndexType BasicTriMeshData::tri_arr [] = {
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
const IndexType BasicTriMeshData::vert_to_el_num [] = {
  6, 4, 4, 4, 4, 4, 4, 6
};
const IndexType BasicTriMeshData::el_nbr_rel [] = {
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
// clang-format on

struct BasicTetMeshData
{
  //cube divided into six tets

  std::vector<double> points;
  std::vector<IndexType> elem;

  static const double point_arr[];
  static const IndexType tet_arr[];
  static const IndexType vert_to_el_num[];
  static const IndexType el_nbr_rel[];

  BasicTetMeshData()
    : points(point_arr, point_arr + 24)
    , elem(tet_arr, tet_arr + 24)
  { }

  static int numVertices() { return 8; }
  static int numTetrahedra() { return 6; }
};

// clang-format off
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
const IndexType BasicTetMeshData::tet_arr [] = {
  3,2,4,0,
  3,1,4,0,
  3,6,2,4,
  3,6,7,4,
  3,5,1,4,
  3,5,7,4
};
const IndexType BasicTetMeshData::vert_to_el_num [] = {
  2, 2, 2, 6, 6, 2, 2, 2
};
const IndexType BasicTetMeshData::el_nbr_rel [] = {
  2,-1,1,-1,
  4,-1,0,-1,
  -1,-1,0,3,
  -1,-1,5,2,
  -1,-1,1,5,
  -1,-1,3,4
};
// clang-format on

/* Check if vertex with ID \a vert_id is in boundary of
 * element with ID \a elem_id in mesh \a ia_mesh
 */
template <typename IAMeshType>
bool isInBoundary(const IAMeshType& ia_mesh, IndexType vert_id, IndexType elem_id)
{
  typename IAMeshType::IndexArray ev = ia_mesh.getVerticesInElement(elem_id);
  const int sz = ev.size();
  for(int j = 0; j < sz; ++j)
  {
    if(ev[j] == vert_id)
    {
      return true;
    }
  }

  return false;
}

/* Check if element with ID \a el_1 is adjacent to
 * element with ID \a el_2 in mesh \a ia_mesh
 */
template <typename IAMeshType>
bool isAdjacent(const IAMeshType& ia_mesh, IndexType el_1, IndexType el_2)
{
  typename IAMeshType::IndexArray ee = ia_mesh.getElementNeighbors(el_2);
  const int sz = ee.size();
  for(int j = 0; j < sz; ++j)
  {
    if(ee[j] == el_1)
    {
      return true;
    }
  }

  return false;
}

}  // end anonymous namespace

TEST(slam_IA, empty_mesh)
{
  SLIC_INFO("Testing creating an empty IA mesh...");

  const int TDIM = 2;
  const int SDIM = 3;
  using IAMeshType = slam::IAMesh<TDIM, SDIM, PointType>;

  IAMeshType ia_mesh;  //empty mesh

  EXPECT_TRUE(ia_mesh.isValid(true));
  EXPECT_TRUE(ia_mesh.isEmpty());
  EXPECT_EQ(0, ia_mesh.getNumberOfElements());
  EXPECT_EQ(0, ia_mesh.getNumberOfVertices());

  // Check boundary relation on non-existent element 0 (w/ warning)
  {
    int sz = ia_mesh.getVerticesInElement(0).size();
    EXPECT_EQ(0, sz);
  }

  // Check coboundary relation on non-existent vertex 0 (w/ warning)
  {
    int sz = ia_mesh.getElementsWithVertex(0).size();
    EXPECT_EQ(0, sz);
  }

  // Check adjacency relation on non-existent element 0 (w/ warning)
  {
    int sz = ia_mesh.getElementNeighbors(0).size();
    EXPECT_EQ(0, sz);
  }

  // Removing invalid vertices/elements is a no-op (w/ warning)
  ia_mesh.removeElement(0);
  ia_mesh.removeVertex(0);
}

TEST(slam_IA, basic_tri_mesh)
{
  SLIC_INFO("Testing constructing basic triangle mesh...");

  const int TDIM = 2;
  const int SDIM = 3;
  using IAMeshType = slam::IAMesh<TDIM, SDIM, PointType>;

  using IndexArray = IAMeshType::IndexArray;

  const int vert_per_elem = IAMeshType::VERTS_PER_ELEM;
  EXPECT_EQ(3, vert_per_elem);

  const int coord_per_vert = IAMeshType::COORDS_PER_VERT;
  EXPECT_EQ(3, coord_per_vert);

  const int numAdjacentElems = vert_per_elem;

  BasicTriMeshData basic_mesh_data;
  IAMeshType ia_mesh(basic_mesh_data.points, basic_mesh_data.elem);

  EXPECT_TRUE(ia_mesh.isValid(true));
  EXPECT_FALSE(ia_mesh.isEmpty());

  const int numVerts = ia_mesh.getNumberOfVertices();
  EXPECT_EQ(basic_mesh_data.numVertices(), numVerts);

  const int numElems = ia_mesh.getNumberOfElements();
  EXPECT_EQ(basic_mesh_data.numTriangles(), numElems);

  // Test Element boundary relation
  for(int e_i = 0; e_i < numElems; ++e_i)
  {
    IndexArray r = ia_mesh.getVerticesInElement(e_i);
    const int numEV = r.size();
    EXPECT_EQ(vert_per_elem, numEV);
    for(int i = 0; i < numEV; i++)
    {
      EXPECT_EQ(basic_mesh_data.elem[e_i * vert_per_elem + i], r[i]);
    }
  }

  // Test Vertex co-boundary relation
  for(int v_i = 0; v_i < numVerts; ++v_i)
  {
    IndexArray r = ia_mesh.getElementsWithVertex(v_i);
    const int numVE = r.size();
    EXPECT_EQ(basic_mesh_data.vert_to_el_num[v_i], numVE);
    for(int i = 0; i < numVE; i++)
    {
      int el = r[i];
      EXPECT_TRUE(isInBoundary(ia_mesh, v_i, el))
        << "Vertex co-boundary relation indicates that " << v_i
        << " should be in boundary of element " << el << " but it is not";
    }
  }

  // Test Element adjacency relation
  for(int e_i = 0; e_i < numElems; ++e_i)
  {
    IndexArray r = ia_mesh.getElementNeighbors(e_i);
    const int numEE = r.size();
    EXPECT_EQ(numAdjacentElems, numEE);
    for(int i = 0; i < numEE; ++i)
    {
      int e_adj = r[i];
      EXPECT_TRUE(isAdjacent(ia_mesh, e_adj, e_i))
        << "Element adjacency relation indicates that " << e_adj
        << " should be adjacent to element " << e_i << " but it is not";
    }
  }

  EXPECT_TRUE(ia_mesh.isManifold(true));

  ia_mesh.print_all();

  SLIC_INFO("Done");
}

TEST(slam_IA, dynamically_build_tri_mesh)
{
  SLIC_INFO("Testing dynamically modifying a triangle mesh...");

  BasicTriMeshData basic_mesh_data;
  const int TDIM = 2;
  const int SDIM = 3;
  using IAMeshType = slam::IAMesh<TDIM, SDIM, PointType>;
  using IndexArray = IAMeshType::IndexArray;

  const int vert_per_elem = IAMeshType::VERTS_PER_ELEM;
  EXPECT_EQ(3, vert_per_elem);

  const int coord_per_vert = IAMeshType::COORDS_PER_VERT;
  EXPECT_EQ(3, coord_per_vert);

  const int numAdjacentElems = vert_per_elem;

  // Build the mesh from nothing
  IAMeshType ia_mesh;  //empty mesh
  EXPECT_TRUE(ia_mesh.isValid(true));
  EXPECT_TRUE(ia_mesh.isEmpty());

  // Add the vertices
  for(int v_i = 0; v_i < basic_mesh_data.numVertices(); ++v_i)
  {
    EXPECT_EQ(v_i, ia_mesh.getNumberOfVertices());

    PointType pt(&basic_mesh_data.points[v_i * coord_per_vert]);
    ia_mesh.addVertex(pt);
    EXPECT_TRUE(ia_mesh.isValid(true));
  }

  const int numVerts = ia_mesh.getNumberOfVertices();
  EXPECT_EQ(basic_mesh_data.numVertices(), numVerts);
  EXPECT_FALSE(ia_mesh.isEmpty());

  // Add the elements
  for(int e_i = 0; e_i < basic_mesh_data.numTriangles(); ++e_i)
  {
    EXPECT_EQ(e_i, ia_mesh.getNumberOfElements());

    int startIdx = e_i * vert_per_elem;
    ia_mesh.addElement(basic_mesh_data.elem[startIdx + 0],
                       basic_mesh_data.elem[startIdx + 1],
                       basic_mesh_data.elem[startIdx + 2]);

    EXPECT_TRUE(ia_mesh.isValid(true));
  }
  const int numTris = ia_mesh.getNumberOfElements();
  EXPECT_EQ(basic_mesh_data.numTriangles(), numTris);
  EXPECT_FALSE(ia_mesh.isEmpty());

  // Check that the element boundary relation is correct
  for(int e_i = 0; e_i < numTris; ++e_i)
  {
    IndexArray r = ia_mesh.getVerticesInElement(e_i);
    const int sz = r.size();
    EXPECT_EQ(vert_per_elem, sz);

    for(int j = 0; j < sz; ++j)
    {
      const int expVert = basic_mesh_data.elem[e_i * vert_per_elem + j];
      EXPECT_EQ(expVert, r[j]);
    }
  }

  // Check that the vertex co-boundary relation is correct
  for(int v_i = 0; v_i < numVerts; ++v_i)
  {
    IndexArray r = ia_mesh.getElementsWithVertex(v_i);
    const int numVE = r.size();
    EXPECT_EQ(basic_mesh_data.vert_to_el_num[v_i], numVE);

    for(int i = 0; i < numVE; i++)
    {
      int el = r[i];
      EXPECT_TRUE(isInBoundary(ia_mesh, v_i, el))
        << "Vertex co-boundary relation indicates that " << v_i
        << " should be in boundary of element " << el << " but it is not";
    }
  }

  // Check the element adjacencies are correct
  for(int e_i = 0; e_i < numTris; ++e_i)
  {
    IndexArray r = ia_mesh.getElementNeighbors(e_i);
    const int numEE = r.size();
    EXPECT_EQ(numAdjacentElems, numEE);
    for(int i = 0; i < numEE; ++i)
    {
      int el_adj = r[i];
      EXPECT_TRUE(isAdjacent(ia_mesh, el_adj, e_i))
        << "Element adjacency relation indicates that " << el_adj
        << " should be adjacent to element " << e_i << " but it is not";
    }
  }
}

TEST(slam_IA, tri_mesh_remove_verts_and_elems)
{
  SLIC_INFO("Testing removing elements and vertices from a triangle mesh...");

  const int TDIM = 2;
  const int SDIM = 3;
  using IAMeshType = slam::IAMesh<TDIM, SDIM, PointType>;

  const int vert_per_elem = IAMeshType::VERTS_PER_ELEM;

  BasicTriMeshData basic_mesh_data;
  IAMeshType ia_mesh(basic_mesh_data.points, basic_mesh_data.elem);

  // Check that we begin with the correct number of verts and tris
  EXPECT_TRUE(ia_mesh.isValid(true));
  EXPECT_EQ(basic_mesh_data.numTriangles(), ia_mesh.getNumberOfElements());

  EXPECT_EQ(basic_mesh_data.numVertices(), ia_mesh.getNumberOfVertices());

  //removing elements bigger than 7
  ia_mesh.removeElement(8);
  ia_mesh.removeElement(9);
  ia_mesh.removeElement(10);
  ia_mesh.removeElement(11);

  EXPECT_FALSE(ia_mesh.isEmpty());

  EXPECT_TRUE(ia_mesh.isValid(true));
  EXPECT_EQ(basic_mesh_data.numTriangles() - 4, ia_mesh.getNumberOfElements());

  //removing vertex 1, which the removed elements contain
  ia_mesh.removeVertex(1);

  EXPECT_TRUE(ia_mesh.isValid(true));
  EXPECT_FALSE(ia_mesh.isEmpty());

  EXPECT_EQ(basic_mesh_data.numTriangles() - 4, ia_mesh.getNumberOfElements());

  EXPECT_EQ(basic_mesh_data.numVertices() - 1, ia_mesh.getNumberOfVertices());

  // Check the validity of the set entries
  for(int v_i = 0; v_i < ia_mesh.getNumberOfVertices(); ++v_i)
  {
    EXPECT_EQ(ia_mesh.isValidVertexEntry(v_i), v_i != 1);
  }

  for(int e_i = 0; e_i < ia_mesh.getNumberOfElements(); ++e_i)
  {
    EXPECT_EQ(ia_mesh.isValidElementEntry(e_i), e_i < 8);
  }

  // Check incidence and adjacency relations
  // TODO: This should not use private relation data
  for(int el_i = 0; el_i < ia_mesh.ee_rel.size(); el_i++)
  {
    EXPECT_EQ(ia_mesh.ee_rel[el_i].size(), vert_per_elem);
    for(int j = 0; j < ia_mesh.ee_rel[el_i].size(); j++)
    {
      int orig_nbr = basic_mesh_data.el_nbr_rel[el_i * vert_per_elem + j];
      if(orig_nbr > 7 || el_i > 7)
      {
        EXPECT_EQ(ia_mesh.ee_rel[el_i][j],
                  (int)IAMeshType::ElementAdjacencyRelation::INVALID_INDEX);
      }
      else
      {
        EXPECT_EQ(ia_mesh.ee_rel[el_i][j], orig_nbr);
      }
    }
  }

  //removing vertex 0
  ia_mesh.removeVertex(0);

  ia_mesh.isValid(true);

  //check that the elements with the removed vertices 0 are also removed.
  for(int i = 0; i < (int)basic_mesh_data.elem.size() / vert_per_elem; i++)
  {
    bool bDeleted = false;
    for(int j = 0; j < vert_per_elem; j++)
    {
      bDeleted |= basic_mesh_data.elem[i * vert_per_elem + j] == 0;
    }
    EXPECT_EQ(ia_mesh.isValidElementEntry(i), !bDeleted && i <= 7);
  }
}

TEST(slam_IA, tri_mesh_remove_elem_and_compact)
{
  SLIC_INFO("Testing removing an element and compacting a triangle mesh...");

  const int TDIM = 2;
  const int SDIM = 3;
  using IAMeshType = slam::IAMesh<TDIM, SDIM, PointType>;

  BasicTriMeshData basic_mesh_data;
  IAMeshType ia_mesh(basic_mesh_data.points, basic_mesh_data.elem);

  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_EQ(basic_mesh_data.numTriangles(), ia_mesh.getNumberOfElements());
  EXPECT_EQ(basic_mesh_data.numVertices(), ia_mesh.getNumberOfVertices());

  ia_mesh.removeElement(4);
  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_EQ(basic_mesh_data.numTriangles() - 1, ia_mesh.getNumberOfElements());
  EXPECT_EQ(basic_mesh_data.numVertices(), ia_mesh.getNumberOfVertices());

  ia_mesh.compact();
  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_EQ(basic_mesh_data.numTriangles() - 1, ia_mesh.getNumberOfElements());
  EXPECT_EQ(basic_mesh_data.numVertices(), ia_mesh.getNumberOfVertices());
}

TEST(slam_IA, tri_mesh_remove_vert_and_compact)
{
  SLIC_INFO("Testing removing a vertex and compacting a triangle mesh...");

  const int TDIM = 2;
  const int SDIM = 3;
  using IAMeshType = slam::IAMesh<TDIM, SDIM, PointType>;

  BasicTriMeshData basic_mesh_data;
  IAMeshType ia_mesh(basic_mesh_data.points, basic_mesh_data.elem);

  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_EQ(basic_mesh_data.numTriangles(), ia_mesh.getNumberOfElements());
  EXPECT_EQ(basic_mesh_data.numVertices(), ia_mesh.getNumberOfVertices());

  ia_mesh.removeVertex(3);
  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_EQ(basic_mesh_data.numTriangles() - 4, ia_mesh.getNumberOfElements());
  EXPECT_EQ(basic_mesh_data.numVertices() - 1, ia_mesh.getNumberOfVertices());

  ia_mesh.compact();
  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_EQ(basic_mesh_data.numTriangles() - 4, ia_mesh.getNumberOfElements());
  EXPECT_EQ(basic_mesh_data.numVertices() - 1, ia_mesh.getNumberOfVertices());
}

TEST(slam_IA, basic_tet_mesh)
{
  SLIC_INFO("Testing constructing basic tetrahedral mesh...");

  const int vert_per_elem = 4;
  const int coord_per_vert = 3;
  using IAMeshType = slam::IAMesh<vert_per_elem - 1, coord_per_vert, PointType>;
  using IndexArray = IAMeshType::IndexArray;

  BasicTetMeshData basic_mesh_data;
  IAMeshType ia_mesh(basic_mesh_data.points, basic_mesh_data.elem);

  EXPECT_EQ((int)ia_mesh.VERTS_PER_ELEM, vert_per_elem);
  EXPECT_EQ((int)ia_mesh.COORDS_PER_VERT, coord_per_vert);

  EXPECT_TRUE(ia_mesh.isValid(true));

  EXPECT_EQ((int)basic_mesh_data.elem.size() / vert_per_elem,
            ia_mesh.getNumberOfElements());
  EXPECT_EQ((int)basic_mesh_data.points.size() / coord_per_vert,
            ia_mesh.getNumberOfVertices());

  for(int el_i = 0; el_i < ia_mesh.getNumberOfElements(); el_i++)
  {
    IndexArray r = ia_mesh.getVerticesInElement(el_i);
    EXPECT_EQ((int)r.size(), vert_per_elem);
    for(int i = 0; i < (int)r.size(); i++)
    {
      EXPECT_EQ(r[i], basic_mesh_data.elem[el_i * vert_per_elem + i]);
    }
  }

  for(int vert_i = 0; vert_i < ia_mesh.getNumberOfVertices(); vert_i++)
  {
    IndexArray r = ia_mesh.getElementsWithVertex(vert_i);
    EXPECT_EQ((int)r.size(), basic_mesh_data.vert_to_el_num[vert_i]);
    for(int i = 0; i < (int)r.size(); i++)
    {
      int el = r[i];
      bool bContains = false;
      for(int j = 0; j < vert_per_elem; j++)
      {
        if(basic_mesh_data.elem[el * vert_per_elem + j] == vert_i)
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

TEST(slam_IA, dynamically_build_tet_mesh)
{
  SLIC_INFO("Testing dynamically modifying a tetrahedral mesh...");

  BasicTetMeshData basic_mesh_data;
  const int vert_per_elem = 4;
  const int coord_per_vert = 3;

  using IAMeshType = slam::IAMesh<vert_per_elem - 1, coord_per_vert, PointType>;
  IAMeshType ia_mesh;  //empty mesh

  EXPECT_TRUE(ia_mesh.isValid(true));

  //build the mesh from nothing

  //adding the vertices
  for(int vert_i = 0; vert_i < (int)basic_mesh_data.points.size();
      vert_i += coord_per_vert)
  {
    PointType pt(&basic_mesh_data.points[vert_i]);
    ia_mesh.addVertex(pt);
    EXPECT_TRUE(ia_mesh.isValid(true));
  }

  //adding the elements
  for(int elem_i = 0; elem_i < (int)basic_mesh_data.elem.size();
      elem_i += vert_per_elem)
  {
    ia_mesh.addElement(basic_mesh_data.elem[elem_i],
                       basic_mesh_data.elem[elem_i + 1],
                       basic_mesh_data.elem[elem_i + 2],
                       basic_mesh_data.elem[elem_i + 3]);

    EXPECT_TRUE(ia_mesh.isValid(true));
  }

  //check the ev_rel entries are correct
  for(int i = 0; i < ia_mesh.ev_rel.size(); i++)
  {
    EXPECT_EQ(ia_mesh.ev_rel[i].size(), vert_per_elem);
    for(int j = 0; j < ia_mesh.ev_rel[i].size(); j++)
    {
      EXPECT_EQ(ia_mesh.ev_rel[i][j],
                basic_mesh_data.elem[i * vert_per_elem + j]);
    }
  }

  //check the ve_rel entries are correct
  for(int v_i = 0; v_i < ia_mesh.ve_rel.size(); v_i++)
  {
    EXPECT_EQ(ia_mesh.ve_rel[v_i].size(), 1);
    int e_i = ia_mesh.ve_rel[v_i][0];
    bool bContains = false;
    for(int j = 0; j < vert_per_elem; j++)
    {
      bContains |= basic_mesh_data.elem[e_i * vert_per_elem + j];
    }
    EXPECT_TRUE(bContains);
  }

  //check the ee_rel entries are correct
  for(int i = 0; i < ia_mesh.ee_rel.size(); i++)
  {
    EXPECT_EQ(ia_mesh.ee_rel[i].size(), vert_per_elem);
    for(int j = 0; j < ia_mesh.ee_rel[i].size(); j++)
    {
      EXPECT_EQ(ia_mesh.ee_rel[i][j],
                basic_mesh_data.el_nbr_rel[i * vert_per_elem + j]);
    }
  }

  //removing elements bigger than 3
  ia_mesh.removeElement(4);
  ia_mesh.removeElement(5);

  //removing vertex 1, which only the removed elements contain
  ia_mesh.removeVertex(5);

  ia_mesh.isValid(true);

  EXPECT_EQ(ia_mesh.getNumberOfElements(),
            (int)basic_mesh_data.elem.size() / vert_per_elem - 2);

  EXPECT_EQ(ia_mesh.getNumberOfVertices(),
            (int)basic_mesh_data.points.size() / coord_per_vert - 1);

  //check the validity of the set entries
  for(int i = 0; i < ia_mesh.getNumberOfElements(); i++)
  {
    EXPECT_EQ(ia_mesh.isValidElementEntry(i), i < 4);
  }

  //check the ee_rel entries are correct after removing the elements
  for(int el_i = 0; el_i < ia_mesh.ee_rel.size(); el_i++)
  {
    EXPECT_EQ(ia_mesh.ee_rel[el_i].size(), vert_per_elem);
    for(int j = 0; j < ia_mesh.ee_rel[el_i].size(); j++)
    {
      int orig_nbr = basic_mesh_data.el_nbr_rel[el_i * vert_per_elem + j];
      if(orig_nbr < 4 && el_i < 4)
      {
        EXPECT_EQ(ia_mesh.ee_rel[el_i][j], orig_nbr);
      }
      else
      {
        EXPECT_EQ(ia_mesh.ee_rel[el_i][j],
                  (int)IAMeshType::ElementAdjacencyRelation::INVALID_INDEX);
      }
    }
  }

  //removing vertex 6
  ia_mesh.removeVertex(6);

  //check that the elements with the removed vertices are also removed.
  for(int i = 0; i < (int)basic_mesh_data.elem.size() / vert_per_elem; i++)
  {
    bool bDeleted = false;
    for(int j = 0; j < vert_per_elem; j++)
    {
      bDeleted |= basic_mesh_data.elem[i * vert_per_elem + j] == 6;
    }
    EXPECT_EQ(ia_mesh.isValidElementEntry(i), !bDeleted && i < 4);
  }

  SLIC_INFO("Done");
}

TEST(slam_IA, compact_mesh)
{
  const int TDIM = 2;
  const int SDIM = 3;
  using IAMeshType = slam::IAMesh<TDIM, SDIM, PointType>;
  using IndexType = IAMeshType::IndexType;

  IAMeshType mesh;

  mesh.addVertex(PointType(0, 0, 1));
  mesh.addVertex(PointType(0, 1, 0));
  mesh.addVertex(PointType(1, 0, 0));
  mesh.addVertex(PointType(1, 1, 1));
  mesh.addElement(0, 1, 2);
  mesh.addElement(1, 2, 3);

  mesh.removeVertex(0);

  EXPECT_TRUE(mesh.isValid(true));

  auto v_before = mesh.getNumberOfVertices();
  auto e_before = mesh.getNumberOfElements();
  EXPECT_EQ(IndexType(3), v_before);
  EXPECT_EQ(IndexType(1), e_before);

  mesh.compact();

  auto v_after = mesh.getNumberOfVertices();
  auto e_after = mesh.getNumberOfElements();

  EXPECT_TRUE(mesh.isValid(true));

  EXPECT_EQ(v_before, v_after);
  EXPECT_EQ(e_before, e_after);
}

TEST(slam_IA, tet_mesh_remove_elem_and_compact)
{
  SLIC_INFO("Testing removing an element and compacting a tetrahedral mesh...");

  const int TDIM = 3;
  const int SDIM = 3;
  using IAMeshType = slam::IAMesh<TDIM, SDIM, PointType>;

  BasicTetMeshData basic_mesh_data;
  IAMeshType ia_mesh(basic_mesh_data.points, basic_mesh_data.elem);

  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_EQ(basic_mesh_data.numTetrahedra(), ia_mesh.getNumberOfElements());
  EXPECT_EQ(basic_mesh_data.numVertices(), ia_mesh.getNumberOfVertices());

  ia_mesh.removeElement(4);
  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_EQ(basic_mesh_data.numTetrahedra() - 1, ia_mesh.getNumberOfElements());
  EXPECT_EQ(basic_mesh_data.numVertices(), ia_mesh.getNumberOfVertices());

  ia_mesh.compact();
  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_EQ(basic_mesh_data.numTetrahedra() - 1, ia_mesh.getNumberOfElements());
  EXPECT_EQ(basic_mesh_data.numVertices(), ia_mesh.getNumberOfVertices());
}

TEST(slam_IA, tet_mesh_remove_vert_and_compact)
{
  SLIC_INFO("Testing removing a vertex and compacting a tetrahedral mesh...");

  const int TDIM = 3;
  const int SDIM = 3;
  using IAMeshType = slam::IAMesh<TDIM, SDIM, PointType>;

  BasicTetMeshData basic_mesh_data;
  IAMeshType ia_mesh(basic_mesh_data.points, basic_mesh_data.elem);

  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_EQ(basic_mesh_data.numTetrahedra(), ia_mesh.getNumberOfElements());
  EXPECT_EQ(basic_mesh_data.numVertices(), ia_mesh.getNumberOfVertices());

  ia_mesh.removeVertex(2);
  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_EQ(basic_mesh_data.numTetrahedra() - 2, ia_mesh.getNumberOfElements());
  EXPECT_EQ(basic_mesh_data.numVertices() - 1, ia_mesh.getNumberOfVertices());

  ia_mesh.compact();
  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_EQ(basic_mesh_data.numTetrahedra() - 2, ia_mesh.getNumberOfElements());
  EXPECT_EQ(basic_mesh_data.numVertices() - 1, ia_mesh.getNumberOfVertices());
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
