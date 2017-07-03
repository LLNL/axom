// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 * \file TopologicalMesh.cpp
 *
 * \brief Example of using SLAM to build topological data structure
 *
 */

#include "axom/core.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"

#include <list>
#include <set>
#include <stdlib.h>

namespace slamTopologicalMesh
{
using IndexType = axom::IndexType;
using DataType = double;

using Point2D = axom::primal::Point2D;
using Point3D = axom::primal::Point3D;
using Vector2D = axom::primal::Vector2D;
using Triangle2D = axom::primal::Triangle<double, 2>;

struct BasicTriMeshData
{
  //This is a regular cube
  std::vector<DataType> points = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0,
                                  0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0,
                                  0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0};
  std::vector<IndexType> tri = {0, 6, 4, 0, 2, 6, 0, 3, 2, 0, 1, 3,
                                2, 7, 6, 2, 3, 7, 4, 6, 7, 4, 7, 5,
                                0, 4, 5, 0, 5, 1, 1, 5, 7, 1, 7, 3};
};

struct Basic2DTriMeshData
{
  //This is a square
  std::vector<DataType> points = {0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0};
  std::vector<IndexType> tri = {0, 1, 2, 1, 3, 2};
};

struct BasicTetMeshData
{
  //cube divided into tets

  std::vector<DataType> points = {
    -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0, 1.0,
    1.0,  -1.0, -1.0, 1.0,  -1.0, 1.0, 1.0,  1.0, -1.0, 1.0,  1.0, 1.0};

  std::vector<IndexType> tet = {3, 2, 4, 0, 3, 1, 4, 0, 3, 6, 2, 4,
                                3, 6, 7, 4, 3, 5, 1, 4, 3, 5, 7, 4};
};

struct Flat2DTriMeshData
{
  //This is a square
  std::vector<DataType> points =
    {-1, -1, -1, 1, 1, -1, 1, 1, -2, -2, -2, 2, 2, -2, 2, 2};
  std::vector<IndexType> tri = {

    1, 3, 5, 4, 0, 5, 0, 1, 5, 0, 2, 1, 0, 4, 2,
    3, 7, 5, 3, 2, 6, 3, 1, 2, 3, 6, 7, 2, 4, 6};
};

//typedef axom::slam::Map< Point2D >            XYZField;
/*
  void calculate_normal(axom::slam::IAMesh<3> & tri_mesh){

    //calculate normal of each element as a map
    XYZField element_to_normal_map = XYZField(&tri_mesh.element_set);
    for(IndexType zidx = 0; zidx < tri_mesh.element_set.size(); zidx++)
    {

      Point p0 = tri_mesh.vcoord_map[ tri_mesh.ev_rel[zidx][0] ];
      Point p1 = tri_mesh.vcoord_map[ tri_mesh.ev_rel[zidx][1] ];
      Point p2 = tri_mesh.vcoord_map[ tri_mesh.ev_rel[zidx][2] ];
      Point p01 = p1 - p0;
      Point p02 = p2 - p0;
      //Do we assume the triangle is not degenerate?
      Point pnormal = cross( p01, p02);
      element_to_normal_map[zidx] = pnormal;
      //SLIC_INFO("normal "<< pnormal.m_x <<" " << pnormal.m_y << " " <<pnormal.m_z);
    }

    //calculate normal of each vertex as a map
    XYZField vertex_to_normal_map = XYZField(&tri_mesh.element_set);
    for(IndexType nidx = 0; nidx < tri_mesh.vertex_set.size(); nidx++)
    {
      Point vertex_normal(0,0,0);
      std::vector<IndexType> element_subset = tri_mesh.getElementWithVertex(nidx);
      for(unsigned int zidx = 0; zidx < element_subset.size(); zidx++)
      {
        Point& element_normal = element_to_normal_map[element_subset[zidx]];
        vertex_normal += element_normal;
      }
      vertex_normal = normalize(vertex_normal);
      vertex_to_normal_map[nidx] = vertex_normal;
      SLIC_INFO("normal "<< vertex_normal.m_x <<" " << vertex_normal.m_y << " " << vertex_normal.m_z);
    }

  }
*/

IndexType find_containing_element(const axom::slam::IAMesh<3, 2>& mesh,
                                  const Point2D& query_pt,
                                  IndexType element_i)
{
  SLIC_INFO("find_containing_element " << element_i);
  //assumes the point is definitely contained within a triangle
  //TODO fix that assumption ^
  SLIC_INFO("Query Pt " << query_pt);

  if(!mesh.isValidElementEntry(element_i))
  {
    if(mesh.getNumberOfElements() == 0)
    {  //There are no elements in the mesh? How did this happen?
      SLIC_WARNING("Warning: no element in the mesh.");
      return -1;
    }

    SLIC_WARNING("Warning: Element "
                 << element_i
                 << " is not a valid element. Replacing with a valid one.");

    //find the first valid element
    for(int i = 0; true; i++)
    {
      if(mesh.isValidElementEntry(i))
      {
        element_i = i;
        break;
      }
    }
  }

  while(1)
  {
    SLIC_INFO("step into element" << element_i);
    std::vector<IndexType> verts = mesh.getVerticesInElement(element_i);

    Triangle2D tri(mesh.getVertexPoint(verts[0]),
                   mesh.getVertexPoint(verts[1]),
                   mesh.getVertexPoint(verts[2]));

    axom::primal::Point3D bary_co = tri.barycentricCoords(query_pt);
    SLIC_INFO("bary_co " << bary_co);

    //Find the most negative
    IndexType i = 0;
    if(bary_co[0] > bary_co[1]) i = 1;
    if(bary_co[i] > bary_co[2]) i = 2;

    if(bary_co[i] >= 0)
    {                    //smaller than i -> outside of the triangle
      return element_i;  //return if inside or on triangle
    }

    std::vector<IndexType> zlist = mesh.getElementNeighbor(element_i);
    IndexType next_el = element_i = zlist[2 - i];

    if(next_el < 0)
    {
      //Either there is a hole in the mesh, or the point is outside of the mesh
      return element_i;
    }
  }
}

// this is moved to primal
/*
  bool is_point_in_circle(const Point2D & p0, const Point2D & p1, const Point2D & p2, const Point2D & q){
    double det = axom::numerics::determinant(
        1.0, p0[0], p0[1], p0[0]*p0[0] + p0[1]*p0[1],
        1.0, p1[0], p1[1], p1[0]*p1[0] + p1[1]*p1[1],
        1.0, p2[0], p2[1], p2[0]*p2[0] + p2[1]*p2[1],
        1.0, q[0],  q[1],  q[0]*q[0]   + q[1]*q[1] );
    return det < 0;
  }*/

/* Find the elements around a given element whose delaunay circle contains the query point.
   */
std::vector<IndexType> find_violating_elements(const axom::slam::IAMesh<3, 2>& mesh,
                                               const Point2D& query_pt,
                                               IndexType q_element_i)
{
  //SLIC_INFO("find_violating_elements from element " << q_element_i );
  IndexType element_i = q_element_i;
  std::vector<IndexType> ret;
  std::list<IndexType> element_list_to_check;
  element_list_to_check.push_back(element_i);

  std::set<IndexType> checked_elements;
  checked_elements.insert(element_i);
  checked_elements.insert(-1);

  //starting from element_i, which contains the point, try its neighbors, and the neighbor's neighbors
  //for point in circle test

  while(element_list_to_check.size() > 0)
  {
    element_i = element_list_to_check.front();
    element_list_to_check.pop_front();
    if(element_i < 0) continue;

    std::vector<IndexType> verts = mesh.getVertexInElement(element_i);
    Point2D p0 = mesh.getVertexPoint(verts[0]);
    Point2D p1 = mesh.getVertexPoint(verts[1]);
    Point2D p2 = mesh.getVertexPoint(verts[2]);
    //bool is_in_circle = is_point_in_circle(p0, p1, p2, query_pt );
    bool is_in_circle = axom::primal::point_in_circle(query_pt, p0, p1, p2);
    if(is_in_circle)
    {
      ret.push_back(element_i);
      std::vector<IndexType> nbr_elem_list = mesh.getElementNeighbor(element_i);
      for(int i = 0; i < (int)nbr_elem_list.size(); i++)
      {
        IndexType nbr_element = nbr_elem_list[i];

        std::set<IndexType>::iterator it = checked_elements.find(nbr_element);
        if(it == checked_elements.end())
        {
          element_list_to_check.push_back(nbr_element);
          checked_elements.insert(nbr_element);
        }
      }
    }
  }
  return ret;
}
}  // namespace slamTopologicalMesh

/* This is for visualization purpose */
int printc = 1;
void print_out_mesh(const axom::slam::IAMesh<3, 2>& ia_mesh2)
{
  ia_mesh2.print_all();
  using namespace slamTopologicalMesh;
  for(int i = 0; i < ia_mesh2.vcoord_map.actual_size(); i++)
  {
    Point2D p = ia_mesh2.vcoord_map[i];
    std::cout << p[0] << "," << p[1] << ";" << std::endl;
  }
  std::cout << "];\n";
  std::cout << "tri=[";
  for(int i = 0; i < ia_mesh2.element_set.actual_size(); i++)
  {
    if(ia_mesh2.element_set.isValidEntry(i))
    {
      std::vector<int> nlist = ia_mesh2.getVertexInElement(i);
      std::cout << nlist[0] << "," << nlist[1] << "," << nlist[2] << ";"
                << std::endl;
    }
  }
  std::cout << "];\n";
}

int main(/*int argc, char** argv*/)
{
  using namespace slamTopologicalMesh;

  axom::slic::SimpleLogger logger;

  SLIC_INFO("Topological Mesh example");

  // Create an IA triangle mesh
  BasicTriMeshData basic_mesh_data;
  axom::slam::IAMesh<3, 3> ia_mesh(basic_mesh_data.points, basic_mesh_data.tri);

  //test element->vertex function
  std::vector<IndexType> r;
  SLIC_INFO("list of vertexs in element 0");
  r = ia_mesh.getVertexInElement(0);
  for(unsigned int i = 0; i < r.size(); i++)
  {
    SLIC_INFO(r[i]);
  }

  SLIC_INFO("list of elements with vertex 0");
  r = ia_mesh.getElementWithVertex(0);
  for(unsigned int i = 0; i < r.size(); i++)
  {
    SLIC_INFO(r[i]);
  }

  //SLIC_INFO("Calculate normals for the vertexs on the triangle mesh:");
  //calculate_normal(ia_mesh);

  //try it for tetrahedrons
  /*
  SLIC_INFO("Creating tetrahedron mesh");
  BasicTetMeshData basic_tet_mesh;
  axom::slam::IAMesh<4> ia_tetmesh(basic_tet_mesh.points, basic_tet_mesh.tet);

  SLIC_INFO("vertexs in element 0");
  r = ia_tetmesh.getVertexInElement(0);
  for(unsigned int i=0; i<r.size(); i++){
    SLIC_INFO(r[i]);
  }

  SLIC_INFO("elements with vertex 1");
  r = ia_tetmesh.getElementWithVertex(1);
  for(unsigned int i=0; i<r.size(); i++){
    SLIC_INFO(r[i]);
  }*/

  //randomly generate triangle mesh in 2d space
  Flat2DTriMeshData basic_tri_mesh_data;
  axom::slam::IAMesh<3, 2> ia_mesh2(basic_tri_mesh_data.points,
                                    basic_tri_mesh_data.tri);

  //print_out_mesh(ia_mesh2);

  //Adding a number of random points
  srand(1234);
  int num_points = 10;
  IndexType starting_element = 0;  //for the random walk
  for(int pt_i = 0; pt_i < num_points; pt_i++)
  {
    double x = (double)rand() / (double)RAND_MAX * 2.0 - 1.0;
    double y = (double)rand() / (double)RAND_MAX * 2.0 - 1.0;
    double pt_coord[2] = {x, y};
    Point2D new_pt(pt_coord, 2);

    //(assumes a mesh is already existing...)
    //find element containing triangle
    IndexType element_i =
      find_containing_element(ia_mesh2, new_pt, starting_element);

    //find the delaunay cavity...
    std::vector<IndexType> elements_to_remove =
      find_violating_elements(ia_mesh2, new_pt, element_i);
    std::cout << "violating elements:";
    for(unsigned int i = 0; i < elements_to_remove.size(); i++)
      std::cout << elements_to_remove[i] << " ";
    std::cout << std::endl;

    //remove each elements from the mesh, keeping a list of cavity edges and their connected elements
    typedef std::pair<IndexType, IndexType> IndexPairType;
    typedef std::pair<IndexPairType, IndexType> CavityPairType;
    typedef std::map<IndexPairType, IndexType> CavityMapType;
    CavityMapType cavity_edges;
    for(unsigned int i = 0; i < elements_to_remove.size(); i++)
    {
      IndexType element_i = elements_to_remove[i];
      std::vector<IndexType> verts_in_this_element =
        ia_mesh2.getVertexInElement(element_i);
      for(unsigned int j = 0; j < verts_in_this_element.size(); j++)
      {
        IndexType n1 = verts_in_this_element[(j * 2) % 3];
        IndexType n2 = verts_in_this_element[(1 + j * 2) % 3];
        CavityMapType::iterator iter = cavity_edges.find(IndexPairType(n2, n1));
        if(iter == cavity_edges.end())
        {
          cavity_edges.insert(CavityPairType(IndexPairType(n1, n2), element_i));
        }
        else
        {  //already exists
          cavity_edges.erase(iter);
        }
      }
      ia_mesh2.removeElement(element_i);
    }

    SLIC_INFO("After removing elements");
    //print_out_mesh(ia_mesh2);
    //std::cout<<"pt = [" << pt_coord[0]<<" "<<pt_coord[1] <<"];" << std::endl;

    //new triangles from the cavity edges

    IndexType new_pt_i = ia_mesh2.addVertex(new_pt);
    for(CavityMapType::iterator iter = cavity_edges.begin();
        iter != cavity_edges.end();
        iter++)
    {
      IndexType n1 = iter->first.first;
      IndexType n2 = iter->first.second;
      //IndexType nbr_element_i = iter->second;

      starting_element = ia_mesh2.addElement(n1, n2, new_pt_i);
      SLIC_INFO("New cell: " << starting_element << " with " << n1 << ", " << n2
                             << ", " << new_pt_i);
    }
    SLIC_INFO("After adding elements");
    //print_out_mesh(ia_mesh2);
    //std::cout<<"pt = [" << pt_coord[0]<<" "<<pt_coord[1] <<"];"<< std::endl;

    if(!ia_mesh2.isValid(true))
      SLIC_INFO("Something wrong with Delaunay on IA.");
  }

  print_out_mesh(ia_mesh2);

  SLIC_INFO("Done!");

  return 0;
}
