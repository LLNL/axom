// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file Delaunay.hpp
 *
 * \brief Construct delaunay triangulation by inserting points one by one.
 * A bounding box of the points needs to be defined first via startWithBoundary(...)
 */

#ifndef QUEST_DELAUNAY_H_
#define QUEST_DELAUNAY_H_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"  //for writing out to VTK file

#include <list>
#include <vector>
#include <set>
#include <stdlib.h>

namespace axom
{
namespace quest
{
using DataType = double;
using IndexType = int;

using PointType = primal::Point<DataType, 2>;
using IATriMeshType = slam::IAMesh<2, 2, PointType>;
using Triangle2D = axom::primal::Triangle<double, 2>;

class Delaunay
{
public:
  Delaunay() { }
  void startWithBoundary(DataType xmin, DataType xmax, DataType ymin, DataType ymax);
  void insertPoint(const PointType& pt);
  void printMesh();
  void writeToVTKFile(const char* filename);
  void removeBoundary();

private:
  IndexType findContainingElement(PointType pt);
  std::vector<IndexType> findViolatingElements(PointType pt, IndexType i);
  void createCavity();
  void delaunayBall();

  IATriMeshType mesh;
  IndexType newest_element_i = 0;
  bool has_boundary = false;
  int num_removed_elements_since_last_compact = 0;
};

void Delaunay::startWithBoundary(DataType xmin,
                                 DataType xmax,
                                 DataType ymin,
                                 DataType ymax)
{
  std::vector<DataType> points = {xmin, ymin, xmin, ymax, xmax, ymin, xmax, ymax};

  std::vector<IndexType> tri = {0, 2, 1, 3, 1, 2};

  mesh = IATriMeshType(points, tri);

  has_boundary = true;
}

void Delaunay::insertPoint(const PointType& new_pt)
{
  //Make sure startWithBoundary(...) is called first
  SLIC_ASSERT_MSG(
    has_boundary,
    "Error: Need a predefined boundary box prior to adding points.");

  //Make sure the new point is inside the boundary box
  SLIC_ASSERT_MSG(
    mesh.vcoord_map[0][0] < new_pt[0] && mesh.vcoord_map[3][0] > new_pt[0] &&
      mesh.vcoord_map[0][1] < new_pt[1] && mesh.vcoord_map[3][1] > new_pt[1],
    "Error: new point is outside of the boundary box.");

  //mesh.print_all();
  //SLIC_INFO("DT: Inserting point "<<new_pt[0] <<"," <<new_pt[1]);

  IndexType element_i = findContainingElement(new_pt);

  //SLIC_INFO("DT: Containing element: "<<element_i );

  std::vector<IndexType> elements_to_remove =
    findViolatingElements(new_pt, element_i);

  /*
  std::cout<<"DT: violating elements:" ;
  for(unsigned int i=0; i<elements_to_remove.size(); i++)
    std::cout<<elements_to_remove[i] <<" ";
  std::cout<<std::endl;
  //*/

  SLIC_ASSERT_MSG(elements_to_remove.size() > 0,
                  "Error: New point is not contained in the mesh");

  //remove each elements from the mesh, keeping a list of cavity edges and their connected elements
  typedef std::pair<IndexType, IndexType> IndexPairType;
  typedef std::pair<IndexPairType, IndexType> CavityPairType;
  typedef std::map<IndexPairType, IndexType> CavityMapType;
  CavityMapType cavity_edges;
  std::vector<IndexPairType> insertOrder;  //to keep track of the insertion order.
    // This is to prevent multiple (temporary) holes in the mesh, which will result
    // in wrong functionality of IA's getElementWithVertex(v_id)
  for(unsigned int i = 0; i < elements_to_remove.size(); i++)
  {
    IndexType element_i = elements_to_remove[i];
    std::vector<IndexType> verts_in_this_element =
      mesh.getVertexInElement(element_i);
    for(unsigned int j = 0; j < verts_in_this_element.size(); j++)
    {
      IndexType n1 = verts_in_this_element[(j * 2) % 3];
      IndexType n2 = verts_in_this_element[(1 + j * 2) % 3];
      CavityMapType::iterator iter = cavity_edges.find(IndexPairType(n2, n1));
      if(iter == cavity_edges.end())
      {
        cavity_edges.insert(CavityPairType(IndexPairType(n1, n2), element_i));
        insertOrder.push_back(IndexPairType(n1, n2));
      }
      else
      {  //already exists
        cavity_edges.erase(iter);
      }
    }
    mesh.removeElement(element_i);
    num_removed_elements_since_last_compact++;
  }

  //SLIC_INFO("DT: After removing elements");
  //print_out_mesh(mesh);
  //std::cout<<"pt = [" << pt_coord[0]<<" "<<pt_coord[1] <<"];" << std::endl;

  //new triangles from the cavity edges

  IndexType new_pt_i = mesh.addVertex(new_pt);

  /*// Replaced by the code below that reduces number of 2+ holes in the mesh
  for(CavityMapType::iterator iter = cavity_edges.begin(); iter != cavity_edges.end(); iter++)
  {
    IndexType n1 = iter->first.first;
    IndexType n2 = iter->first.second;

    IndexType starting_element = mesh.addElement(n1, n2, new_pt_i);
    //SLIC_INFO("DT: New cell: " << starting_element <<" with "<< n1 <<", " << n2 << ", "<<new_pt_i );

    this->newest_element_i= starting_element;
  }*/

  std::vector<IndexType> new_cells;

  for(unsigned int i = 0; i < insertOrder.size(); i++)
  {
    CavityMapType::iterator iter = cavity_edges.find(insertOrder[i]);
    if(iter != cavity_edges.end())
    {
      IndexType n1 = iter->first.first;
      IndexType n2 = iter->first.second;

      IndexType new_element = mesh.addElement(n1, n2, new_pt_i);
      //SLIC_INFO("DT: New cell: " << starting_element <<" with "<< n1 <<", " << n2 << ", "<<new_pt_i );
      new_cells.push_back(new_element);
      this->newest_element_i = new_element;
    }
  }

  //SLIC_INFO("DT: After adding elements");
  //mesh.print_all();

  // Here is some ugly code that checks and fixes element neighbors that are incorrect.
  // This is because IA::getElementWithVertex() only guarantee to work on manifold mesh, which
  // can result in incorrect ee_rel.
  bool all_element_valid = false;
  //Check for incorrect element
  while(!all_element_valid)
  {
    all_element_valid = true;

    for(unsigned int i = 0; i < new_cells.size(); i++)
    {
      IndexType new_el = new_cells[i];

      int invalid_nbr_count = 0;
      for(int j = 0; j < mesh.VERTS_PER_ELEM; j++)
      {
        invalid_nbr_count += mesh.ee_rel[new_el][j] < 0;
      }

      if(invalid_nbr_count > 0)
      {
        int num_boundary_box_vert = 0;
        for(int j = 0; j < mesh.VERTS_PER_ELEM; j++)
        {
          num_boundary_box_vert += mesh.ev_rel[new_el][j] < 4;
        }
        if(invalid_nbr_count == 1 &&
           num_boundary_box_vert == mesh.VERTS_PER_ELEM - 1)
        {
          //it's a correct boundary element
        }
        else
        {
          //Reinsert this element
          all_element_valid = false;
          axom::slam::IndexListType vlist = mesh.getVertexInElement(new_el);
          mesh.removeElement(new_el);
          new_cells[i] = mesh.addElement(vlist[0], vlist[1], vlist[2]);
        }
      }
    }

  }  //end while, end ugly code

  //Make sure ee_rel has only 4 invalid
  int invalid_nbr_count = 0;
  for(int i = 0; i < mesh.ee_rel.size(); i++)
  {
    if(mesh.ee_rel.isValidEntry(i))
      for(int j = 0; j < 3; j++)
      {
        invalid_nbr_count += mesh.ee_rel[i][j] < 0;
      }
  }
  SLIC_ASSERT(invalid_nbr_count == 4);

  //call compact() if there are too many invalid points
  if(num_removed_elements_since_last_compact > 64 &&
     num_removed_elements_since_last_compact > (mesh.element_set.size() / 2))
  {
    mesh.compact();
    num_removed_elements_since_last_compact = 0;
  }

  /*if(! mesh.isValid(true) )
    SLIC_INFO("Something wrong with Delaunay on IA.");*/
}

IndexType Delaunay::findContainingElement(PointType query_pt)
{
  IndexType element_i = newest_element_i;

  //SLIC_INFO("DT: find_containing_element " << element_i );
  //assumes the point is definitely contained within a triangle
  //TODO fix that assumption ^
  //SLIC_INFO("Query Pt " << query_pt);

  if(!mesh.isValidElementEntry(element_i))
  {
    if(mesh.isEmpty() == 0)
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
    //SLIC_INFO("DT: step into element"<<element_i);
    std::vector<IndexType> verts = mesh.getVertexInElement(element_i);

    Triangle2D tri(mesh.getVertexPoint(verts[0]),
                   mesh.getVertexPoint(verts[1]),
                   mesh.getVertexPoint(verts[2]));

    axom::primal::Point3D bary_co = tri.barycentricCoords(query_pt);
    //SLIC_INFO("bary_co " << bary_co);

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

std::vector<IndexType> Delaunay::findViolatingElements(PointType query_pt,
                                                       IndexType element_i)
{
  //SLIC_INFO("find_violating_elements from element " << q_element_i );

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
    PointType p0 = mesh.getVertexPoint(verts[0]);
    PointType p1 = mesh.getVertexPoint(verts[1]);
    PointType p2 = mesh.getVertexPoint(verts[2]);
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

void Delaunay::writeToVTKFile(const char* filename)
{
  mesh.compact();
  newest_element_i = mesh.element_set.size() - 1;

  axom::mint::UnstructuredMesh<MINT_TRIANGLE> mint_mesh(2);
  for(int i = 0; i < mesh.vertex_set.size(); i++)
  {
    if(!mesh.vertex_set.isValidEntry(i)) continue;
    mint_mesh.insertNode(mesh.vcoord_map[i][0], mesh.vcoord_map[i][1]);
  }
  for(int i = 0; i < mesh.ev_rel.size(); i++)
  {
    if(!mesh.ev_rel.isValidEntry(i)) continue;
    std::vector<int> vec = mesh.ev_rel[i];
    const int* ptr = vec.data();
    mint_mesh.insertCell(ptr, MINT_TRIANGLE, 3);
  }
  mint_mesh.toVtkFile(filename);
}

void Delaunay::removeBoundary()
{
  if(has_boundary)
  {
    //remove the boundary box, which will be the first 4 points.
    mesh.removeVertex(0);
    mesh.removeVertex(1);
    mesh.removeVertex(2);
    mesh.removeVertex(3);

    has_boundary = false;
  }
}

}  // end namespace quest
}  // end namespace axom

#endif  //  QUEST_DELAUNAY_H_
