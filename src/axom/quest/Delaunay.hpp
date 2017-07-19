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

#define TOPO_DIM 2
#define VERT_PER_ELEMENT 3
using PointType = primal::Point<DataType, 2>;
using IATriMeshType = slam::IAMesh<TOPO_DIM, 2, PointType>;
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

void Delaunay::printMesh() { mesh.print_all(); }

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

  // typedef for map structures for cavity edges and new elements
  //cavity_edge_map: <vlist, <vlist, element_idx> >
  // cavity_edge_map is a map from sorted face vertices to pair of face vertices and element idx.
  // It keeps track of the edges that forms the cavity, which will be edges of the new element.
  //new_element_face_map: <vlist, <element_idx,face_idx> >
  // new_element_face_map is a map from sorted face vertices to pair of element and its face idx.
  // It keep tracks of the new elements added, to correctly construct element->element relation

  typedef axom::slam::IndexListType IndexListType;
  typedef std::pair<IndexType, IndexType> IndexPairType;

  typedef std::map<IndexListType, IndexListType> CavityMapType;
  typedef std::pair<IndexListType, IndexListType> CavityMapPairType;

  typedef std::pair<IndexListType, IndexPairType> NewElemFacePairType;
  typedef std::map<IndexListType, IndexPairType> NewElemFaceMapType;

  CavityMapType cavity_edge_map;
  NewElemFaceMapType new_element_face_map;

  // Take out the elements one by one, keeping track of the cavity edges
  for(unsigned int i = 0; i < elements_to_remove.size(); i++)
  {
    IndexType element_i = elements_to_remove[i];

    //For each face of the element, insert its face into the map (or remove if duplicate)
    for(unsigned int face_i = 0; face_i < VERT_PER_ELEMENT; face_i++)
    {
      IndexListType face_vlist = mesh.getElementFace(element_i, face_i);

      IndexListType face_vlist_sorted(
        face_vlist);  //sort the list to make the key for the map
      std::sort(face_vlist_sorted.begin(), face_vlist_sorted.end());

      std::pair<CavityMapType::iterator, bool> ret =
        cavity_edge_map.insert(CavityMapPairType(face_vlist_sorted, face_vlist));
      if(!ret.second)
      {
        cavity_edge_map.erase(ret.first);
      }
    }

    mesh.removeElement(element_i);
    num_removed_elements_since_last_compact++;
  }

  //SLIC_INFO("DT: After removing elements");
  //mesh.print_all();

  //Add the new point
  IndexType new_pt_i = mesh.addVertex(new_pt);

  //Add new triangles from the cavity edges
  for(CavityMapType::iterator cav_edge_iter = cavity_edge_map.begin();
      cav_edge_iter != cavity_edge_map.end();
      cav_edge_iter++)
  {
    //SLIC_INFO("new element " << new_el <<": "<<cav_edge_iter->second.first[0]<<" "<<cav_edge_iter->second.first[1] << " " <<new_pt_i);

    IndexType new_el = mesh.addElement(cav_edge_iter->second[0],
                                       cav_edge_iter->second[1],
                                       new_pt_i);

    // Check the new element's face neighbor, because it can be wrong sometimes
    // due to the mesh being temporarily non-manifold.
    for(int fi = 0; fi < VERT_PER_ELEMENT; fi++)
    {
      IndexListType vlist_sorted = mesh.getElementFace(new_el, fi);
      std::sort(vlist_sorted.begin(), vlist_sorted.end());

      NewElemFaceMapType::iterator iter = new_element_face_map.find(vlist_sorted);

      if(iter == new_element_face_map.end())
      {
        new_element_face_map.insert(
          NewElemFacePairType(vlist_sorted, IndexPairType(new_el, fi)));
      }
      else
      {
        IndexType nbr_el = iter->second.first;
        IndexType nbr_fi = iter->second.second;

        if(mesh.ee_rel[new_el][fi] < 0)
        {
          //SLIC_INFO("found neighbor " << iter->second.first <<"-"<<nbr_fi);

          mesh.ee_rel.modify(new_el, fi, nbr_el);
          mesh.ee_rel.modify(nbr_el, nbr_fi, new_el);
        }

        new_element_face_map.erase(iter);
      }
    }
  }

  //SLIC_INFO("DT: After adding elements");
  //mesh.print_all();

  //Make sure ee_rel has only 4 invalid
  int invalid_nbr_count = 0;
  for(int i = 0; i < mesh.ee_rel.size(); i++)
  {
    if(mesh.ee_rel.isValidEntry(i))
    {
      for(int j = 0; j < VERT_PER_ELEMENT; j++)
      {
        invalid_nbr_count += mesh.ee_rel[i][j] < 0;
      }
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

  if(!mesh.isValid(true))
  {
    SLIC_INFO("IA mesh is invalid after adding new point");
  }
}

IndexType Delaunay::findContainingElement(PointType query_pt)
{
  //SLIC_INFO("DT: find_containing_element " << element_i );
  //SLIC_INFO("Query Pt " << query_pt);

  SLIC_ASSERT(!mesh.isEmpty());

  //find the last valid element to use as starting element
  IndexType element_i;
  for(int i = mesh.element_set.size() - 1; true; i--)
  {
    if(mesh.isValidElementEntry(i))
    {
      element_i = i;
      break;
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
    {                    // bary_val bigger than zero -> inside triangle
      return element_i;  //return if inside or on triangle
    }

    std::vector<IndexType> zlist = mesh.getElementNeighbor(element_i);
    IndexType next_el = element_i = zlist[(i + 1) % 3];

    if(next_el < 0)
    {
      SLIC_WARNING(
        "Either there is a hole in the mesh, or the point is outside of the "
        "mesh.");
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
