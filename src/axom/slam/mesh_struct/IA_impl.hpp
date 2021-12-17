// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SLAM_IA_IMPL_H_
#define SLAM_IA_IMPL_H_

/*
 * \file IA_impl.hpp
 *
 * \brief Contains the implementation of the IAMesh class and helper functions
 */

#include "axom/slam/policies/SizePolicies.hpp"
#include "axom/slam/ModularInt.hpp"

#include "axom/fmt.hpp"

#include <vector>
#include <map>

namespace axom
{
namespace slam
{
/**
 * Helper function on std::vectors
 */

namespace /*anonymous*/
{
// Checks if v is in the list s
template <typename T>
bool is_subset(T v, const std::vector<T>& s)
{
  const int SZ = s.size();
  for(int i = 0; i < SZ; ++i)
  {
    if(s[i] == v)
    {
      return true;
    }
  }
  return false;
}

}  //end anonymous namespace

template <unsigned int TDIM, unsigned int SDIM, typename P>
typename IAMesh<TDIM, SDIM, P>::ElementAndFaceIdxType
IAMesh<TDIM, SDIM, P>::ElemNbrFinder(V2EMapType& vertpair_to_elem_map,
                                     IndexType element_i,
                                     IndexType side_i)
{
  // NOTE: V2EMapType maps a sorted tuple of vertex IDs to a face on a given
  //       mesh element. It is used to find the element index of the opposite
  //       face within the mesh

  IndexArray vlist = getElementFace(element_i, side_i);
  std::sort(vlist.begin(), vlist.end());

  ElementAndFaceIdxType zs_pair(element_i, side_i);

  auto map_ret2 = vertpair_to_elem_map.insert(std::make_pair(vlist, zs_pair));
  if(!map_ret2.second)  //if this pair is in the map, we've found our match
  {
    auto orig_pair = map_ret2.first->second;
    vertpair_to_elem_map.erase(map_ret2.first);
    return orig_pair;
  }

  //No matching pair is found. Return an invalid pair
  return ElementAndFaceIdxType(-1, -1);
}

template <unsigned int TDIM, unsigned int SDIM, typename P>
void IAMesh<TDIM, SDIM, P>::print_all() const
{
  axom::fmt::memory_buffer out;
  axom::fmt::format_to(out,
                       "IA mesh: {} mesh in {}d with {} valid vertices (of {}) "
                       "and {} valid elements (of {})\n",
                       (TDIM == 2 ? "triangle" : "tetrahedral"),
                       SDIM,
                       vertex_set.numberOfValidEntries(),
                       vertex_set.size(),
                       element_set.numberOfValidEntries(),
                       element_set.size());

  //print out the sets, relation, and map
  axom::fmt::format_to(out,
                       "  element_set ({}/{}): [{}]\n",
                       element_set.numberOfValidEntries(),
                       element_set.size(),
                       axom::fmt::join(element_set, ", "));

  axom::fmt::format_to(out,
                       "  vertex_set ({}/{}): [{}]\n",
                       vertex_set.numberOfValidEntries(),
                       vertex_set.size(),
                       axom::fmt::join(vertex_set, ", "));

  {
    const int sz = ev_rel.size();
    std::vector<std::string> strs(sz);
    for(auto pos : element_set.positions())
    {
      strs[pos] = element_set.isValidEntry(pos)
        ? axom::fmt::format("{}: {}", pos, ev_rel[pos])
        : axom::fmt::format("{}: {{}}", pos);
    }
    axom::fmt::format_to(out,
                         "  ev_rel ({}/{}): [{}]\n",
                         ev_rel.numberOfValidEntries(),
                         sz,
                         axom::fmt::join(strs, "; "));
  }

  {
    const int sz = ve_rel.size();
    std::vector<std::string> strs(sz);
    for(auto pos : vertex_set.positions())
    {
      strs[pos] = vertex_set.isValidEntry(pos)
        ? axom::fmt::format("{}: {}", pos, ve_rel[pos])
        : axom::fmt::format("{}: {{}}", pos);
    }
    axom::fmt::format_to(out,
                         "  ve_rel ({}/{}): [{}]\n",
                         ve_rel.numberOfValidEntries(),
                         sz,
                         axom::fmt::join(strs, "; "));
  }

  {
    const int sz = ee_rel.size();
    std::vector<std::string> strs(sz);
    for(auto pos : element_set.positions())
    {
      strs[pos] = element_set.isValidEntry(pos)
        ? axom::fmt::format("{}: {}", pos, ee_rel[pos])
        : axom::fmt::format("{}: {{}}", pos);
    }
    axom::fmt::format_to(out,
                         "  ee_rel ({}/{}): [{}]\n",
                         ee_rel.numberOfValidEntries(),
                         sz,
                         axom::fmt::join(strs, "; "));
  }

  {
    const int sz = vcoord_map.size();
    std::vector<std::string> strs(sz);
    for(auto pos : vertex_set.positions())
    {
      strs[pos] = vcoord_map.isValidEntry(pos)
        ? axom::fmt::format("{}: {}", pos, vcoord_map[pos])
        : axom::fmt::format("{}: --", pos);
    }
    axom::fmt::format_to(out,
                         "  vertex coord ({}/{}): [{}]\n",
                         vcoord_map.numberOfValidEntries(),
                         sz,
                         axom::fmt::join(strs, "; "));
  }

  SLIC_INFO(axom::fmt::to_string(out));
}

/*********************************************************************************/

template <unsigned int TDIM, unsigned int SDIM, typename P>
IAMesh<TDIM, SDIM, P>::IAMesh()
  : vertex_set(0)
  , element_set(0)
  , ev_rel(&element_set, &vertex_set)
  , ve_rel(&vertex_set, &element_set)
  , ee_rel(&element_set, &element_set)
  , vcoord_map(&vertex_set)
{ }

template <unsigned int TDIM, unsigned int SDIM, typename P>
IAMesh<TDIM, SDIM, P>::IAMesh(const IAMesh& m)
{
  operator=(m);
}

template <unsigned int TDIM, unsigned int SDIM, typename P>
IAMesh<TDIM, SDIM, P>& IAMesh<TDIM, SDIM, P>::operator=(const IAMesh& m)
{
  if(&m != this)
  {
    vertex_set = m.vertex_set;
    element_set = m.element_set;
    ev_rel = ElementBoundaryRelation(&element_set, &vertex_set);
    ve_rel = VertexCoboundaryRelation(&vertex_set, &element_set);
    ee_rel = ElementAdjacencyRelation(&element_set, &element_set);
    vcoord_map = PositionMap(&vertex_set);

    ev_rel.data() = m.ev_rel.data();
    ve_rel.data() = m.ve_rel.data();
    ee_rel.data() = m.ee_rel.data();
    vcoord_map.data() = m.vcoord_map.data();
  }

  return *this;
}

template <unsigned int TDIM, unsigned int SDIM, typename P>
IAMesh<TDIM, SDIM, P>::IAMesh(std::vector<double>& points,
                              std::vector<IndexType>& tri)
  : vertex_set(points.size() / COORDS_PER_VERT)
  , element_set(tri.size() / VERTS_PER_ELEM)
  , ev_rel(&element_set, &vertex_set)
  , ve_rel(&vertex_set, &element_set)
  , ee_rel(&element_set, &element_set)
  , vcoord_map(&vertex_set)
{
  // Relation, element to vertex boundary relation
  for(IndexType idx = 0; idx < element_set.size() * VERTS_PER_ELEM; ++idx)
  {
    ev_rel.insert(idx / VERTS_PER_ELEM, tri[idx]);
  }
  SLIC_ASSERT_MSG(
    ev_rel.isValid(),
    "Error creating (dynamic) relation from elements to vertices!");

  // The map, vertex to coordinates
  for(IndexType idx = 0; idx < vertex_set.size(); ++idx)
  {
    vcoord_map[idx] = Point(&(points[idx * COORDS_PER_VERT]));
  }
  SLIC_ASSERT_MSG(vcoord_map.isValid(true),
                  "Error creating map from vertex to coords!");

  //Vertex element relation. 1->1 mapping only 1 element per vertex.
  for(IndexType zIdx = 0; zIdx < element_set.size(); ++zIdx)
  {
    for(IndexType idx = 0; idx < (int)ev_rel[zIdx].size(); ++idx)
    {
      ve_rel.modify(ev_rel[zIdx][idx], 0, zIdx);
    }
  }
  SLIC_ASSERT_MSG(
    ve_rel.isValid(true),
    "Error creating (dynamic) relation from vertices to elements!\n");

  //Before making element to element relation, construct the data.
  // For every cell, find the union of triangles for each pair of vertices
  IndexBuf::BufferType element_element_vec =
    index_buffer.addBuffer("element_element_vector",
                           element_set.size() * VERTS_PER_ELEM);
  for(unsigned int i = 0; i < element_element_vec.size(); i++)
  {
    //initialize to no neighbor element
    element_element_vec[i] = ElementBoundaryRelation::INVALID_INDEX;
  }

  V2EMapType vertpair_to_elem_map;

  for(IndexType element_i = 0; element_i < element_set.size(); element_i++)
  {
    for(IndexType side_i = 0; side_i < VERTS_PER_ELEM; side_i++)
    {
      ElementAndFaceIdxType nst =
        ElemNbrFinder(vertpair_to_elem_map, element_i, side_i);

      IndexType other_element_idx = nst.first;
      IndexType other_side_idx = nst.second;

      if(!element_set.isValidEntry(other_element_idx)) continue;

      int idx0 = element_i * VERTS_PER_ELEM + side_i;
      element_element_vec[idx0] = other_element_idx;

      int idx1 = other_element_idx * VERTS_PER_ELEM + other_side_idx;
      element_element_vec[idx1] = element_i;
    }
  }

  //Element adjacency relation along facets
  for(int i = 0; i < element_set.size(); i++)
  {
    for(int j = 0; j < VERTS_PER_ELEM; j++)
    {
      ee_rel.modify(i, j, element_element_vec[i * VERTS_PER_ELEM + j]);
    }
  }
  SLIC_ASSERT_MSG(
    ee_rel.isValid(true),
    "Error creating (dynamic) relation from elements to elements!");
}

template <unsigned int TDIM, unsigned int SDIM, typename P>
typename IAMesh<TDIM, SDIM, P>::IndexArray
IAMesh<TDIM, SDIM, P>::getVerticesInElement(IndexType element_idx) const
{
  IndexArray ret;
  if(!ev_rel.isValidEntry(element_idx))
  {
    SLIC_WARNING("Attempting to retrieve data with an invalid element");
    return ret;
  }

  typename ElementBoundaryRelation::RelationSubset rvec = ev_rel[element_idx];
  ret.resize(rvec.size());
  for(int i = 0; i < rvec.size(); i++)
  {
    ret[i] = rvec[i];
  }

  return ret;
}

template <unsigned int TDIM, unsigned int SDIM, typename P>
typename IAMesh<TDIM, SDIM, P>::IndexArray
IAMesh<TDIM, SDIM, P>::getElementsWithVertex(IndexType vertex_idx) const
{
  IndexArray ret;

  if(!ve_rel.isValidEntry(vertex_idx))
  {
    //this vertex is not connected to any elements
    SLIC_WARNING_IF(
      !vertex_set.isValidEntry(vertex_idx),
      "Attempting to retrieve data with an invalid vertex id: " << vertex_idx);
    return ret;
  }

  IndexType starting_element_idx = ve_rel[vertex_idx][0];

  ret.push_back(starting_element_idx);
  IndexArray element_traverse_queue;
  element_traverse_queue.push_back(starting_element_idx);

  while(element_traverse_queue.size() > 0)
  {
    IndexType element_i = element_traverse_queue.back();
    element_traverse_queue.pop_back();

    typename ElementAdjacencyRelation::RelationSubset surrounded_elements =
      ee_rel[element_i];

    for(int i = 0; i < surrounded_elements.size(); i++)
    {
      IndexType ele_i = surrounded_elements[i];

      if(ele_i == ElementAdjacencyRelation::INVALID_INDEX ||
         is_subset(ele_i, ret))
        continue;

      // If this element contains the vertex in question, add it to
      // the return vector and enqueue to check neighbors.
      for(int j = 0; j < ev_rel[ele_i].size(); j++)
      {
        IndexType vert_i = ev_rel[ele_i][j];
        if(vert_i == vertex_idx)
        {
          ret.push_back(ele_i);
          element_traverse_queue.push_back(ele_i);
          break;
        }
      }
    }
  }

  /*
     //This code is for confirmation purpose, since it's a known issue that if
     // the mesh is non-manifold, this list may be incomplete.
     // This code traverses the whole relation list to find the complete
     // list of elements in the star of this vertex.
     std::vector<IndexType> elem_list;
     for(int elem_i=0; elem_i<(int)ev_rel.size(); elem_i++){
     for(int j=0; j<ev_rel[elem_i].size(); j++){
      if( ev_rel[elem_i][j] == vertex_idx ){
        elem_list.push_back(elem_i);
        break;
      }
     }
     }

     //check if the neighbor-walk list is missing elements from complete list
     for(int i=0; i<(int)elem_list.size(); i++){
       IndexType elem = elem_list[i];
       if( !is_subset(elem, ret) && ee_rel.isValidEntry(elem)){
        SLIC_INFO(
          "!!!! Element "<< elem
          << " is missing from nbr list of vertex "<<vertex_idx);
       }
     }
     //*/

  //return elem_list; //complete list

  return ret;
}

template <unsigned int TDIM, unsigned int SDIM, typename P>
typename IAMesh<TDIM, SDIM, P>::IndexArray IAMesh<TDIM, SDIM, P>::getElementFace(
  IndexType element_idx,
  IndexType face_idx) const
{
  using CTSize = slam::policies::CompileTimeSize<IndexType, VERTS_PER_ELEM>;
  slam::ModularInt<CTSize> mod_face(face_idx);

  IndexArray ret;

  if(!element_set.isValidEntry(element_idx))
  {
    SLIC_WARNING(
      "Attempting to retrieve data with an invalid element: " << element_idx);

    return ret;
  }

  SLIC_ASSERT_MSG(0 <= face_idx && face_idx < VERTS_PER_ELEM,
                  "Face index is invalid.");

  SLIC_ASSERT(ev_rel[element_idx].size() == VERTS_PER_ELEM);

  for(int i = 0; i < VERTS_PER_ELEM - 1; i++)
  {
    ret.push_back(ev_rel[element_idx][mod_face + i]);
  }

  return ret;
}

template <unsigned int TDIM, unsigned int SDIM, typename P>
typename IAMesh<TDIM, SDIM, P>::IndexArray
IAMesh<TDIM, SDIM, P>::getElementNeighbors(IndexType element_idx) const
{
  IndexArray ret;

  if(!ee_rel.isValidEntry(element_idx))
  {
    //this element is invalid
    SLIC_WARNING("Attempting to retrieve data with an invalid element.");
    return ret;
  }

  typename ElementBoundaryRelation::RelationSubset rvec = ee_rel[element_idx];
  ret.resize(rvec.size());
  for(int i = 0; i < rvec.size(); i++)
  {
    ret[i] = rvec[i];
  }

  return ret;
}

template <unsigned int TDIM, unsigned int SDIM, typename P>
const typename IAMesh<TDIM, SDIM, P>::Point& IAMesh<TDIM, SDIM, P>::getVertexPoint(
  IndexType vertex_idx) const
{
  SLIC_ASSERT(isValidVertexEntry(vertex_idx));

  return vcoord_map[vertex_idx];
}

template <unsigned int TDIM, unsigned int SDIM, typename P>
void IAMesh<TDIM, SDIM, P>::removeVertex(IndexType vertex_idx)
{
  if(!vertex_set.isValidEntry(vertex_idx))
  {
    SLIC_WARNING("Attempting to remove an invalid vertex");
    return;
  }

  //check if any element uses this vertex. If so, remove them too.
  IndexArray attached_elements = getElementsWithVertex(vertex_idx);
  for(int i = 0; i < (int)attached_elements.size(); i++)
  {
    removeElement(attached_elements[i]);
  }

  vertex_set.remove(vertex_idx);
  ve_rel.remove(vertex_idx);
  //Note: once the set entry is removed, its corresponding
  // map entry is assumed to be invalid
}

template <unsigned int TDIM, unsigned int SDIM, typename P>
void IAMesh<TDIM, SDIM, P>::removeElement(IndexType element_idx)
{
  if(!element_set.isValidEntry(element_idx))
  {
    SLIC_WARNING("Attempting to remove an invalid element");
    return;
  }

  //check if ve_rel needs to be updated for the vertices of the removed cell
  for(int i = 0; i < ev_rel[element_idx].size(); i++)
  {
    IndexType vertex_i = ev_rel[element_idx][i];

    IndexArray element_with_this_vertex = getElementsWithVertex(vertex_i);

    if(element_with_this_vertex.size() == 1)
    {  //the element being removed is the last element using this vertex
      ve_rel.modify(vertex_i, 0, VertexCoboundaryRelation::INVALID_INDEX);
    }
    else
    {
      //there are other elements using this vertex.
      // update ve_rel to not use the removed element
      if(ve_rel[vertex_i][0] == element_idx)
      {
        bool modified = false;
        for(unsigned int j = 0; j < element_with_this_vertex.size(); j++)
        {
          if(element_with_this_vertex[j] != element_idx &&
             element_set.isValidEntry(element_with_this_vertex[j]))
          {
            ve_rel.modify(vertex_i, 0, element_with_this_vertex[j]);
            modified = true;
            break;
          }
        }
        SLIC_ASSERT(modified);
      }
    }
  }

  //erase this element's data
  element_set.remove(element_idx);
  ev_rel.remove(element_idx);

  //erase neighbor element's data
  for(int zi = 0; zi < ee_rel[element_idx].size(); zi++)
  {
    IndexType nbr_element_idx = ee_rel[element_idx][zi];
    if(nbr_element_idx < 0) continue;
    for(int i = 0; i < ee_rel[nbr_element_idx].size(); i++)
    {
      if(ee_rel[nbr_element_idx][i] == element_idx)
      {
        ee_rel.modify(nbr_element_idx, i, ElementBoundaryRelation::INVALID_INDEX);
        break;
      }
    }
  }
  ee_rel.remove(element_idx);
}

template <unsigned int TDIM, unsigned int SDIM, typename P>
typename IAMesh<TDIM, SDIM, P>::IndexType IAMesh<TDIM, SDIM, P>::addVertex(
  const Point& p)
{
  IndexType vertex_idx = vertex_set.insert();
  vcoord_map.insert(vertex_idx, p);
  ve_rel.insert(vertex_idx, VertexCoboundaryRelation::INVALID_INDEX);

  return vertex_idx;
}

template <unsigned int TDIM, unsigned int SDIM, typename P>
typename IAMesh<TDIM, SDIM, P>::IndexType IAMesh<TDIM, SDIM, P>::addElement(
  IndexType n0,
  IndexType n1,
  IndexType n2,
  IndexType n3)
{
  SLIC_ASSERT(VERTS_PER_ELEM <= 4);
  IndexType nlist[] = {n0, n1, n2, n3};
  return addElement(nlist);
}

template <unsigned int TDIM, unsigned int SDIM, typename P>
typename IAMesh<TDIM, SDIM, P>::IndexType IAMesh<TDIM, SDIM, P>::addElement(
  const IndexType* nlist)
{
  // Implementation note:
  //   This function reconstructs the vertex-element relation
  //    on each vertex ID of the new element
  // Can we optimize this function?

  for(int i = 0; i < VERTS_PER_ELEM; i++)
  {
    if(!vertex_set.isValidEntry(nlist[i]))
      SLIC_WARNING(
        "Trying to add an element with invalid vertex index:" << nlist[i]);
  }

  IndexType element_idx = element_set.insert();

  for(int i = 0; i < VERTS_PER_ELEM; i++) ev_rel.insert(element_idx, nlist[i]);

  //make sure the space is allocated in ee_rel
  ee_rel.insert(element_idx, ElementAdjacencyRelation::INVALID_INDEX);

  V2EMapType vertpair_to_elem_map;

  //First add each face of this new element into the map
  for(int side_i = 0; side_i < VERTS_PER_ELEM; side_i++)
  {
    ElementAndFaceIdxType zs_pair =
      ElemNbrFinder(vertpair_to_elem_map, element_idx, side_i);
    SLIC_ASSERT(zs_pair.first == -1);
  }

  //Make a list of elements that shares at least 1 vertex of the new element
  std::set<IndexType> elem_list;
  for(int n = 0; n < VERTS_PER_ELEM; n++)
  {
    IndexArray ele_list_short = getElementsWithVertex(nlist[n]);
    for(unsigned int i = 0; i < ele_list_short.size(); i++)
    {
      elem_list.insert(ele_list_short[i]);
    }
  }

  //Check if any of the elements share a face with the new element.
  // If so, modify ee_rel to reflect that.
  for(std::set<IndexType>::iterator it = elem_list.begin(); it != elem_list.end();
      it++)
  {
    IndexType otherElementIdx = *it;
    if(otherElementIdx < 0 || otherElementIdx == element_idx) continue;
    for(IndexType s = 0; s < VERTS_PER_ELEM; s++)
    {
      IndexType otherSideIdx = s;
      // insert the pair
      ElementAndFaceIdxType zs_pair =
        ElemNbrFinder(vertpair_to_elem_map, otherElementIdx, otherSideIdx);

      //If zs_pair returned is the new element, save this nbr to set later
      IndexType foundElementIdx = zs_pair.first;
      IndexType foundSideIdx = zs_pair.second;

      if(foundElementIdx == element_idx)
      {
        // if there is already a neighbor on the save list, this mesh is not a
        // manifold.  Example: Having an edge with 3 faces...

        SLIC_ASSERT(ee_rel[otherElementIdx][otherSideIdx] < 0);

        ee_rel.modify(foundElementIdx, foundSideIdx, otherElementIdx);
        ee_rel.modify(otherElementIdx, otherSideIdx, foundElementIdx);

        //put new element pair back in queue to check if mesh is manifold
        ElemNbrFinder(vertpair_to_elem_map, foundElementIdx, foundSideIdx);
      }
    }
  }

  //update ve_rel
  for(int i = 0; i < VERTS_PER_ELEM; i++)
  {
    if(!ve_rel.isValidEntry(nlist[i]))
    {
      ve_rel.modify(nlist[i], 0, element_idx);
    }
  }

  return element_idx;
}

template <unsigned int TDIM, unsigned int SDIM, typename P>
void IAMesh<TDIM, SDIM, P>::fixVertexNeighborhood(
  IndexType vertex_idx,
  const std::vector<IndexType>& new_elements)
{
  using IndexPairType = std::pair<IndexType, IndexType>;
  using FaceVertMapType = std::map<IndexArray, IndexPairType>;
  using FaceVertPairType = std::pair<IndexArray, IndexPairType>;
  FaceVertMapType vert_map;

  const int SZ = new_elements.size();
  for(int i = 0; i < SZ; i++)
  {
    IndexType el = new_elements[i];

    for(int face_i = 0; face_i < VERTS_PER_ELEM; face_i++)
    {
      IndexArray fv_list = getElementFace(el, face_i);

      //only concerned with faces that contain the vertex in question
      if(!is_subset(vertex_idx, fv_list)) continue;

      //sort vertices on this face
      std::sort(fv_list.begin(), fv_list.end());

      std::pair<FaceVertMapType::iterator, bool> ret =
        vert_map.insert(FaceVertPairType(fv_list, IndexPairType(el, face_i)));

      if(!ret.second)  //found a matching face
      {
        IndexType nbr_elem = ret.first->second.first;
        IndexType nbr_face_i = ret.first->second.second;

        ee_rel.modify(el, face_i, nbr_elem);
        //ee_rel[el][face_i] = nbr_elem;

        ee_rel.modify(nbr_elem, nbr_face_i, el);
        //ee_rel[nbr_elem][nbr_face_i] = el;

        vert_map.erase(ret.first);
      }
    }
  }
}

/* Remove all the invalid entries in the IA structure*/
template <unsigned int TDIM, unsigned int SDIM, typename P>
void IAMesh<TDIM, SDIM, P>::compact()
{
  //Construct an array that maps original set indices to new compacted indices
  IndexArray vertex_set_map(vertex_set.size(), -1);
  IndexArray element_set_map(element_set.size(), -1);

  int v_count = 0;
  for(auto i = 0; i < vertex_set.size(); ++i)
  {
    if(vertex_set.isValidEntry(i))
    {
      vertex_set_map[i] = v_count++;
    }
  }

  int e_count = 0;
  for(auto i = 0; i < element_set.size(); ++i)
  {
    if(element_set.isValidEntry(i))
    {
      element_set_map[i] = e_count++;
    }
  }

  //update the relations
  typename ElementBoundaryRelation::RelationVec& ev_rel_data = ev_rel.data();
  for(auto i = 0; i < ev_rel.size(); ++i)
  {
    int new_entry_index = element_set_map[i];
    if(new_entry_index < 0) continue;
    for(auto j = 0; j < ev_rel[i].size(); ++j)
    {
      int val = ev_rel_data[i * ev_rel[i].size() + j];
      if(val != ElementBoundaryRelation::INVALID_INDEX)
      {
        val = vertex_set_map[val];
      }
      ev_rel_data[new_entry_index * ev_rel[i].size() + j] = val;
    }
  }
  ev_rel_data.resize(e_count * VERTS_PER_ELEM);

  typename ElementBoundaryRelation::RelationVec& ve_rel_data = ve_rel.data();
  for(int i = 0; i < (int)ve_rel.size(); i++)
  {
    int new_entry_index = vertex_set_map[i];
    if(new_entry_index < 0) continue;

    int val = ve_rel_data[i];
    if(val != VertexCoboundaryRelation::INVALID_INDEX)
    {
      val = element_set_map[val];
    }
    ve_rel_data[new_entry_index] = val;
  }
  ve_rel_data.resize(v_count);

  typename ElementBoundaryRelation::RelationVec& ee_rel_data = ee_rel.data();
  for(int i = 0; i < (int)ee_rel.size(); i++)
  {
    int new_entry_index = element_set_map[i];
    if(new_entry_index < 0) continue;
    for(int j = 0; j < (int)ee_rel[i].size(); j++)
    {
      int val = ee_rel_data[i * ee_rel[i].size() + j];
      if(val != ElementAdjacencyRelation::INVALID_INDEX)
      {
        val = element_set_map[val];
      }
      ee_rel_data[new_entry_index * ee_rel[i].size() + j] = val;
    }
  }
  ee_rel_data.resize(e_count * VERTS_PER_ELEM);

  //Update the map
  for(int i = 0; i < (int)vertex_set_map.size(); i++)
  {
    int new_entry_index = vertex_set_map[i];
    if(new_entry_index < 0) continue;
    vcoord_map[new_entry_index] = vcoord_map[i];
  }
  vcoord_map.resize(v_count);

  //update the sets
  vertex_set = VertexSet(v_count);
  element_set = ElementSet(e_count);
}

template <unsigned int TDIM, unsigned int SDIM, typename P>
bool IAMesh<TDIM, SDIM, P>::isEmpty() const
{
  return vertex_set.size() == 0;
}

template <unsigned int TDIM, unsigned int SDIM, typename P>
bool IAMesh<TDIM, SDIM, P>::isManifold(bool verboseOutput) const
{
  /* TODO check that...
   * - Each valid non-boundary entry in the ee_rel should have a counterpart
   *   that points to it.
   * - [More expensive] That we can reconstruct the star of every vertex and
   *   that that the star is either a ball or a half-ball
   * (alternatively, that the link is a sphere or a disk)
   * - Other things that ensures the mesh is manifold?
   */

  std::stringstream errSstr;

  bool bValid = true;

  if(!isValid(verboseOutput)) return false;

  // Each valid vertex should have a valid entry in ve_rel.
  for(IndexType i = 0; i < vertex_set.size(); ++i)
  {
    if(vertex_set.isValidEntry(i) && !ve_rel.isValidEntry(i))
    {
      if(verboseOutput)
      {
        errSstr << "\n\t vertex " << i
                << " is not connected to any elements.\n\t";
      }
      bValid = false;
    }
  }

  if(verboseOutput && !bValid)
  {
    SLIC_DEBUG(errSstr.str());
  }

  return bValid;
}

template <unsigned int TDIM, unsigned int SDIM, typename P>
bool IAMesh<TDIM, SDIM, P>::isValid(bool verboseOutput) const
{
  std::stringstream errSstr;

  bool bValid = true;

  bValid &= vertex_set.isValid(verboseOutput);
  bValid &= element_set.isValid(verboseOutput);
  bValid &= ev_rel.isValid(verboseOutput);
  bValid &= ve_rel.isValid(verboseOutput);
  bValid &= ee_rel.isValid(verboseOutput);
  bValid &= vcoord_map.isValid(verboseOutput);

  //Check that sizes for vertices match
  if(vertex_set.size() != ve_rel.size() || vertex_set.size() != vcoord_map.size())
  {
    if(verboseOutput)
    {
      errSstr << "\n\t vertex set and relation size don't match.\n\t";
      errSstr << "vertex size: " << vertex_set.size() << "\n\t"
              << "ve_rel size: " << ve_rel.size() << "\n\t"
              << "vcoord size: " << vcoord_map.size();
    }
    bValid = false;
  }

  //Check that sizes for elements match
  if(element_set.size() != ev_rel.size() || element_set.size() != ee_rel.size())
  {
    if(verboseOutput)
    {
      errSstr << "\n\t element set and relation size don't match.";
      errSstr << "element_set size: " << element_set.size() << "\n\t"
              << "ev_rel size: " << ev_rel.size() << "\n\t"
              << "ee_rel size: " << ee_rel.size();
    }
    bValid = false;
  }

  //Check that all ev_rel are valid if the element_set is valid.
  for(IndexType pos = 0; pos < element_set.size(); ++pos)
  {
    if(element_set.isValidEntry(pos))
    {
      for(IndexType rpos = 0; rpos < ev_rel[pos].size(); rpos++)
      {
        if(ev_rel[pos][rpos] == ElementBoundaryRelation::INVALID_INDEX)
        {
          if(verboseOutput)
          {
            errSstr
              << "\n\t* Element->Vertex relation contains an invalid entry"
              << " for a valid element \n\t pos: " << pos << ", entry: " << rpos
              << ".";
          }
          bValid = false;
        }
      }
    }
  }

  //Check that valid entries in relation/map map to valid entries in set
  for(IndexType pos = 0; pos < vertex_set.size(); ++pos)
  {
    if(ve_rel.isValidEntry(pos) && !vertex_set.isValidEntry(pos))
    {
      if(verboseOutput)
      {
        errSstr
          << "\n\t * Relation contains a valid entry with an invalid set entry"
          << " at pos " << pos << ".";
      }
      bValid = false;
    }
  }

  if(verboseOutput && !bValid)
  {
    SLIC_DEBUG(errSstr.str());
  }

  return bValid;
}

}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_IA_IMPL_H_
