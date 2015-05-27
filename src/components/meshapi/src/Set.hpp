/**
 * \file Set.h
 *
 * \brief Basic API for a set of entities in a simulation
 *
 */

#ifndef MESHAPI_SET_H_
#define MESHAPI_SET_H_

#include <cstddef>
#include <vector>

#include "Utilities.hpp"

namespace asctoolkit{
namespace meshapi{


  /**
   * \class Set
   *
   * \brief An indexed set (a tuple) of entities in a simulation
   *
   * A container class for a set of entities in a simulation.
   * Each entity has an index.
   *
   * Examples of sets include:
   * <ol>
   *  <li> Mesh elements: vertices, edges, faces, cells
   *  <li> Subzonal elements: sides, corners, finite element degrees of freedom
   *  <li> Boundary elements: external surfaces and springs
   *  <li> Elements of a space partition: e.g. Domains in a block structured mesh, leaf nodes of an octree/kd-tree
   *  <li> AMR bricks / tiles
   *  <li> Thread ids, MPI ranks, warps, thread groups, etc...
   *  <li> particles
   *  <li> boundary conditions
   *  <li> ...
   * </ol>
   *
   * Examples of subsets include:
   * <ol>
   *  <li> Regions
   *  <li> Ghost cells -- send, receive
   *  <li> Boundary cells -- external surface
   * </ol>
   *
   * Note: Elements of a set do not necessarily need explicit indices.
   * E.g. if we have a contiguous range of elements (or slices of contiguous ranges),
   * they can be implicitly encoded.
   *
   * Thus, we can have
   * <ol>
   *  <li> Implicit indexes -- all we need here is a size operator
   *  <li> Sliced indices -- here we need the dimension and the striding
   *  <li> Explicit indices -- for a subset, we need the indices with respect to some other indexing scheme
   * </ol>
   *
   * The interface is for constant access to the elements.
   */
  class Set
  {
  public:
    typedef MeshIndexType                 SetIndex;
    typedef std::vector<SetIndex>            ArrType;

    typedef ArrType::const_iterator       ArrCIter;
    typedef std::pair<ArrCIter, ArrCIter> ArrCIterPair;

    typedef ArrType::iterator             ArrIter;
    typedef std::pair<ArrIter, ArrIter>   ArrIterPair;

    typedef ArrType::size_type            size_type;


  public:
      // Set () {}
      virtual ~Set () {}


      /**
       * \brief Unchecked random access to the entities of the set (const version)
       * @param idx The index of the desired element
       * @return A const reference to the encoded index
       * \pre idx must be less than the number of elements in the set ( size() )
       */
      virtual SetIndex operator[](SetIndex idx) const  =0;

      /**
       * \brief Get the number of entities in the set
       * @return The number of entities in the set.
       */
      virtual size_type size() const      =0;

#if 0
      /**
       * @return An iterator to the beginning of the entities
       */
      virtual ArrIter  begin()            =0;

      /**
       * @return A const iterator to the beginning of the entities
       */
      virtual ArrCIter begin() const      =0;

      /**
       * @return An iterator to the end of the entities.
       */
      virtual ArrIter  end()              =0;

      /**
       * @return A const iterator to the end of the entities
       */
      virtual ArrCIter end() const        =0;

      /**
        * @return A pair of begin/end iterators
       */
      virtual ArrIterPair  range()        =0;

      /**
        * @return A pair of const begin/end iterators (const version)
       */
      virtual ArrCIterPair range() const  =0;

      /**
       * \brief Determines if the Set is a Subset of another set.
       * @return true if the set is a subset of another set, otherwise false.
       */
      virtual bool isSubset() const       =0;

      /**
       * @return A pointer to the parent set.  NULL if there is no parent
       */
      virtual Set* parentSet()            =0;
#endif
  };

} // end namespace meshapi
} // end namespace asctoolkit

#endif //  MESHAPI_SET_H_
