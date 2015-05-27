/**
 * \file IndirectionSet.h
 *
 * \brief Basic API for a set of entities in a simulation
 *
 */

#ifndef MESHAPI_INDIRECTION_SET_H_
#define MESHAPI_INDIRECTION_SET_H_

#include <cstddef>
#include <vector>

#include "Utilities.hpp"

namespace asctoolkit{
namespace meshapi{


  /**
   * \class IndirectionSet
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
   * Below is an initial implementation for a set with explicit indexes (encoded here using a vector).
   *
   * The interface is for constant access to the elements.
   */
  class IndirectionSet: public Set
  {
  public:
    typedef unsigned int                  Index;
    typedef std::vector<Index>            ArrType;

    typedef ArrType::const_iterator       ArrCIter;
    typedef std::pair<ArrCIter, ArrCIter> ArrCIterPair;

    typedef ArrType::iterator             ArrIter;
    typedef std::pair<ArrIter, ArrIter>   ArrIterPair;

    typedef ArrType::size_type            size_type;


  public:
      IndirectionSet () : m_parentSet(NULL) {}
      ~IndirectionSet () {}


      /**
       * \brief A function to initialize the indices of the set
       * \note Not yet implemented.
       * @param begin An iterator to the beginning of the desired range
       * @param end An iterator to the end of the desired range
       */
      template<typename IteratorType>
      void setupEntities(IteratorType begin, IteratorType end) { throw NotImplementedException(); }

      /**
       * \brief Unchecked random access to the entities of the set
       * @param idx The index of the desired element
       * @return A reference to the encoded index
       * \pre idx must be less than the number of elements in the set ( size() )
       */
      Index&       operator[](size_type idx)        { return m_entities[idx];}

      /**
       * \brief Unchecked random access to the entities of the set (const version)
       * @param idx The index of the desired element
       * @return A const reference to the encoded index
       * \pre idx must be less than the number of elements in the set ( size() )
       */
      Index const& operator[](size_type idx) const  { return m_entities[idx];}

      /**
       * \brief Checked random access to the entities of the set
       * \throws A std::out_of_range exception if idx >= the number of entities in the set
       * @param idx The index of the desired element
       * @return A reference to the encoded index
       * \pre idx must be less than the number of elements in the set ( size() )
       */
      Index&       at(size_type idx);

      /**
       * \brief Checked random access to the entities of the set (const version)
       * \throws A std::out_of_range exception if idx >= the number of entities in the set
       * @param idx The index of the desired element
       * @return A const reference to the encoded index
       * \pre idx must be less than the number of elements in the set ( size() )
       */
      Index const& at(size_type idx) const;

      /**
       * \brief Get the number of entities in the set
       * @return The number of entities in the set.
       */
      size_type size() const      { return m_entities.size(); }

      /**
       * @return An iterator to the beginning of the entities
       */
      ArrIter  begin()            { return m_entities.begin(); }

      /**
       * @return A const iterator to the beginning of the entities
       */
      ArrCIter begin() const      { return m_entities.begin(); }

      /**
       * \note This duplicates the functionality of begin() to ensure we get a const iterator
       * @return A const iterator to the beginning of the elements.
       */
      ArrCIter cbegin() const     { return m_entities.begin(); }

      /**
       * @return An iterator to the end of the entities.
       */
      ArrIter  end()              { return m_entities.end(); }

      /**
       * @return A const iterator to the end of the entities
       */
      ArrCIter end() const        { return m_entities.end(); }

      /**
       * \note This duplicates the functionality of end() to ensure we get a const iterator
       * @return A const iterator to the end of the entities
       */
      ArrCIter cend() const       { return m_entities.end(); }

      /**
        * @return A pair of begin/end iterators
       */
      ArrIterPair  range()        { return std::make_pair(begin(), end()); }

      /**
        * @return A pair of const begin/end iterators (const version)
       */
      ArrCIterPair range() const  { return std::make_pair(begin(), end()); }


      /**
       * \brief Determines if the Set is a Subset of another set.
       * @return true if the set is a subset of another set, otherwise false.
       */
      bool isSubset() const       { return m_parentSet == NULL; }

      /**
       * @return A pointer to the parent set.  NULL if there is no parent
       */
      Set* parentSet()            { return m_parentSet;}

  private:

      ArrType   m_entities;
      Set*      m_parentSet;
  };


  /**
   * \brief Two IndirectionSets are equal if they have the same cardinality and all indexes are the same (and in the same order)
   * \note Two sets of different types are (currently) considered to be unequal
   */
  inline bool operator==(IndirectionSet const& firstSet, IndirectionSet const& otherSet)
  {
      if(firstSet.size() != otherSet.size() )
          return false;
      for(SetIndex idx=0; idx< firstSet.size(); ++idx)
      {
          if( firstSet[idx] != otherSet[idx] )
              return false;
      }
      return true;
  }

  /**
   * \brief Two IndirectionSets are equal if they have the same cardinality and all indexes are the same (and in the same order)
   * \note Two sets of different types are (currently) considered to be unequal
   */
  inline bool operator!=(IndirectionSet const& firstSet, IndirectionSet const& otherSet) { return !(firstSet==otherSet); }




} // end namespace meshapi
} // end namespace asctoolkit

#endif //  MESHAPI_INDIRECTION_SET_H_
