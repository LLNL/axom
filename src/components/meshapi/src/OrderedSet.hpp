/**
 * \file OrderedSet.h
 *
 * \brief Basic API for an ordered set of entities in a simulation
 *
 */

#ifndef MESHAPI_ORDERED_SET_H_
#define MESHAPI_ORDERED_SET_H_

#include <cstddef>
#include <vector>

#include <boost/iterator/counting_iterator.hpp>

#include "Utilities.hpp"

namespace asctoolkit{
namespace meshapi{


    class OrderedSet {
    public:
      typedef MeshIndexType     Index;
      typedef MeshSizeType      size_type;

      typedef boost::counting_iterator<Index> iterator;

    public:
      OrderedSet(size_type size): m_size(size) {}

      size_type size()  { return m_size; }

      iterator  begin()  { return iterator(0); }
      iterator  end()    { return iterator(m_size); }

      Index     operator[](Index idx) { return idx;}
      Index     at(Index idx)   ;


      void reset(size_type new_size) { throw NotImplementedException(); }

    private:
      int m_size;

    };


} // end namespace meshapi
} // end namespace asctoolkit

#endif //  MESHAPI_ORDERED_SET_H_
