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
#include "Set.hpp"

namespace asctoolkit{
namespace meshapi{


    class OrderedSet  : public Set
    {
    public:
      typedef Set::SetIndex     SetIndex;
      typedef Set::size_type    size_type;

      typedef boost::counting_iterator<SetIndex> iterator;
      typedef std::pair<iterator,iterator> iterator_pair;

    public:
      OrderedSet(size_type size = size_type()): m_size(size) {}

      size_type size()  const { return m_size; }

      iterator  begin() const  { return iterator(0); }
      iterator  end()   const  { return iterator(m_size); }
      iterator_pair  range() const  { return std::make_pair(begin(), end()); }

      SetIndex     operator[](SetIndex idx) const { return idx;}
      SetIndex     at(SetIndex idx)   const ;


      void reset(size_type) { throw NotImplementedException(); }

      bool isValid(bool verboseOutput = false) const;
    private:
      SetIndex m_size;

    };


} // end namespace meshapi
} // end namespace asctoolkit

#endif //  MESHAPI_ORDERED_SET_H_
