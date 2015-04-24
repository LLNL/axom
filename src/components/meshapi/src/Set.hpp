/*
 * Set.h
 *
 *  Created on: Apr 23, 2015
 *      Author: weiss27
 */

#ifndef MESHAPI_SET_H_
#define MESHAPI_SET_H_

#include <cstddef>
#include <vector>

namespace asctoolkit{
namespace meshapi{

  class Set
  {
  public:
    typedef unsigned int             Index;
    typedef std::vector<Index>       ArrType;
    typedef ArrType::const_iterator  ArrCIter;
    typedef ArrType::iterator        ArrIter;
    typedef ArrType::size_type       size_type;


  public:
      Set () : m_parentSet(NULL) {}
      ~Set () {}


      Index operator[](size_type idx) { return m_entities[idx];}

      size_type size() const  { return m_entities.size(); }

      ArrIter  begin()        { return m_entities.begin(); }
      ArrCIter begin() const  { return m_entities.begin(); }

      ArrIter  end()          { return m_entities.end(); }
      ArrCIter end() const    { return m_entities.end(); }


      bool isSubset() const { return m_parentSet == NULL; }

  private:

      ArrType   m_entities;
      Set*      m_parentSet;
  };

} // end namespace meshapi
} // end namespace asctoolkit

#endif //  MESHAPI_SET_H_
