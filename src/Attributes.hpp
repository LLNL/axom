/*
 * Attributes.hpp
 *
 *  Created on: Dec 1, 2014
 */

#ifndef ATTRIBUTES_HPP_
#define ATTRIBUTES_HPP_

#include<string>
#include "Types.hpp"


namespace DataStore
{

/**
 * just putting something in here cause I don't know what it should look like. Do we want different value types?
 */
class Attribute
{
public:
  Attribute() = delete;

  Attribute( const std::string& key ):
    m_key(key),
    m_value(0.0)
  {}

  Attribute( const std::string& key, const real64 value ):
    m_key(key),
    m_value(value)
  {}

  ~Attribute() {}

  const std::string& Name() const
  { return m_key; }


private:
  std::string m_key;
  real64 m_value;

};
}

#endif /* ATTRIBUTES_HPP_ */
