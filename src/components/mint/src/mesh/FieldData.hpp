/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifndef FIELDDATA_HPP_
#define FIELDDATA_HPP_

#include "axom/Macros.hpp"
#include "mint/Field.hpp"
#include "mint/FieldVariable.hpp"
#include "mint/config.hpp"

// C/C++ includes
#include <map>
#include <string>
#include <vector>

namespace axom
{
namespace mint
{


class FieldData
{

public:

  /*!
   * \brief Default constructor. Creates an empty FieldData instance.
   */
  FieldData()
  {}

  /*!
   * \brief Destructor.
   */
  virtual ~FieldData()
  { clear(); }

  /*!
   * \brief Checks if the field with the given name exists.
   * \param [in] name the name of the field to check.
   * \return status true if the field exists, else, false.
   */
  inline bool hasField( const std::string& name ) const;

  /*!
   * \brief Returns the number of fields of this FieldData instance.
   * \return N the number of fiels in this instance.
   * \post N == 0 \iff this->empty() == true.
   */
  inline int getNumberOfFields() const;

  /*!
   * \brief Returns the ith field of this FieldData instance.
   * \param [in] i the index of the field in query.
   * \return f pointer to the field in query.
   * \pre i >= 0 && i < this->getNumberOfFields()
   * \post f == AXOM_NULLPTR \iff i < 0 || i >= this->getNumberOfFieds()
   */
  inline Field* getField( int i );

  /*!
   * \brief Returns the ith field of this FieldData instance as a constant
   * pointer.
   * \param [in] i the index of the field in query.
   * \return f constant pointer to the field in query.
   * \pre i >= 0 && i < this->getNumberOfFields()
   * \post f == AXOM_NULLPTR \iff i < 0 || i >= this->getNumberOfFieds()
   */
  inline const Field* getField( int i ) const;

  /*!
   * \brief Returns the field with the given name.
   * \param [in] name the name of the field in query.
   * \return f pointer to the field in query.
   * \pre this->hasField( name )==true.
   * \post f == AXOM_NULLPTR \iff this->hasField( name )==false.
   */
  inline Field* getField( const std::string& name );

  /*!
   * \brief Returns the field with the given name as a constant pointer.
   * \param [in] name the name of the field in query.
   * \return f constant pointer to the field in query.
   * \pre this->hasField( name )==true.
   * \post f == AXOM_NULLPTR \iff this->hasField( name )==false.
   */
  inline const Field* getField( const std::string& name ) const;

  /*!
   * \brief Deletes all fields associated with this FieldData instance.
   * \post this->empty() == true.
   */
  inline void clear();

  /*!
   * \brief Checks if this FieldData instance is empty.
   * \return status true if empty, else, false.
   */
  inline bool empty() const
  { return( m_container.empty() ); }

  inline void addField( Field* f )
  {
    SLIC_ASSERT(  f != AXOM_NULLPTR );
    SLIC_ASSERT(  this->hasField( f->getName() )==false );

    m_fields.push_back( f->getName() );
    m_container[ f->getName() ] = f;

    SLIC_ASSERT( m_fields.size() == m_container.size() );
  }

  template < typename FieldType >
  inline Field* addField( const std::string& name, IndexType size,
                          IndexType capacity, int num_components,
                          double resize_ratio );


  inline void resize( IndexType size );


  inline void reserve( IndexType capacity );


  inline void setResizeRatio( double ratio );

private:

  // TODO: Revise this. We also need the ability to remove fields (?)
  // Sidre has an class MapCollection where the vector holds the objects and
  // the map holds the indices. We should look into using it here. It supports
  // removal.
  std::vector< std::string > m_fields;
  std::map< std::string, Field* > m_container;

  DISABLE_COPY_AND_ASSIGNMENT(FieldData);
  DISABLE_MOVE_AND_ASSIGNMENT(FieldData);
};


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
inline bool FieldData::hasField( const std::string& name ) const
{
  bool status = false;

  if ( m_container.find( name ) != m_container.end() )
  {
    status = true;
  }

  return status;
}

//------------------------------------------------------------------------------
inline int FieldData::getNumberOfFields() const
{
  SLIC_ASSERT( m_fields.size() == m_container.size() );
  return static_cast< int >( m_fields.size() );
}

//------------------------------------------------------------------------------
inline Field* FieldData::getField( int i )
{
  SLIC_ASSERT( i >= 0 && i < this->getNumberOfFields() );

  if ( i < 0 || i >= this->getNumberOfFields() )
  {
    return AXOM_NULLPTR;
  }

  return this->getField( m_fields[ i ] );
}

//------------------------------------------------------------------------------
inline const Field* FieldData::getField( int i ) const
{
  SLIC_ASSERT( i >= 0 && i < this->getNumberOfFields() );

  if ( i < 0 || i >= this->getNumberOfFields() )
  {
    return AXOM_NULLPTR;
  }

  return this->getField( m_fields[ i ] );
}

//------------------------------------------------------------------------------
inline Field* FieldData::getField( const std::string& name )
{
  SLIC_ASSERT( this->hasField( name ) );
  return m_container.at( name );
}

//------------------------------------------------------------------------------
const Field* FieldData::getField( const std::string& name ) const
{
  SLIC_ASSERT( this->hasField( name ) );
  return m_container.at( name );
}

//------------------------------------------------------------------------------
inline void FieldData::clear()
{
  typename std::map< std::string, Field* >::iterator it;
  for ( it = m_container.begin() ; it != m_container.end() ; ++it )
  {
    delete it->second;
    it->second = AXOM_NULLPTR;
  }
  m_container.clear();
}

//------------------------------------------------------------------------------
template < typename FieldType >
inline Field* FieldData::addField( const std::string& name, IndexType size,
                                   IndexType capacity, int num_components,
                                   double resize_ratio )
{
  if ( hasField( name ) )
  {
    Field* f = getField( name );
    if ( f->getType() != field_of< FieldType >::type )
    {
      SLIC_WARNING( "Field with name " << name << " already exists but it " <<
                    "has a different type." );
    }

    return f;
  }

  Field* f = new FieldVariable< FieldType >( name, size, capacity,
                                             num_components, resize_ratio );
  m_fields.push_back( f->getName() );
  m_container[ f->getName() ] = f;

  SLIC_ASSERT( m_fields.size() == m_container.size() );
  return f;
}

//------------------------------------------------------------------------------
inline void FieldData::resize( IndexType size )
{
  typename std::map< std::string, Field* >::iterator it;
  for ( it = m_container.begin() ; it != m_container.end() ; ++it )
  {
    it->second->resize( size );
  }
}

//------------------------------------------------------------------------------
inline void FieldData::reserve( IndexType capacity )
{
  typename std::map< std::string, Field* >::iterator it;
  for ( it = m_container.begin() ; it != m_container.end() ; ++it )
  {
    it->second->reserve( capacity );
  }
}

//------------------------------------------------------------------------------
inline void FieldData::setResizeRatio( double ratio )
{
  typename std::map< std::string, Field* >::iterator it;
  for ( it = m_container.begin() ; it != m_container.end() ; ++it )
  {
    it->second->setResizeRatio( ratio );
  }
}


} /* namespace mint */
} /* namespace axom */

#endif /* FIELDDATA_HPP_ */
