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

#ifndef MINT_FIELDDATA_HPP_
#define MINT_FIELDDATA_HPP_

// Axom includes
#include "axom/Macros.hpp"        // for Axom macros

// Mint includes
#include "mint/config.hpp"           // for mint compile time definitions
#include "mint/Field.hpp"            // for mint::Field definition
#include "mint/FieldVariable.hpp"    // for mint::FieldVariable definition
#include "mint/FieldAssociation.hpp" // for mint::FieldAssociation enum

// C/C++ includes
#include <iterator>               // for std::advance
#include <string>                 // for std::string
#include <map>                    // for std::map
#include <vector>                 // for std::vector

namespace axom
{

// Axom Forward declarations
namespace sidre
{
class Group;
}

namespace mint
{

/*!
 * \class FieldData
 *
 * \brief Provides a container for storing fields associated with a specified
 *  mesh topology and methods to create, access and remove fields from the
 *  container.
 *
 *  A FieldData object may store fields associated with a single mesh topology,
 *  e.g., node-, cell- , face- or edge-centered. The FieldData object provides
 *  the means to create, access, modify and remove fields from the container.
 *
 *  Each field in the container is a FieldVariable instance which may
 *  store scalar, vector, tensor fields etc. as explained in more detail in
 *  the FieldVariable documentation.
 *
 *  When Mint is compiled with Sidre support, the FieldData object can be
 *  bound to a particular sidre::Group to store fields in a Sidre hierarchy
 *  according to the <a href="http://llnl-conduit.readthedocs.io/en/latest/">
 *  mesh blueprint </a>.
 *
 *  \note When using Sidre as the back-end storage, the action of creating or
 *   removing fields from the FieldData object is reflected in the associated
 *   Sidre hierarchy by default.
 *
 * \see Field
 * \see FieldVariable
 */
class FieldData
{
public:

/// \name Constructors
/// @{

  /*!
   * \brief Default constructor. Disabled.
   */
  FieldData( ) = delete;

  /*!
   * \brief Creates an empty FieldData instance with the specified mesh
   *  topology association.
   *
   * \param [in] association the mesh topology association, e.g., NODE_CENTERED
   *
   * \pre association >= NODE_CENTERED && association < NUM_FIELD_ASSOCIATIONS
   * \post this->getAssociation() == association
   * \post this->empty() == true
   * \post this->hasSidreGroup() == false
   * \post this->getNumFields() == 0
   */
  explicit FieldData( int association );

#ifdef MINT_USE_SIDRE
  /*!
   * \brief Constructs a FieldData instance which uses Sidre as the back-end
   *  data-store for storing fields and is bound to the specified group in the
   *  Sidre hierarchy.
   *
   * \param [in] association the mesh topology association, e.g., NODE_CENTERED
   * \param [in] field_group pointer to the Group containing fields.
   *
   * \note The user-supplied Sidre group should conform to the specifications of
   *  outlined in the <a href="http://llnl-conduit.readthedocs.io/en/latest/">
   *  mesh blueprint </a> for storing fields. The key requirements are:
   *  <ul>
   *   <li> Each field is a sub-group under the given group. </li>
   *   <li> Each sub-group consists of the following views:
   *        <ul>
   *          <li> <b> association </b>: encodes the mesh topology </li>
   *          <li> <b> volume_dependent </b>: volume scaling type </li>
   *          <li> <b> values </b>: holds the raw buffer of the field </li>
   *        </ul>
   *   </li>
   *  </ul>
   *
   * \note If the supplied Sidre group is not empty, the resulting FieldData
   *  instance will be populated with the fields that match the specified
   *  association.
   *
   * \pre association >= NODE_CENTERED && association < NUM_FIELD_ASSOCIATIONS
   * \pre field_group != AXOM_NULLPTR
   * \post this->getAssociation() == association
   * \post this->hasSidreGroup() == true
   */
  FieldData( int association, sidre::Group* field_group);
#endif

/// @}

  /*!
   * \brief Destructor.
   */
  ~FieldData() { clear(); }

/// \name Query Methods
/// @{

  /*!
   * \brief Returns the mesh topology association of fields in this container.
   * \return association the field association.
   * \see FieldAssociation.hpp
   */
  inline int getAssociation() const
  { return m_association; }

  /*!
   * \brief Checks if this FieldData instance is empty.
   * \return status true if empty, else, false.
   */
  inline bool empty() const
  { return ( m_name2idx.empty() ); }

  /*!
   * \brief Checks if the field with the given name exists.
   * \param [in] name the name of the field to check.
   * \return status true if the field exists, else, false.
   */
  inline bool hasField( const std::string& name ) const
  { return ( m_name2idx.find( name ) != m_name2idx.end() ); }

  /*!
   * \brief Returns the number of fields of this FieldData instance.
   * \return N the number of fiels in this instance.
   * \post N == 0 \iff this->empty() == true.
   */
  inline int getNumFields( ) const
  { return static_cast< int >( m_name2idx.size() ); }

  /*!
    * \brief Checks if a Sidre Group is associated with this FieldData instance.
    * \return status true if the FieldData is associated with Sidre, else, false.
    */
  inline bool hasSidreGroup( ) const;

/// @}

/// \name Methods to create new Fields
/// @{

  /*!
   * \brief Creates a new field with the given name consisting of the specified
   *  number of tuples and number of components.
   *
   * \param [in] name the user-supplied name to identify this field.
   * \param [in] num_tuples the number of tuples in the field.
   * \param [in] num_components number of components per tuple (optional)
   * \param [in] storeInSidre indicates whether to store the field in the
   *  associated Sidre group of this FieldData instance (optional).
   * \param [in] capacity initial max capacity for the field (optional).
   *
   * \tparam T the underlying data type of the field, e.g., double, int, etc.
   *
   * \note num_components is an optional argument. It defaults to '1' if a value
   *  is not explicitly specified by the caller.
   *
   * \note The 'storeInSidre' boolean argument that indicates whether the field
   *  data should be stored in Sidre. This is an optional argument which
   *  defaults to true and is only relevant when the code is compiled with Sidre
   *  support and the FieldData instance is associated with a Sidre Group. For
   *  persistent data, e.g., state variables etc., Sidre should own the data.
   *  However, for temporary variables it may be desirable to not modify the
   *  Sidre hierarchy.
   *
   * \return ptr pointer to the buffer associated with this field.
   *
   * \pre hasName( name ) == false
   * \pre num_components >= 1
   * \post ptr != AXOM_NULLPTR
   *
   */
  template < typename T >
  inline T* createField( const std::string& name,
                         IndexType num_tuples,
                         IndexType num_components=1,
                         bool storeInSidre=true,
                         IndexType capacity=Array< T >::USE_DEFAULT );

  /*!
   * \brief Creates a new field from a supplied external buffer which consists
   *  of the specified number of tuples and components.
   *
   * \param [in] name the user-supplied name of this field
   * \param [in] data supplied external buffer
   * \param [in] num_tuples the number of tuples in the field.
   * \param [in] num_components the numbere of components per tuple (optional).
   *
   * \tparam T the underlying data type of the field, e.g., double, int, etc.
   *
   * \note num_components is an optional argument. It defaults to '1' if a value
   *  is not explicitly specified by the caller.
   *
   * \note The supplied pointer must point to a buffer that is able to hold at
   *  least \f$ num\_tuples \times num\_components \f$
   *
   * \note An external field is not inserted in the Sidre tree hierarchy.
   *
   * \note Once an external field is created, subsequent calls to the FieldData
   *  methods to "resize()" and "reserve()" will fail.
   *
   * \return ptr pointer to the buffer associated with this field.
   *
   * \pre hasName( name ) == false
   * \pre data != AXOM_NULLPTR
   * \post ptr != AXOM_NULLPTR
   * \post ptr == data
   */
  template < typename T >
  inline T* createField( const std::string& name,
                         T* data,
                         IndexType num_tuples,
                         IndexType num_components=1 );
/// @}

/// \name Methods to remove Fields
/// @{

  /*!
   * \brief Removes the field with the given name.
   *
   * \param [in] name the name of the field to remove.
   *
   * \pre name.empty() == false
   * \pre hasField( name ) == true
   *
   */
  void removeField( const std::string& name );

  /*!
   * \brief Removes the ith field from the container.
   *
   * \param [in] i the index of the field to remove from the container
   *
   * \pre i >= 0 && i < getNumFields()
   */
  void removeField( int i );

/// @}

/// \name Field Object access methods
/// @{

  /*!
   * \brief Returns the ith field of this FieldData instance.
   *
   * \param [in] i the index of the field in query.
   * \return f pointer to the field in query.
   *
   * \pre i >= 0 && i < this->getNumFields()
   * \post f == AXOM_NULLPTR \iff i < 0 || i >= this->getNumberOfFieds()
   */
  /// @{

  inline Field* getField( int i )
  {
    auto it = m_name2idx.begin( );
    std::advance( it, i );
    return getField( it->first );
  }

  inline const Field* getField( int i ) const
  {
    auto it = m_name2idx.begin( );
    std::advance( it, i );
    return getField( it->first );
  }

  /// @}

  /*!
   * \brief Returns the field with the given name.
   *
   * \param [in] name the name of the field in query.
   * \return f pointer to the field in query.
   *
   * \pre this->hasField( name )==true.
   * \post f == AXOM_NULLPTR \iff this->hasField( name )==false.
   */
  /// @{

  inline Field* getField( const std::string& name )
  {
    SLIC_ASSERT( hasField( name ) );
    const int fldId = getFieldIndex( name );
    return (fldId != INVALID_FIELD_INDEX) ? getFieldAt(fldId) : AXOM_NULLPTR;
  }

  inline const Field* getField( const std::string& name ) const
  {
    SLIC_ASSERT( hasField( name ) );
    const int fldId = getFieldIndex( name );
    return (fldId != INVALID_FIELD_INDEX) ? getFieldAt(fldId) : AXOM_NULLPTR;
  }

  /// @}

/// @}

/// \name Field pointer access methods
/// @{

  /*!
   * \brief Returns pointer to the buffer of the field with the given name.
   *
   * \param [in] name the name of the field in query.
   * \param [out] num_tuples the number of tuples in the field (optional)
   * \param [out] num_components the number of components per tuple (optional).
   *
   * \return ptr pointer to the buffer of the specified field.
   *
   * \pre hasField( name ) == true
   * \post ptr != AXOM_NULLPTR
   */
  /// @{
  template < typename T >
  inline T* getFieldPtr( const std::string& name );

  template < typename T >
  inline T* getFieldPtr( const std::string& name, IndexType& num_tuples );

  template < typename T >
  inline T* getFieldPtr( const std::string& name,
                         IndexType& num_tuples,
                         IndexType& num_components );
  /// @}

  /*!
   * \brief Returns const pointer to the buffer of the field with the given name.
   *
   * \param [in] name the name of the field in query.
   * \param [out] num_tuples the number of tuples in the field (optional)
   * \param [out] num_components the number of components per tuple (optional)
   *
   * \return ptr pointer to the buffer of the specified field.
   *
   * \pre hasField( name ) == true.
   * \post ptr != AXOM_NULLPTR
   */
  /// @{
  template < typename T >
  inline const T* getFieldPtr( const std::string& name ) const;

  template < typename T >
  inline const T* getFieldPtr( const std::string& name,
                               IndexType& num_tuples ) const;

  template < typename T >
  inline const T* getFieldPtr( const std::string& name,
                               IndexType& num_tuples,
                               IndexType& num_components ) const;
  /// @}

/// @}

/// \name Size/Capacity Modifiers
/// @{

  /*!
   * \brief Changes the num-tuples of all fields in this FieldData instance.
   *
   * \param [in] newNumTuples the new number of tuples.
   *
   * \warning If the FieldData instance contains a FieldVariable that cannot be
   *  re-sized, e.g., it points to an external buffer, this method will abort
   *  with an error.
   *
   * \see FieldVariable
   */
  void resize( IndexType newNumTuples );

  /*!
   * \brief Changes the tuple capacity of all fields in this FieldData instance.
   *
   * \param [in] newCapacity the new max tuple capacity.
   *
   * \warning If the FieldData instance contains a FieldVariable that cannot be
   *  re-sized, e.g., it points to an external buffer, this method will abort
   *  with an error.
   *
   * \see FieldVariable
   */
  void reserve( IndexType newCapacity );

/// @}

private:
  static constexpr int INVALID_FIELD_INDEX = -1;

/// \name Private helper methods
/// @{

  /*!
   * \brief Deletes all fields associated with this FieldData instance.
   *
   * \note If the FieldData instance is using Sidre as the back-end data-store,
   *  this method does not remove the associated from the Sidre hierarchy.
   *
   * \post this->empty() == true.
   */
  void clear();

  /*!
   * \brief Returns the internal index of the field with the given name.
   *
   * \param [in] name the name of the field in query.
   * \return index the corresponding field index
   *
   * \warning INVALID_FIELD_INDEX will be returned if a field with the supplied
   *  name does not exit in this container instance.
   *
   * \pre hasField( name ) == true
   * \post index >= 0 && index < m_fields.size()
   * \post m_fields[ index ] != AXOM_NULLPTR
   */
  inline int getFieldIndex( const std::string& name ) const
  { return ( hasField( name ) ? m_name2idx.at(name) : INVALID_FIELD_INDEX ); }

  /*!
   * \brief Returns the Field object at the given index.
   *
   * \param [in] fieldIndex the index of the field in query.
   * \return f pointer to the field object at the given index.
   *
   * \pre fieldIndex >= 0 && fieldIndex < m_fields.size()
   * \post f != AXOM_NULLPTR
   */
  inline Field* getFieldAt( int fieldIndex ) const
  {
    // sanity checks
    SLIC_ASSERT( (fieldIndex >= 0) &&
                 (fieldIndex < static_cast< int >( m_fields.size() ) ) );
    SLIC_ASSERT( m_fields[ fieldIndex ] != AXOM_NULLPTR );
    return ( m_fields[ fieldIndex ] );
  }

  /*!
   * \brief Returns the index at which a new field can be stored.
   *
   * \return index the location at which a new field may be stored.
   * \post m_fields[ index ] == AXOM_NULLPTR
   */
  inline int getNewFieldIndex( );

  /*!
   * \brief Returns a string representation of the association name.
   *
   * \return name the string association name.
   * \post name.empty() == false
   */
  inline std::string getAssociationName( );

  /*!
   * \brief Removes the field at the given field index.
   *
   * \param [in] i the index of the field to remove.
   *
   * \pre i >= 0 && i < m_fields.size()
   * \pre m_fields[ i ]  != AXOM_NULLPTR
   */
  void removeFieldAt( int i );

/// @}

  int m_association;

  std::vector< Field* > m_fields;
  std::vector< int > m_nextidx;
  std::map< std::string, int > m_name2idx;

#ifdef MINT_USE_SIDRE
  sidre::Group* m_fields_group;
#endif

  DISABLE_COPY_AND_ASSIGNMENT(FieldData);
  DISABLE_MOVE_AND_ASSIGNMENT(FieldData);
};


//------------------------------------------------------------------------------
//  IMPLEMENTATION OF TEMPLATE & IN-LINE METHODS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
inline int FieldData::getNewFieldIndex( )
{
  int nextIdx = INVALID_FIELD_INDEX;

  if ( m_nextidx.empty() )
  {
    // no empty slots, append to the vector
    m_fields.push_back( AXOM_NULLPTR );
    nextIdx = m_fields.size( )-1;
  } // END if
  else
  {
    // use an empty slot
    nextIdx = m_nextidx.back();
    m_nextidx.pop_back();

    // ensure the object on that slot is NULL
    SLIC_ASSERT( m_fields[ nextIdx ] == AXOM_NULLPTR );
  } // END else

  return nextIdx;
}

//------------------------------------------------------------------------------
inline std::string FieldData::getAssociationName( )
{
  // TODO: note, currently edge/face data are not supported by the blueprint

  std::string name = "";
  switch ( m_association )
  {
  case NODE_CENTERED:
    name = "vertex";
    break;
  case CELL_CENTERED:
    name = "element";
    break;
  case EDGE_CENTERED:
    name = "edge";
    break;
  case FACE_CENTERED:
    name = "face";
    break;
  default:
    SLIC_ERROR( "undefined field association [" << m_association << "]" );
  } // END switch

  return ( name );
}

//------------------------------------------------------------------------------
inline bool FieldData::hasSidreGroup( ) const
{
#ifdef MINT_USE_SIDRE
  return ( m_fields_group != AXOM_NULLPTR );
#else
  return false;
#endif
}

//------------------------------------------------------------------------------
template < typename T >
inline T* FieldData::getFieldPtr( const std::string& name )
{
  IndexType num_tuples     = 0;
  IndexType num_components = 0;
  return ( getFieldPtr< T >( name, num_tuples, num_components) );
}

//------------------------------------------------------------------------------
template < typename T >
inline T* FieldData::getFieldPtr( const std::string& name,
                                  IndexType& num_tuples )
{
  IndexType num_components = 0;
  return ( getFieldPtr< T >( name, num_tuples, num_components) );
}

//------------------------------------------------------------------------------
template < typename T >
inline T* FieldData::getFieldPtr( const std::string& name,
                                  IndexType& num_tuples,
                                  IndexType& num_components )
{
  SLIC_ERROR_IF( !hasField(name), "field [" << name << "] does not exist!" );

  mint::Field* f = getField( name );
  SLIC_ASSERT( f != AXOM_NULLPTR );

  num_tuples     = f->getNumTuples();
  num_components = f->getNumComponents( );
  SLIC_ASSERT( num_components >= 1 );

  return mint::Field::getDataPtr< T >( f );
}

//------------------------------------------------------------------------------
template < typename T >
inline const T* FieldData::getFieldPtr( const std::string& name ) const
{
  IndexType num_tuples     = 0;
  IndexType num_components = 0;
  return ( getFieldPtr< T >( name, num_tuples, num_components ) );
}

//------------------------------------------------------------------------------
template < typename T >
inline const T* FieldData::getFieldPtr( const std::string& name,
                                        IndexType& num_tuples ) const
{
  IndexType num_components = 0;
  return ( getFieldPtr< T >( name, num_tuples, num_components ) );
}

//------------------------------------------------------------------------------
template < typename T >
inline const T* FieldData::getFieldPtr( const std::string& name,
                                        IndexType& num_tuples,
                                        IndexType& num_components ) const
{
  SLIC_ERROR_IF( !hasField(name), "field [" << name << "] does not exist!" );

  const mint::Field* f = getField( name );
  SLIC_ASSERT( f != AXOM_NULLPTR );

  num_tuples     = f->getNumTuples();
  num_components = f->getNumComponents( );
  SLIC_ASSERT( num_components >= 1 );

  return mint::Field::getDataPtr< T >( f );
}

//------------------------------------------------------------------------------
template < typename T >
inline T* FieldData::createField( const std::string& name,
                                  IndexType num_tuples,
                                  IndexType num_components,
                                  bool storeInSidre,
                                  IndexType capacity )
{
  SLIC_ERROR_IF( hasField( name ), "Field [" << name << "] already exists!");

  int fldIdx = getNewFieldIndex( );
  SLIC_ASSERT( ( fldIdx >= 0 ) &&
               ( fldIdx < static_cast< int >( m_fields.size() ) ) );
  SLIC_ASSERT( m_fields[ fldIdx ] == AXOM_NULLPTR );

  m_name2idx[ name ] = fldIdx;

  if ( capacity == Array< T > ::USE_DEFAULT )
  {
    // by default, the total capacity is the total number of tuples specified.
    capacity = num_tuples;
  }

  mint::Field* newField = AXOM_NULLPTR;
#ifdef MINT_USE_SIDRE
  // create the field on sidre
  if ( hasSidreGroup() && storeInSidre )
  {
    SLIC_ERROR_IF( m_fields_group->hasGroup( name ),
             "Field [" << name << "] already exists in the Sidre tree!" );

    sidre::Group *field = m_fields_group->createGroup( name );
    field->createView( "association" )->setString( getAssociationName() );
    field->createView( "volume_dependent" )->setString( "true" );

    // TODO: how should we bind this to the topology?
    field->createView( "type" )->setString( "topo" );

    sidre::View* values = field->createView( "values" );
    newField = new mint::FieldVariable< T >( name, values, num_tuples,
                                             num_components, capacity );
  } // END if
  else
  {
    newField = new mint::FieldVariable< T >( name, num_tuples, num_components,
                                             capacity );
  } // END else
#else
  newField = new mint::FieldVariable< T >( name, num_tuples, num_components,
                                           capacity );
#endif

  SLIC_ASSERT( newField != AXOM_NULLPTR );
  m_fields[ fldIdx ] = newField;

  return ( mint::Field::getDataPtr< T >( m_fields[ fldIdx ] ) );
}

//------------------------------------------------------------------------------
template < typename T >
inline T* FieldData::createField( const std::string& name,
                                  T* data,
                                  IndexType num_tuples,
                                  IndexType num_components )
{
  SLIC_ERROR_IF( data==AXOM_NULLPTR, "supplied buffer is NULL" );
  SLIC_ERROR_IF( hasField( name ), "Field [" << name << "] already exists!" );

  int fldIdx = getNewFieldIndex( );
  SLIC_ASSERT( ( fldIdx >= 0 ) &&
               ( fldIdx < static_cast< int >( m_fields.size() ) ) );
  SLIC_ASSERT( m_fields[ fldIdx ] == AXOM_NULLPTR );

  m_name2idx[ name ] = fldIdx;
  m_fields[ fldIdx ] =
    new mint::FieldVariable< T >( name, data, num_tuples, num_components );

  return ( mint::Field::getDataPtr< T >( m_fields[ fldIdx ] ) );
}

} /* namespace mint */
} /* namespace axom */

#endif /* FIELDDATA_HPP_ */
