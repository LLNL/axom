/*
 * DatastoreInterface.hpp
 *
 *  Created on: Nov 17, 2014
 *      Author: it is everyones fault
 */

#ifndef DATASTOREINTERFACE_HPP_
#define DATASTOREINTERFACE_HPP_


#include <limits.h>
#include <stdint.h>
#include <vector>
#include <unordered_map>
#include <set>

#include "Types.hpp"
#include "Attributes.hpp"
#include "DataGroup.hpp"



namespace DataStore
{

/// placeholder for the datatypes
class DataType{};


/**
 * Statically-shaped arrays of data types
 *
 *  Single and multi-dimensional arrays without a mesh topology association (energy group bounds, LEOS table data).
 *  Scalar-valued mesh-based fields defined by a mesh centering (zone, node, face), such as density, temperature, etc.
 *  Scalar mesh fields with depth; i.e., multiple copies of related data for energy groups, quadrature points, time steps, etc.
 *  Vector-valued mesh fields (velocity, magnetic field), tensors (stress), etc.
 *  Vector-valued, tensor-valued, etc. mesh fields with depth.
 *
 *  By "statically-defined shape", we mean that the array dimensions are defined once and the arrays are not resized,
 *  except when significant aspects of the problem change (e.g., turn physics on or off). Typically, such arrays are
 *  described to the datastore and inserted once. They may be retrieved many times and their values may change, but the
 *  original data description remains valid and is sufficient for operations like serialization, etc.
 *  Typically, the choice of multi-dimensional array data layouts (e.g., multi-dimensional array of scalar field values
 *   or a single-dimensional array of vector field values) is based on access patterns used in a code for performance
 *   or physics algorithm reasons. Both data organization schemes need to be supported to accommodate all use cases.
 *   A general shape descriptor for array data is required for allocation and access (including slices). The following
 *   shape descriptor information was suggested
 */
class Dcoll1
{
public:
  Dcoll1();
  template< typename...ARGS > Dcoll1( ARGS...args );
  Dcoll1(const Dcoll1& source);
  Dcoll1(const Dcoll1&& source);
  Dcoll1& operator=(const Dcoll1& rhs);
  Dcoll1& operator=(const Dcoll1&& rhs);

  template< typename...ARGS >
  void reshape( ARGS...args );

  DataShape& getDescriptor();
  const DataShape& getDescriptor() const;



private:
  DataShape m_descriptor;
};




/**
 * This is just an interface. These functions should be going in various places in the data structure. I assume that
 * there will be something that oversees the datastore as a structure, and then calls functions at the "DataNode" level
 * to extract the data.
 */


  DataGroup* CreateDataStore( const std::string& name ) ;
  DataGroup* GetDataStore( const std::string& name );
  void DeleteDataStore( const std::string& name );



  /**
   *
   * @name Declare
   *
   */
  ///@{
  DataObject* CreateDataObject( DataGroup* const dataGroup, const std::string& name );

  ///@}





  /**
   *
   * @name Allocation
   *
   */
  ///@{

  inline DataObject* Allocate( DataObject* const object ) {}


  ///@}


  /**
   *
   * @name Data Insertion
   *
   */
  ///@{
  template< typename TYPE >
  void insert( const std::string& pathandName, TYPE& data )
  {

  }

  template< typename TYPE >
  void insert( const std::string& pathandName, TYPE* const data, const std::size_t length  )
  {

  }

  ///@}

  /**
   *
   * @name Add Attributes
   *
   */
  ///@{
  DataObject& SetAttribute( DataObject* const dataObject, const Attribute& att );
  DataObject& SetAttribute( const std::string& pathAndName, const Attribute& att );


  ///@}

  DataObject* SetDataShape( DataObject* const dataObject, const DataShape& shape );

  DataShape* GetDataShape( DataObject* const obj );

  /**
   *
   * @name QueryGroup
   *
   * \brief Group of functions to query data from the Datastore
   */
  ///@{

  DataGroup* GetDataGroup( const DataGroup* group, const Attribute& attribute );
  inline DataObject* GetDataObject( DataGroup* const dg, const std::string name ){}

  template< typename TYPE >
  typename std::enable_if<std::is_pointer<TYPE>::value,TYPE>::type GetData( DataObject* const obj )
  {
    return obj->GetData<TYPE>();
  }


  /**
   *
   * @param name
   * @return
   *  Dquery-1  Access object by name
   *
   *  Must Have
   */
  DataObject& DataQuery1( const std::string& pathAndName );
  DataObject& DataQuery1( const DataGroup& group, const std::string& name );

  /**
   *
   * @param index
   * @return
   *  Dquery-2  Access object by index  User may have a one time query on name that gives them the integer index of a
   *  data object. That index could be used for faster access rather than a name look up.
   *
   *  May have
   *
   *  This may require integer keys, with a look-up table for the name query?
   *
   *  Depending on the actual key type. If we are storing the objects in an array, then the native key type is an
   *  integer, and the name lookup (D-query-1) requires a table query to get the key, then this Dquery-2 simply assumes
   *  that the caller already has the valid integer key.
   *
   *
   */
  DataObject&  DataQuery2( const sizet& index );
  DataObject&  DataQuery2( const DataGroup& group, const sizet& index );

  /**
   *
   * @param att
   * @return
   *  Dquery-3  Access set of objects by attribute    Must Have
   *
   */
  DataGroup DataQuery3( const Attribute& att );

  /**
   *
   * @param type
   * @return
   *
   *  Dquery-4  Access set of objects by type
   *
   *  Must Have
   *
   */
  DataGroup DataQuery4( const DataType& type );

  /**
   *
   * @return
   *
   *  Dquery-5  Return the number of objects that  have a particular attribute. A code may allocate temporary buffer
   *  space based on the number of objects it needs to communicate of a given type.
   *
   *  Must Have
   *
   */
  int DataQuery5( const Attribute& att );
  int DataQuery5( const std::string& path, const Attribute& att );
  int DataQuery5( const DataGroup& group,  const Attribute& att );

  /**
   *
   * @return
   *
   *  Dquery-6  Answer whether an object with a given name has a particular attribute.  A code needs to ask whether
   *  certain objects are persistent, temporary, etc.
   *
   *  Must Have
   *
   */
  bool DataQuery6( const std::string& path, const Attribute& att );
  bool DataQuery6( const DataGroup& group, const std::string& name, const Attribute& att );


  /**
   *
   * @return
   *  Dquery-7  Support boolean expression queries that combine multiple attributes, attribute and type information, etc.
   *
   *  Some codes provide functionality to enable queries like: give me all zonal variables that are advected, where
   *  "zonal" is a predefined centering type and "advected" is a algorithm-specific string attribute convention applied
   *  in the code.
   *  Other codes provide functionality to get all "advected" fields that are also "intensive" using an expression
   *  syntax such as "advected & intensive".
   *
   *  Must Have
   *
   *  Depending on how this is done, it could involve one function call or several (get all zonals, then iterate over
   *  them to pick out the ones that have the advected attribute). The second option seems clumsy.
   *
   *  It may be a good idea to put reasonable constraints on the queries that are supported to avoid an overly complex
   *  implementation and prevent users from being to clever for their own good.
   *  We discussed whether & (and), | (or), and ~ (not) operators were sufficient. It may be that only & and | are
   *  really required. Providing ~ may coax users into doing things that can be done better (or more efficiently) by
   *  another means.
   *
   *  We agreed to ignore the "not" operator in the initial implementation and see if it arises as a serious use case.
   *
   *
   */
  DataGroup DataQuery7( );

  ///@}



}

#endif /* DATASTOREINTERFACE_HPP_ */
