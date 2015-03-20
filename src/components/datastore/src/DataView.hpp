/**
 * DataView.hpp
 *
 */

#ifndef DATAVIEW_HPP_
#define DATAVIEW_HPP_

#include <map>
#include <set>
#include <vector>

#include "conduit/conduit.h"

using conduit::Node;
using conduit::Schema;
using conduit::DataType;

namespace DataStoreNS
{

class DataGroup;
class DataStore;
class DataBuffer;
/**
 * \class
 *
 * \brief A class to manages interface to an actual piece of data, whether it be owned and allocated by the
 * or passed in from the outside.
 *
 * Requirements for this class are:
 *    - contain pointer to the owning entity of this .
 *    - contain collection of attributes
 *    - description of shape of data (if applicable)
 *    - description of type of data (if applicable)
 *
 */
class DataView
{

private:

  std::string m_name;

  DataGroup* m_parentGroup;

  /// pointer to the DataBuffer
  DataBuffer* m_dataBuffer;

  Node   m_node;
  Schema m_schema;
  Node   m_desc;

  // Cyrus's Note: we may still need m_viewStart, but it seems like keeping 
  // the buffer pointer could be better

public:


  DataView( const std::string& name,
            DataGroup* const parentGroup,
            DataBuffer* const dataBuffer );

  DataView( const std::string& name,
            DataGroup* const parentGroup,
            DataStore* const dataStore );


  /// copy constructor
  DataView(const DataView& source );


  /// destructor
  ~DataView();


  DataView* Init(const Schema &schema);
  DataView* Init(const DataType &dtype);
  
  DataView* ApplyDescriptor();

  void ReconcileWithBuffer();
  
  std::string GetName() const
  {return m_name;}

  DataView* SetDescriptor(const Schema &schema);  
  
  DataView* SetDescriptor(const DataType &dtype);
  
  /**
   *
   * @return m_schema
   */
  const Schema &GetDescriptor() const
  { return m_schema; }
  
  
  /// TODO: dangerous const issue needs to be resolved ??
  Node &GetNode()
  { return m_node; }  
  

  DataBuffer *GetBuffer()
  { return m_dataBuffer; }  
  
  
  
   void Print() const;  
   void Print(Node &n) const;
  ///@}


};




} /* namespace DataStore */
#endif /* DATAVIEW_HPP_ */
