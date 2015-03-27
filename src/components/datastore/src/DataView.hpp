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
public:
    friend class DataGroup;
    
    /// destructor is public
    /// a user detaches a view, they need some way to destroy it
    ~DataView();


    /// if there is a 1-1 relationship between this view and its buffer
    /// this will force a description & an allocate on the underlying
    /// buffer, otherwise an exception will be thrown
    DataView* Allocate();
    DataView* Allocate(const Schema &schema);
    DataView* Allocate(const DataType &dtype);
  
    /// sets the description that will be used to view data
    DataView* Declare(const Schema &schema);  
    DataView* Declare(const DataType &dtype);
  
    /// applies the description to buffer to init GetNode()
    DataView* Apply();
    DataView* Apply(const Schema &schema);  
    DataView* Apply(const DataType &dtype);
  
    /**
     *
     * @return m_schema
     */
    const Schema &GetDescriptor() const
    { return m_schema; }
  
  
    /// TODO: dangerous const issue needs to be resolved ??
    Node &GetNode()
    {return m_node; }  

    Node const& GetNode() const
    { return m_node; }
  

    DataBuffer *GetBuffer()
    { return m_dataBuffer; }  
  
    std::string GetName() const
    {return m_name;}
  
    DataGroup* GetParentGroup()
    {return m_parentGroup;}
  
  
    void Print() const;  
    void Info(Node &n) const;
    ///@}

private:
    // These are private b/c we expect Views to only be 
    // create from the context of a Group
    DataView( const std::string& name,
              DataGroup* const parentGroup,
              DataBuffer* const dataBuffer );

    DataView( const std::string& name,
              DataGroup* const parentGroup);


    /// copy constructor
    /// if this is public, we could get into some bookkeeping messes
    /// for example what do we do with the parent group pointer?
    DataView(const DataView& source );


    /// this view's name
    std::string m_name;

    /// this views parent group
    DataGroup*  m_parentGroup;

    /// pointer to the DataBuffer
    DataBuffer* m_dataBuffer;

    /// conduit schema used as descriptor
    Schema      m_schema;
    /// conduit node used to access the data
    Node        m_node;
};




} /* namespace DataStore */
#endif /* DATAVIEW_HPP_ */
