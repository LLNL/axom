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
  
  
    bool   HasBuffer() const
        { return m_buffer != nullptr;}

    bool   IsOpaque() const
        {return m_opaque;}
    
    void*  GetOpaque() const;
  
    /**
     *
     * @return m_schema
     */
    const Schema &GetDescriptor() const
    { return m_schema; }

    std::string GetName() const
    {return m_name;}
  
    /// note: in most cases, we want to use the const version of the node
    Node &GetNode()
    {return m_node; }  

    Node const& GetNode() const
    { return m_node; }

    /// for now, we assume a dataview always has a buffer and group
 
    DataBuffer *GetBuffer()
    { return m_buffer; }  
 
     DataBuffer const *GetBuffer() const
     { return m_buffer; }
 
    DataGroup* GetParent()
    {return m_group;}

    DataGroup const* GetParent() const
    {return m_group;}


    /// TODO: Bad name?
    bool Applied() const
    {return m_applied;}

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
              DataGroup* const parentGroup,
              void *opaque);



    /// copy constructor
    /// if this is public, we could get into some bookkeeping messes
    /// for example what do we do with the parent group pointer?
    DataView(const DataView& source );

    
    /// destructor is private
    ~DataView();


    /// this view's name
    std::string m_name;

    /// this views parent group
    DataGroup*  m_group;

    /// pointer to the DataBuffer
    DataBuffer* m_buffer;

    /// conduit schema used as descriptor
    Schema      m_schema;
    
    /// conduit node used to access the data
    Node        m_node;
    
    /// used to tell if the view is fully inited 
    /// may be a bad name
    bool        m_applied;
    
    /// bookkeeping for now, we should absorb this meta data into
    /// a more general descriptor
    bool        m_opaque;
};




} /* namespace DataStore */
#endif /* DATAVIEW_HPP_ */
