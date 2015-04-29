/**
 * DataView.hpp
 *
 */

#ifndef DATAVIEW_HPP_
#define DATAVIEW_HPP_

#include <map>
#include <set>
#include <vector>
#include "Types.hpp"

#include "conduit/conduit.h"

using conduit::Node;
using conduit::Schema;
using conduit::DataType;

namespace sidre
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
    /// buffer, otherwise an assertion will result
    DataView* allocate();
    DataView* allocate(const Schema &schema);
    DataView* allocate(const DataType &dtype);
  
    /// sets the description that will be used to view data
    DataView* declare(const Schema &schema);  
    DataView* declare(const DataType &dtype);
  
    /// applies the description to buffer to init GetNode()
    DataView* apply();
    DataView* apply(const Schema &schema);  
    DataView* apply(const DataType &dtype);
  
  
    bool   hasBuffer() const
    { return m_buffer != nullptr;}

    bool   isOpaque() const
    {return m_opaque;}
    
    void*  getOpaque() const;
  
    /**
     *
     * @return m_schema
     */
    const Schema &getDescriptor() const
    { return m_schema; }

    std::string getName() const
    {return m_name;}
  
    /// note: in most cases, we want to use the const version of the node
    Node &getNode()
    {return m_node; }  

    Node const& getNode() const
    { return m_node; }

    /// for now, we assume a dataview always has a buffer and group
 
    DataBuffer *getBuffer()
    { return m_buffer; }  
 
     DataBuffer const *getBuffer() const
     { return m_buffer; }
 
    DataGroup* getParent()
    {return m_group;}

    DataGroup const* getParent() const
    {return m_group;}


    /// TODO: Bad name?
    bool isApplied() const
    {return m_applied;}

    void print() const;  
    void info(Node &n) const;
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




} /* namespace sidre */
#endif /* DATAVIEW_HPP_ */
