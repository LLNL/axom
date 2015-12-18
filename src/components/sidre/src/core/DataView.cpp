/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

/*!
 ******************************************************************************
 *
 * \file
 *
 * \brief   Implementation file for DataView class.
 *
 ******************************************************************************
 */


// Associated header file
#include "DataView.hpp"

// Other toolkit project headers
#include "common/CommonTypes.hpp"
#include "slic/slic.hpp"

// SiDRe project headers
#include "DataBuffer.hpp"
#include "DataGroup.hpp"
#include "DataStore.hpp"


namespace asctoolkit
{
namespace sidre
{

/*
 *************************************************************************
 *
 * Allocate data for view, previously declared.
 *
 *************************************************************************
 */
DataView * DataView::allocate()
{
  SLIC_ASSERT_MSG( allocateIsValid(), 
                   "View state does not allow data allocation");

  if ( allocateIsValid() ) 
  {
    if ( m_data_buffer == ATK_NULLPTR ) 
    {
       m_data_buffer = m_owning_group->getDataStore()->createBuffer();
       m_data_buffer->attachView(this);
    }

    if ( m_data_buffer->getNumViews() == 1 )
    {
      TypeID type = static_cast<TypeID>(m_schema.dtype().id());
      SidreLength num_elems = m_schema.dtype().number_of_elements();
      m_data_buffer->allocate(type, num_elems);
      m_state = ALLOCATED;
      apply();  
    }
  } 

  return this;
}

/*
 *************************************************************************
 *
 * Allocate data for view with type and number of elements.
 *
 *************************************************************************
 */
DataView * DataView::allocate( TypeID type, SidreLength num_elems)
{
  SLIC_ASSERT_MSG( allocateIsValid(), 
                   "View state does not allow data allocation");
  SLIC_ASSERT_MSG(num_elems >= 0, "Must allocate number of elements >= 0");

  if ( allocateIsValid() && num_elems >= 0 )
  {
    declare(type, num_elems);
    allocate();
    m_state = ALLOCATED;
    apply();
  }
  return this;
}

/*
 *************************************************************************
 *
 * Allocate data for view described by a Conduit data type object.
 *
 *************************************************************************
 */
DataView * DataView::allocate(const DataType& dtype)
{
  SLIC_ASSERT_MSG( allocateIsValid(), 
                   "View state does not allow data allocation");

  if ( allocateIsValid() )
  {
    declare(dtype);
    allocate();
    m_state = ALLOCATED;
    apply();
  }
  return this;
}

/*
 *************************************************************************
 *
 * Allocate data for view described by a Conduit schema object.
 *
 *************************************************************************
 */
DataView * DataView::allocate(const Schema& schema)
{
  SLIC_ASSERT_MSG( allocateIsValid(), 
                   "View state does not allow data allocation");

  if ( allocateIsValid() )
  {
    declare(schema);
    allocate();
    m_state = ALLOCATED;
    apply();
  }
  return this;
}

/*
 *************************************************************************
 *
 * Reallocate data for view to given number of elements.
 *
 *************************************************************************
 */
DataView * DataView::reallocate(SidreLength num_elems)
{
  SLIC_ASSERT_MSG( allocateIsValid(), 
                   "View state does not allow data allocation");
  SLIC_ASSERT_MSG(num_elems >= 0, "Must re-allocate number of elements >= 0");

  if ( allocateIsValid() && num_elems >= 0 )
  {
    // preserve current type
    TypeID vtype = static_cast<TypeID>(m_schema.dtype().id());
    declare(vtype, num_elems);
    m_data_buffer->reallocate(num_elems);
    m_state = ALLOCATED;
    apply();
  }
  return this;
}

/*
 *************************************************************************
 *
 * Reallocate data for view using a Conduit data type object.
 *
 *************************************************************************
 */
DataView * DataView::reallocate(const DataType& dtype)
{
  SLIC_ASSERT_MSG( allocateIsValid(), 
                   "View state does not allow data allocation");

  if ( allocateIsValid() )
  {
    TypeID type = static_cast<TypeID>(dtype.id());
    TypeID view_type = static_cast<TypeID>(m_schema.dtype().id());
    SLIC_ASSERT_MSG( type == view_type,
		     "Attempting to reallocate with a different type");
    if (type == view_type)
    {
      declare(dtype);
      SidreLength num_elems = dtype.number_of_elements();
      m_data_buffer->reallocate(num_elems);
      m_state = ALLOCATED;
      apply();
    }
  }
  return this;
}

/*
 *************************************************************************
 *
 * Reallocate data for view using a Conduit schema object.
 *
 *************************************************************************
 */
DataView * DataView::reallocate(const Schema& schema)
{
  SLIC_ASSERT_MSG( allocateIsValid(), 
                   "View state does not allow data allocation");

  if ( allocateIsValid() )
  {
    TypeID type = static_cast<TypeID>(schema.dtype().id());
    TypeID view_type = static_cast<TypeID>(m_schema.dtype().id());
    SLIC_ASSERT_MSG( type == view_type,
		     "Attempting to reallocate with a different type");
    if (type == view_type)
    {
      declare(schema);
      SidreLength num_elems = schema.dtype().number_of_elements();
      m_data_buffer->reallocate(num_elems);
      m_state = ALLOCATED;
      apply();
    }
  }
  return this;
}

/*
 *************************************************************************
 *
 * Attach buffer to view.
 *
 *************************************************************************
 */
DataView * DataView::attachBuffer(DataBuffer * buff)
{
  SLIC_ASSERT_MSG( attachBufferIsValid(),
                  "View state does not allow attaching buffer");
  SLIC_ASSERT_MSG( m_data_buffer == ATK_NULLPTR,
                  "Cannot attach buffer to view that already has a buffer");
  SLIC_CHECK( buff != ATK_NULLPTR );

  if ( attachBufferIsValid() && buff != ATK_NULLPTR )
  {
    buff->attachView(this);
    m_data_buffer = buff;
    m_state = BUFFER_ATTACHED;
    m_is_applied = false;

    if ( m_schema.total_bytes() <= m_data_buffer->getTotalBytes() )
    {
      apply(); 
    }
  }
  return this;
}

/*
 *************************************************************************
 *
 * Apply a previously declared data description to data held in the buffer.
 *
 *************************************************************************
 */
DataView * DataView::apply()
{
  SLIC_ASSERT_MSG( applyIsValid(),
                  "View state does not allow apply operation");
  SLIC_CHECK_MSG( !m_schema.dtype().is_empty(),
                  "View has no data description, apply() is a no-op");

  if ( applyIsValid() )
  {
    if ( m_data_buffer == ATK_NULLPTR || m_schema.dtype().is_empty() ) 
    {
      m_is_applied = false;
    }
    else 
    {
      if (m_state == EXTERNAL) 
      { 
        TypeID type = static_cast<TypeID>(m_schema.dtype().id());
        SidreLength num_elems = m_schema.dtype().number_of_elements();
        m_data_buffer->declare(type, num_elems);
      }
        
      m_node.set_external(m_schema, m_data_buffer->getVoidPtr());
      m_is_applied = true;
    }
  }
  return this;
}

/*
 *************************************************************************
 *
 * Apply given # elems, offset, stride description to data view.
 *
 *************************************************************************
 */
DataView * DataView::apply(SidreLength num_elems,
                           SidreLength offset,
                           SidreLength stride)
{
  SLIC_ASSERT_MSG( applyIsValid(),
                  "View state does not allow apply operation");
  SLIC_ASSERT_MSG( m_state != EXTERNAL, 
                  "View state does not allow apply operation");
  SLIC_ASSERT_MSG( m_data_buffer != ATK_NULLPTR,
                  "View needs buffer to get type information");
  SLIC_ASSERT_MSG(num_elems >= 0, "Must give number of elements >= 0");
  SLIC_ASSERT_MSG(offset >= 0, "Must give offset >= 0");

  if ( applyIsValid() && m_state != EXTERNAL &&
       m_data_buffer != ATK_NULLPTR &&
       num_elems >= 0 && offset >= 0)
  {
    DataType dtype = conduit::DataType::default_dtype(m_data_buffer->getTypeID());
    dtype.set_number_of_elements(num_elems);
    dtype.set_offset(offset * dtype.element_bytes() );
    dtype.set_stride(stride * dtype.element_bytes() );

    declare(dtype);
    apply();
  }
  return this;
}

/*
 *************************************************************************
 *
 * Apply given type, # elems, offset, stride desscription to data view.
 *
 *************************************************************************
 */
DataView * DataView::apply(TypeID type, SidreLength num_elems,
                                        SidreLength offset,
                                        SidreLength stride)
{
  SLIC_ASSERT_MSG( applyIsValid(),
                  "View state does not allow apply operation");
  SLIC_ASSERT_MSG(num_elems >= 0, "Must give number of elements >= 0");
  SLIC_ASSERT_MSG(offset >= 0, "Must give offset >= 0");

  if ( applyIsValid() &&
       num_elems >= 0 && offset >= 0)
  {
    DataType dtype = conduit::DataType::default_dtype(type);

    size_t bytes_per_elem = dtype.element_bytes();

    dtype.set_number_of_elements(num_elems);
    dtype.set_offset(offset * bytes_per_elem);
    dtype.set_stride(stride * bytes_per_elem);

    declare(dtype);
    apply();
  }
  return this;
}

/*
 *************************************************************************
 *
 * Apply given type, number of dimensions and shape to data view.
 *
 *************************************************************************
 */
DataView * DataView::apply(TypeID type, int ndims, SidreLength * shape)
{
  SLIC_ASSERT_MSG( applyIsValid(),
                  "View state does not allow apply operation");
  SLIC_ASSERT_MSG(ndims >= 1, "Must give number of dimensions >= 0");
  SLIC_ASSERT_MSG(shape != ATK_NULLPTR, "Pointer to shape cannot be null");

  if ( applyIsValid() && ndims >= 0 && shape != ATK_NULLPTR)
  {
    if (m_shape != ATK_NULLPTR)
    {
	m_shape->resize(ndims);
    }
    else
    {
	m_shape = new std::vector<SidreLength>(ndims);
    }

    SidreLength num_elems = 1;
    for (int i=0; i < ndims; i++)
    {
      num_elems *= shape[i];
      (*m_shape)[i] = shape[i];
    }
    apply(type, num_elems );
  }
  return this;
}

/*
 *************************************************************************
 *
 * Apply a Consuit data type description to data view.
 *
 *************************************************************************
 */
DataView * DataView::apply(const DataType &dtype)
{
  SLIC_ASSERT_MSG( applyIsValid(),
                  "View state does not allow apply operation");

  if ( applyIsValid() )
  {
    declare(dtype);
    apply();
  }
  return this;
}

/*
 *************************************************************************
 *
 * Apply a Conduit Schema to data view.
 *
 *************************************************************************
 */
DataView * DataView::apply(const Schema& schema)
{
  SLIC_ASSERT_MSG( applyIsValid(),
                  "View state does not allow apply operation");
 
  if ( applyIsValid() )
  { 
    declare(schema);
    apply();
  }
  return this;
}

int DataView::getNumDimensions() const
{
  if (m_shape == ATK_NULLPTR)
  {
    return 1;
  }
  else 
  {
    return m_shape->size();
  }
}

int DataView::getShape(int ndims, SidreLength * shape) const
{
  if (m_shape == ATK_NULLPTR)
  {
    if (ndims > 0)
    {
      shape[0] = getNumElements();
      return 1;
    }
    else
    {
      return -1;
    }
  }
  else 
  {
      if (static_cast<unsigned>(ndims) < m_shape->size())
    {
      return -1;
    }
    else
    {
#if 0
      for(std::vector<SidreLength>::iterator it = v.begin(); it != v.end(); ++it)
      {
          *shape++ = it.
      }
#else
      for(std::vector<SidreLength>::size_type i = 0; i != m_shape->size(); i++)
      {
        shape[i] = (*m_shape)[i];
      }
#endif
    }
    return m_shape->size();
  }
}

/*
 *************************************************************************
 *
 * Print JSON description of data view to stdout.
 *
 *************************************************************************
 */
void DataView::print() const
{
  print(std::cout);
}

/*
 *************************************************************************
 *
 * Set data view to hold external data.
 *
 *************************************************************************
 */
DataView * DataView::setExternalDataPtr(void * external_ptr)
{
  SLIC_ASSERT_MSG( setExternalDataPtrIsValid(), 
                  "View state does not allow setting external data pointer");

  if ( setExternalDataPtrIsValid() )
  {

    if ( m_data_buffer == ATK_NULLPTR )
    {
      m_data_buffer = m_owning_group->getDataStore()->createBuffer();
      m_data_buffer->setExternalData(external_ptr);
      m_data_buffer->attachView(this);
    }

    // todo, conduit should provide a check for if uint64 is a
    // good enough type to rep void *
    m_node.set((conduit::uint64)external_ptr);

    // 
    // If view has a data description, apply it.
    //
    apply();

    m_state = EXTERNAL;
  }

  return this;
}

/*
 *************************************************************************
 *
 * Print JSON description of data view to an  ostream.
 *
 *************************************************************************
 */
void DataView::print(std::ostream& os) const
{
  Node n;
  info(n);
  n.to_json_stream(os);
}

/*
 *************************************************************************
 *
 * Copy data view description to given Conduit node.
 *
 *************************************************************************
 */
void DataView::info(Node &n) const
{
  n["name"] = m_name;
  n["schema"] = m_schema.to_json();
  n["node"] = m_node.to_json();
  n["state"] = getStateStringName(m_state);
  n["is_applied"] = m_is_applied;
}

/*
 *************************************************************************
 *
 * PRIVATE ctor for DataView not associated with any data. 
 *
 *************************************************************************
 */
DataView::DataView( const std::string& name,
                    DataGroup * const owning_group)
  :   m_name(name),
  m_owning_group(owning_group),
  m_data_buffer(ATK_NULLPTR),
  m_schema(),
  m_node(),
  m_shape(ATK_NULLPTR),
  m_state(EMPTY),
  m_is_applied(false)
{
}

/*
 *************************************************************************
 *
 * PRIVATE dtor.
 *
 *************************************************************************
 */
DataView::~DataView()
{
  if (m_data_buffer != ATK_NULLPTR)
  {
    m_data_buffer->detachView(this);
  }
  if (m_shape != ATK_NULLPTR)
  {
    delete m_shape;
  }
}

/*
 *************************************************************************
 *
 * PRIVATE method to declare data view with type and number of elements.
 *
 *************************************************************************
 */
DataView * DataView::declare(TypeID type, SidreLength num_elems)
{
  SLIC_ASSERT_MSG(num_elems >= 0, "Must give number of elements >= 0");

  if ( num_elems >= 0 ) 
  {
    DataType dtype = conduit::DataType::default_dtype(type);
    dtype.set_number_of_elements(num_elems);
    m_schema.set(dtype);

    if ( m_state != EXTERNAL ) 
    {
      m_state = DESCRIBED;
    }

    m_is_applied = false;
  }
  return this;
}

/*
 *************************************************************************
 *
 * PRIVATE method to declare data view with a Conduit data type object.
 *
 *************************************************************************
 */
DataView * DataView::declare(const DataType& dtype)
{
  m_schema.set(dtype);

  if ( m_state != EXTERNAL ) 
  {
    m_state = DESCRIBED;
  }

  m_is_applied = false;

  return this;
}

/*
 *************************************************************************
 *
 * PRIVATE method to declare data view with a Conduit schema object.
 *
 *************************************************************************
 */
DataView * DataView::declare(const Schema& schema)
{
  m_schema.set(schema);

  if ( m_state != EXTERNAL ) 
  {
    m_state = DESCRIBED;
  }

  m_is_applied = false;

  return this;
}

/*
 *************************************************************************
 *
 * PRIVATE method returns true if view can allocate data; else false.
 *
 *************************************************************************
 */
bool DataView::allocateIsValid() const
{
  bool alloc_is_valid = false;
  if ( m_state != EXTERNAL && m_state != SCALAR && m_state != STRING )
  {
    alloc_is_valid = ( m_data_buffer == ATK_NULLPTR ||
                       (!m_data_buffer->isExternal() &&
                         m_data_buffer->getNumViews() == 1 ) );
  }
  return alloc_is_valid;
}

/*
 *************************************************************************
 *
 * PRIVATE method returns true if attaching buffer to view is valid; 
 * else false.
 *
 *************************************************************************
 */
bool DataView::attachBufferIsValid() const
{
   return ( m_state == EMPTY || m_state == DESCRIBED );
}

/*
 *************************************************************************
 *
 * PRIVATE method returns true if setting external data pointer on view 
 * is valid; else false.
 *
 *************************************************************************
 */
bool DataView::setExternalDataPtrIsValid() const
{
   return ( m_state == EMPTY || m_state == DESCRIBED || m_state == EXTERNAL );
}

/*
 *************************************************************************
 *
 * PRIVATE method returns true if apply ia a valid operation on view;
 * else false.
 *
 *************************************************************************
 */
bool DataView::applyIsValid() const
{
   return ( m_state != SCALAR && m_state != STRING );
}

/*
 *************************************************************************
 *
 * PRIVATE method returns string name of given view state enum value.
 *
 *************************************************************************
 */
char const * DataView::getStateStringName(State state) const
{
   char const * ret_string = NULL;

   switch ( state ) 
   {
      case EMPTY : 
      {
         ret_string = "EMPTY";
         break;
      }

      case DESCRIBED : 
      {
         ret_string = "DESCRIBED";
         break;
      }

      case ALLOCATED : 
      {
         ret_string = "ALLOCATED";
         break;
      }

      case BUFFER_ATTACHED : 
      {
         ret_string = "BUFFER_ATTACHED";
         break;
      }

      case EXTERNAL : 
      {
         ret_string = "EXTERNAL";
         break;
      }

      case SCALAR : 
      {
         ret_string = "SCALAR";
         break;
      }

      case STRING : 
      {
         ret_string = "STRING";
         break;
      }

      default :
      {
         ret_string = "/0";
      }
   }

   return( ret_string );
}


} /* end namespace sidre */
} /* end namespace asctoolkit */
