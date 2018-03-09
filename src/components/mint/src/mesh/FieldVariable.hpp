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

#ifndef FIELDVARIABLE_HPP_
#define FIELDVARIABLE_HPP_

#include "mint/Field.hpp" // Base class

// axom includes
#include "axom/Macros.hpp"
#include "axom/Types.hpp"
#include "mint/FieldTypes.hpp"
#include "slic/slic.hpp"

// C/C++ includes
#include <cstddef> // for NULL
#include <string>  // for C++ string

namespace axom
{
namespace mint
{

template < typename FieldType >
class FieldVariable : public Field
{
public:

  /*!
   * \brief Creates a FieldVariable instance associated by the given name,
   *  size and number of components per tuple.
   * \param [in] name the name of this FieldVariable.
   * \param [in] size the number of tuples
   * \param [in] num_components the number of components per tuple.
   */
  FieldVariable(const std::string& name, int size, int num_components=1);

  /*!
   * \brief Destructor.
   */
  virtual ~FieldVariable();

  /*!
   * \brief Returns a double pointer to the field data.
   * \return ptr pointer to the field data of type double.
   * \post ptr==AXOM_NULLPTR iff the data is not of type double.
   */
  virtual double* getDoublePtr();

  /*!
   * \brief Returns a constant double pointer to the field data.
   * \return ptr constant pointer to the field data of type double.
   * \post ptr==AXOM_NULLPTR iff the data is not of type double.
   */
  virtual const double* getDoublePtr() const;

  /*!
   * \brief Returns an integer pointer to the field data.
   * \return ptr pointer to the field data of type int.
   * \post ptr==AXOM_NULLPTR iff the data is not of type int.
   */
  virtual int* getIntPtr();

  /*!
   * \brief Returns a constant integer pointer to the field data.
   * \return ptr constant pointer to the field data of type int.
   * \post ptr==AXOM_NULLPTR iff the data is not of type int.
   */
  virtual const int* getIntPtr() const;

private:

  FieldType* m_data;

  /*!
   * \brief FieldVariable constructor. Does nothing. Made private to prevent
   *  its use in application code.
   */
  FieldVariable() { m_data=AXOM_NULLPTR; };

  DISABLE_COPY_AND_ASSIGNMENT(FieldVariable);
  DISABLE_MOVE_AND_ASSIGNMENT(FieldVariable);
};

} /* namespace mint */
} /* namespace axom */

//------------------------------------------------------------------------------
//                  FIELD VARIABLE IMPLEMENTATION
//------------------------------------------------------------------------------
namespace axom
{
namespace mint
{

template < typename FieldType >
FieldVariable< FieldType >::FieldVariable(
  const std::string& name, int size, int nc ) :
  Field( name, size, nc )
{
  SLIC_ASSERT(  size >= 1 );
  SLIC_ASSERT(  nc >= 1 );

  m_type = field_of< FieldType >::type;
  m_data = new FieldType[ size * nc ];
}

//------------------------------------------------------------------------------
template < typename FieldType >
FieldVariable< FieldType >::~FieldVariable( )
{
  if ( m_data != AXOM_NULLPTR )
  {
    delete [] m_data;
    m_data = AXOM_NULLPTR;
  }
}

//------------------------------------------------------------------------------
template < typename FieldType >
double* FieldVariable< FieldType >::getDoublePtr()
{
  if ( m_type == DOUBLE_FIELD_TYPE )
  {
    return reinterpret_cast< double* >( m_data );
  }
  return AXOM_NULLPTR;
}

//------------------------------------------------------------------------------
template < typename FieldType >
const double* FieldVariable< FieldType >::getDoublePtr() const
{
  return const_cast< const double* >(
    const_cast< FieldVariable* >( this )->getDoublePtr() );
}

//------------------------------------------------------------------------------
template < typename FieldType >
int* FieldVariable< FieldType >::getIntPtr()
{
  if ( m_type == INTEGER_FIELD_TYPE )
  {
    return reinterpret_cast< int* >( m_data );
  }
  return AXOM_NULLPTR;
}

//------------------------------------------------------------------------------
template < typename FieldType >
const int* FieldVariable< FieldType >::getIntPtr() const
{
  return const_cast< const int* >(
    const_cast< FieldVariable* >( this )->getIntPtr() );
}

} /* namespace mint */
} /* namespace axom */

#endif /* FIELDVARIABLE_HPP_ */
