/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

#ifndef FIELDVARIABLE_HPP_
#define FIELDVARIABLE_HPP_

#include "mint/Field.hpp" // Base class

// ATK includes
#include "common/ATKMacros.hpp"
#include "common/CommonTypes.hpp"
#include "mint/FieldTypes.hpp"
#include "slic/slic.hpp"

// C/C++ includes
#include <cstddef> // for NULL
#include <string>  // for C++ string

namespace axom {
namespace mint {

template < typename FieldType >
class FieldVariable:public Field
{
public:

  /*!
   ****************************************************************************
   * \brief Creates a FieldVariable instance associated by the given name,
   *  size and number of components per tuple.
   * \param [in] name the name of this FieldVariable.
   * \param [in] size the number of tuples
   * \param [in] num_components the number of components per tuple.
   ****************************************************************************
   */
  FieldVariable(const std::string& name, int size, int num_components=1);

  /*!
   ****************************************************************************
   * \brief Destructor.
   ****************************************************************************
   */
  virtual ~FieldVariable();

  /*!
   ****************************************************************************
   * \brief Returns a double pointer to the field data.
   * \return ptr pointer to the field data of type double.
   * \post ptr==ATK_NULLPTR iff the data is not of type double.
   ****************************************************************************
   */
  virtual double* getDoublePtr();

  /*!
   ****************************************************************************
   * \brief Returns an integer pointer to the field data.
   * \return ptr pointer to the field data of type int.
   * \post ptr==ATK_NULLPTR iff the data is not of type int.
   ****************************************************************************
   */
  virtual int* getIntPtr();

private:

  FieldType* m_data;

  /*!
   ****************************************************************************
   * \brief FieldVariable constructor. Does nothing. Made private to prevent
   *  its use in application code.
   ****************************************************************************
   */
  FieldVariable() { m_data=ATK_NULLPTR; };

  DISABLE_COPY_AND_ASSIGNMENT(FieldVariable);
  DISABLE_MOVE_AND_ASSIGNMENT(FieldVariable);
};

} /* namespace mint */
} /* namespace axom */

//------------------------------------------------------------------------------
//                  FIELD VARIABLE IMPLEMENTATION
//------------------------------------------------------------------------------
namespace axom {
namespace mint {

template < typename FieldType >
FieldVariable< FieldType >::FieldVariable(
  const std::string& name, int size, int nc ):
  Field( name, size, nc )
{
  SLIC_ASSERT(  size >= 1 );
  SLIC_ASSERT(  nc >= 1 );

  m_type = field_of< FieldType >::type;
  m_data = new FieldType[ size*nc ];
}

//------------------------------------------------------------------------------
template < typename FieldType >
FieldVariable< FieldType >::~FieldVariable( )
{
  if ( m_data != ATK_NULLPTR ) {
    delete [] m_data;
    m_data = ATK_NULLPTR;
  }
}

//------------------------------------------------------------------------------
template < typename FieldType >
double* FieldVariable< FieldType >::getDoublePtr()
{
  if ( m_type == DOUBLE_FIELD_TYPE ) {
    return reinterpret_cast< double* >( m_data );
  }
  return ATK_NULLPTR;
}

//------------------------------------------------------------------------------
template < typename FieldType >
int* FieldVariable< FieldType >::getIntPtr()
{
  if ( m_type == INTEGER_FIELD_TYPE ) {
    return reinterpret_cast< int* >( m_data );
  }
  return ATK_NULLPTR;
}

} /* namespace mint */
} /* namespace axom */

#endif /* FIELDVARIABLE_HPP_ */
