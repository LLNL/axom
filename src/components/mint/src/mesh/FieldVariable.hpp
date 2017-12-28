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

#include "mint/Field.hpp"               // Base class

// axom includes
#include "axom/Macros.hpp"
#include "axom/Types.hpp"
#include "mint/Array.hpp"
#include "mint/config.hpp"
#include "mint/FieldTypes.hpp"
#include "slic/slic.hpp"

// C/C++ includes
#include <cstddef>  // for NULL
#include <string>   // for C++ string

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
   *
   * \param [in] name the name of this FieldVariable.
   * \param [in] size the number of tuples.
   * \param [in] capacity the number of tuples to allocate.
   * \param [in] num_components the number of components per tuple.
   */
  FieldVariable( const std::string& name,
                 IndexType num_tuples,
                 int num_components=1 ) :
    Field( name ),
    m_data( num_tuples, num_components )
  {
    SLIC_ASSERT( num_tuples >= 0 );
    SLIC_ASSERT( num_components >= 0 );

    this->m_type = field_of< FieldType >::type;
  }

  /*!
   * \brief Destructor.
   */
  virtual ~FieldVariable() { }

  virtual IndexType size() const { return m_data.size(); }

  virtual int getNumComponents() const { return m_data.getNumComponents(); }

  virtual IndexType getCapacity() const { return m_data.getCapacity(); }

  virtual double getResizeRatio() const { return m_data.getResizeRatio(); }

  /*!
   * \brief Returns a double pointer to the field data.
   * \return ptr pointer to the field data of type double.
   * \post ptr==AXOM_NULLPTR iff the data is not of type double.
   */
  virtual double* getDoublePtr() { return AXOM_NULLPTR; }

  /*!
   * \brief Returns a constant double pointer to the field data.
   * \return ptr constant pointer to the field data of type double.
   * \post ptr==AXOM_NULLPTR iff the data is not of type double.
   */
  virtual const double* getDoublePtr() const { return AXOM_NULLPTR; }

  /*!
   * \brief Returns an integer pointer to the field data.
   * \return ptr pointer to the field data of type int.
   * \post ptr==AXOM_NULLPTR iff the data is not of type int.
   */
  virtual int* getIntPtr() { return AXOM_NULLPTR; }

  /*!
   * \brief Returns a constant integer pointer to the field data.
   * \return ptr constant pointer to the field data of type int.
   * \post ptr==AXOM_NULLPTR iff the data is not of type int.
   */
  virtual const int* getIntPtr() const { return AXOM_NULLPTR; }

  virtual void resize( IndexType num_tuples )
  { m_data.resize( num_tuples ); }

  virtual void reserve( IndexType capacity )
  { m_data.reserve( capacity ); }

  virtual void setResizeRatio( double ratio )
  { m_data.setResizeRatio( ratio ); }

private:

  /*!
   * \brief FieldVariable constructor. Does nothing. Made private to prevent
   *  its use in application code.
   */
  FieldVariable() { }

  Array< FieldType > m_data;

  DISABLE_COPY_AND_ASSIGNMENT(FieldVariable);
  DISABLE_MOVE_AND_ASSIGNMENT(FieldVariable);
};

//------------------------------------------------------------------------------
//                  FIELD VARIABLE IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <>
double* FieldVariable< double >::getDoublePtr();

//------------------------------------------------------------------------------
template <>
const double* FieldVariable< double >::getDoublePtr() const;

//------------------------------------------------------------------------------
template <>
int* FieldVariable< int >::getIntPtr();

//------------------------------------------------------------------------------
template <>
const int* FieldVariable< int >::getIntPtr() const;



} /* namespace mint */
} /* namespace axom */

#endif /* FIELDVARIABLE_HPP_ */
