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

#ifndef FIELD_HPP_
#define FIELD_HPP_

// axom includes
#include "axom/Macros.hpp" // for DISABLE_COPY_AND_ASSIGNMENT
#include "axom/Types.hpp" // for AXOM_NULLPTR
#include "mint/FieldTypes.hpp"
#include "mint/DataTypes.hpp"

// C/C++ includes
#include <string>

namespace axom
{
namespace mint
{

class Field
{

public:

  /*!
   * \brief Destructor.
   */
  virtual ~Field()
  {}

  /*!
   * \brief Returns the type of the field.
   * \return t the type of the field.
   * \post t < NUMBER_OF_FIELD_TYPES
   * \see FieldTypes
   */
  int getType() const
  { return m_type; }

  /*!
   * \brief Returns the name of the field.
   * \return name the name of the field.
   */
  std::string getName() const
  { return m_name; }

  /*!
   * \brief Returns the number of tuples in the field.
   * \return ntuples the number of tuples in the field.
   * \post ntuples >= 0.
   */
  virtual localIndex size() const = 0;

  /*!
   * \brief Returns the number of components per tuple.
   * \return nc the number of components per tuple.
   * \post nc >= 1.
   */
  virtual int getNumComponents() const = 0;

  virtual localIndex getCapacity() const = 0;

  virtual double getResizeRatio() const = 0;

  virtual void resize( localIndex size ) = 0;

  virtual void reserve( localIndex capacity ) = 0;

  virtual void setResizeRatio( double ratio ) = 0;

  //TODO: Need to re-think API here!!!
  /*!
   * \brief Returns a double pointer to the field data.
   * \return ptr pointer to the field data.
   * \post ptr==AXOM_NULLPTR iff the data is not of type double.
   */
  virtual double* getDoublePtr() { return AXOM_NULLPTR; }

  /*!
   * \brief Returns a constant double pointer to the field data.
   * \return ptr constant pointer to the field data.
   * \post ptr==AXOM_NULLPTR iff the data is not of type double.
   */
  virtual const double* getDoublePtr() const { return AXOM_NULLPTR; }

  /*!
   * \brief Returns an int pointer to the field data.
   * \return ptr pointer to the field data.
   * \post ptr==AXOM_NULLPTR iff the is not an integer type.
   */
  virtual int* getIntPtr() { return AXOM_NULLPTR; }

  /*!
   * \brief Returns a constant int pointer to the field data.
   * \return ptr constant pointer to the field data.
   * \post ptr==AXOM_NULLPTR iff the is not an integer type.
   */
  virtual const int* getIntPtr() const { return AXOM_NULLPTR; }

protected:

  Field( const std::string& name ) :
    m_name( name ),
    m_type( UNDEFINED_FIELD_TYPE )
  {}

  std::string m_name;    /*!< the name of the field  */
  int m_type;            /*!< the field type */

private:

  DISABLE_COPY_AND_ASSIGNMENT(Field);
  DISABLE_MOVE_AND_ASSIGNMENT(Field);
};

} /* namespace mint */
} /* namespace axom */

#endif /* FIELD_HPP_ */
