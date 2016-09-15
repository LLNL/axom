/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file
 *
 * \date Sep 19, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#ifndef FIELD_HPP_
#define FIELD_HPP_

// ATK includes
#include "common/ATKMacros.hpp" // for DISABLE_COPY_AND_ASSIGNMENT

// C/C++ includes
#include <string>

namespace mint {

class Field
{
public:

  /*!
   *****************************************************************************
   * \brief Destructor.
   *****************************************************************************
   */
  virtual ~Field();

  /*!
   *****************************************************************************
   * \brief Returns the type of the field.
   * \return t the type of the field.
   * \post t < NUMBER_OF_FIELD_TYPES
   * \see FieldTypes
   *****************************************************************************
   */
  int getType() const { return m_type; };

  /*!
   *****************************************************************************
   * \brief Returns the name of the field.
   * \return name the name of the field.
   *****************************************************************************
   */
  std::string getName() const { return m_name; };

  /*!
   *****************************************************************************
   * \brief Returns the number of tuples in the field.
   * \return ntuples the number of tuples in the field.
   * \post ntuples >= 0.
   *****************************************************************************
   */
  int getNumTuples() const { return m_num_tuples; };

  /*!
   *****************************************************************************
   * \brief Retuns the number of components per tuple.
   * \return nc the number of components per tuple.
   * \post nc >= 1.
   *****************************************************************************
   */
  int getNumComponents() const { return m_num_components; };

  //TODO: Need to re-think API here!!!
  /*!
   *****************************************************************************
   * \brief Returns a double pointer to the field data.
   * \return ptr pointer to the field data.
   * \post ptr==ATK_NULLPTR iff the data is not of type double.
   *****************************************************************************
   */
  virtual double* getDoublePtr();

  /*!
   *****************************************************************************
   * \brief Returns an int pointer to the field data.
   * \return ptr pointer to the field data.
   * \post ptr==ATK_NULLPTR iff the is not an integer type.
   *****************************************************************************
   */
  virtual int* getIntPtr();

protected:

  /*!
   *****************************************************************************
   * \brief Default Constructor. Does nothing. Made it protected to prevent
   *  its use elsewhere.
   *****************************************************************************
   */
  Field();

  /*!
   *****************************************************************************
   * \brief Custom constructor. Creates a Field instance of the given size and
   *  number of components.
   * \param [in] name the name of the field.
   * \param [in] size the size of the field.
   * \param [in] num_components the number of components per element.
   * \note The custom constructor is intended for use in the constructor of
   *  subclasses to properly initialize the properties associated with this
   *  Field instance.
   *****************************************************************************
   */
  Field(const std::string& name, int size, int num_components);

  std::string m_name;    /*!< the name of the field  */
  int m_num_tuples;      /*!< total number of tuples */
  int m_num_components;  /*!< the number of components per tuple */
  int m_type;            /*!< the field type */

private:
  DISABLE_COPY_AND_ASSIGNMENT(Field);
};

} /* namespace mint */

#endif /* FIELD_HPP_ */
