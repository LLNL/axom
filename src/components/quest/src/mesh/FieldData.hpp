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
 * \file FieldData.hpp
 *
 * \date Sep 19, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#ifndef FIELDDATA_HPP_
#define FIELDDATA_HPP_

#include "common/ATKMacros.hpp"

// C/C++ includes
#include <map>
#include <string>
#include <vector>

namespace meshtk {

// Forward Declarations
class Field;

class FieldData
{
public:

   /*!
    ****************************************************************************
    * \brief Default constructor. Creates an empty FieldData instance.
    ****************************************************************************
    */
   FieldData();

   /*!
    ****************************************************************************
    * \brief Destructor.
    ****************************************************************************
    */
   virtual ~FieldData();

   /*!
    ****************************************************************************
    * \brief Checks if the field with the given name exists.
    * \param [in] name the name of the field to check.
    * \return status true if the field exists, else, false.
    ****************************************************************************
    */
   bool hasField( const std::string& name );

   /*!
    ****************************************************************************
    * \brief Adds a field to this FieldData instance.
    * \param [in] f pointer to the field
    * \pre f != NULL
    * \pre this->hasField( f->getName() ) == false.
    * \note When a field is added ownership is transfered to the corresponding
    *  FieldData instance.
    ****************************************************************************
    */
   void addField( Field* f );

   /*!
    ****************************************************************************
    * \brief Returns the number of fields of this FieldData instance.
    * \return N the number of fiels in this instance.
    * \post N==0 \iff this->empty()==true.
    ****************************************************************************
    */
   int getNumberOfFields() const;

   /*!
    ****************************************************************************
    * \brief Returns the ith field of this FieldData instance.
    * \param [in] i the index of the field in query.
    * \return f pointer to the field in query.
    * \pre i >= 0 && i < this->getNumberOfFields()
    * \post f==ATK_NULLPTR \iff i < 0 || i >= this->getNumberOfFieds()
    ****************************************************************************
    */
   Field* getField( int i );

   /*!
    ***************************************************************************
    * \brief Returns the field with the given name.
    * \param [in] name the name of the field in query.
    * \return f pointer to the field in query.
    * \pre this->hasField( name )==true.
    * \post f==ATK_NULLPTR \iff this->hasField( name )==false.
    ***************************************************************************
    */
   Field* getField( const std::string& name );

   /*!
    ****************************************************************************
    * \brief Deletes all fields associated with this FieldData instance.
    * \post this->empty() == true.
    ****************************************************************************
    */
   void clear();

   /*!
    ****************************************************************************
    * \brief Checks if this FieldData instance is empty.
    * \return status true if empty, else, false.
    ****************************************************************************
    */
   bool empty() const;

private:

   // TODO: Revise this. We also need the ability to remove fields (?)
   std::vector< std::string > m_fields;
   std::map< std::string, Field* > m_container;

   DISABLE_COPY_AND_ASSIGNMENT(FieldData);
};

} /* namespace meshtk */

#endif /* FIELDDATA_HPP_ */
