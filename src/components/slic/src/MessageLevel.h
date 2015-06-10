/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

/*!
 *******************************************************************************
 * \file MessageLevel.h
 *
 * \date May 8, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *
 *******************************************************************************
 */

#ifndef MESSAGELEVEL_H_
#define MESSAGELEVEL_H_

#include <string>

namespace asctoolkit {

namespace slic {

namespace message {

/*!
 *******************************************************************************
 * \enum MessageType
 *
 * \brief Enumerates the different types of messaging.
 *
 * \note The ordering of the enumerators reflects the level of severity of
 *  a message.
 *
 * \see Logger
 *******************************************************************************
 */
enum Level {
  Fatal,     //!< FATAL log a non-recoverable event.
  Error,     //!< ERROR log an error that *may* be recoverable.
  Warning,   //!< WARNING log a warning.
  Info,      //!< INFO log information that is useful for users & developers.
  Debug,     //!< DEBUG log information that is useful for developers.

  Num_Levels //!< Num_Levels
};

/*!
 * \brief Array of strings corresponding to the Level enum.
 */
const static std::string MessageLevelName[ Num_Levels ] = {
    "FATAL",
    "ERROR",
    "WARNING",
    "INFO",
    "DEBUG"
};

/*!
 *******************************************************************************
 * \brief Returns the string name represetnation of the given level.
 * \param [in] l the level in query.
 * \return name a string corresponding to the name of the given leve.
 * \pre l >= 0 && l < Num_levels
 * \post name="UNKNOWN-LEVEL" \iff (l < 0) || (l >= Num_Levels)
 *******************************************************************************
 */
static inline std::string getLevelAsString( Level l )
{
  if ( l < 0 || l >= Num_Levels ) {

    return ( "UNKNOWN-LEVEL" );

  }

  return( MessageLevelName[ l ] );
}

} /* namespace loglevel */

} /* namespace slic */

} /* namespace asctoolkit */

#endif /* MESSAGELEVEL_H_ */
