/*!
 *******************************************************************************
 * \file MessageType.h
 *
 * \date May 8, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *
 *******************************************************************************
 */

#ifndef MESSAGETYPE_H_
#define MESSAGETYPE_H_

namespace asctoolkit {

namespace logapi {


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
enum MessageType {
  Fatal,        //!< FATAL log a non-recoverable event.
  Error,        //!< ERROR log an error that *may* be recoverable.
  Warning,      //!< WARNING log a warning.
  Info,         //!< INFO log information that is useful for users & developers.
  Debug,        //!< DEBUG log information that is useful for developers.

  Num_Msg_Types //!< Num_Msg_Types
};

} /* namespace logapi */

} /* namespace asctoolkit */

#endif /* MESSAGETYPE_H_ */
