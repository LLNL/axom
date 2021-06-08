// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file TextEqualityCombiner.hpp
 *
 * \brief This file contains the class implementation of the
 * TextEqualityCombiner.
 *******************************************************************************
 */

#ifndef TEXTEQUALITYCOMBINER_HPP
#define TEXTEQUALITYCOMBINER_HPP

#include "axom/lumberjack/Combiner.hpp"
#include "axom/lumberjack/Message.hpp"

#include <string>

namespace axom
{
namespace lumberjack
{
/*!
 *******************************************************************************
 * \class TextEqualityCombiner
 *
 * \brief Combines Message classes if their Message::text are equal.
 *
 *  This class instance is automatically added to Lumberjack's Lumberjack for
 *  you. If you want it removed call Lumberjack::removeCombiner with the string
 * "TextEqualityCombiner" as it's parameter.
 *
 * \see Combiner Lumberjack
 *******************************************************************************
 */
class TextEqualityCombiner : public Combiner
{
public:
  TextEqualityCombiner() : m_id("TextEqualityCombiner") { }

  /*!
   *****************************************************************************
   * \brief Returns the unique string identifier for this combiner. Used by
   *  Lumberjack to differentiate between other combiners.
   *****************************************************************************
   */
  const std::string id() { return m_id; }

  /*!
   *****************************************************************************
   * \brief Function used by Lumberjack to indicate whether two messages should
   *  be combined.
   *
   * They are not actually combined by this function. Message classes are
   * triggered for combination if both Message::text are equal.
   *
   * \param [in] leftMessage One of the Messages to be compared.
   * \param [in] rightMessage One of the Messages to be compared.
   *****************************************************************************
   */
  bool shouldMessagesBeCombined(const Message& leftMessage,
                                const Message& rightMessage)
  {
    if(leftMessage.text().compare(rightMessage.text()) == 0)
    {
      return true;
    }
    return false;
  }

  /*!
   *****************************************************************************
   * \brief Combines the combinee into the combined Message.
   *
   * The only thing truly combined in this Combiner is the ranks from combinee
   * to combined, since text is already equal.
   *
   * \param [in,out] combined the Message that will be modified.
   * \param [in] combinee the Message that is combined into the other.
   * \param [in] ranksLimit The limit on how many individual ranks are tracked
   *  in the combined Message. Message::rankCount is always incremented.
   *****************************************************************************
   */
  void combine(Message& combined, const Message& combinee, const int ranksLimit)
  {
    combined.addRanks(combinee.ranks(), combinee.count(), ranksLimit);
  }

private:
  std::string m_id;
};

}  // end namespace lumberjack
}  // end namespace axom

#endif
