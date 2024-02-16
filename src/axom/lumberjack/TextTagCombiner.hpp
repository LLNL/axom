// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file TextTagCombiner.hpp
 *
 * \brief This file contains the class implementation of the
 * TextTagCombiner.
 *******************************************************************************
 */

#ifndef TEXTTAGCOMBINER_HPP
#define TEXTTAGCOMBINER_HPP

#include "axom/lumberjack/Combiner.hpp"
#include "axom/lumberjack/Message.hpp"

#include <string>

namespace axom
{
namespace lumberjack
{
/*!
 *******************************************************************************
 * \class TextTagCombiner
 *
 * \brief Combines Message classes if their Message::text and Message::tag
 *        are equal.
 *
 *  This class can be added to Lumberjack's Lumberjack by calling
 *  Lumberjack::addCombiner with a
 *  TextTagCombiner instance as its parameter.
 *
 * \see Combiner Lumberjack
 *******************************************************************************
 */
class TextTagCombiner : public Combiner
{
public:
  TextTagCombiner() { }

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
   * triggered for combination if both Message::text and Message::tag are equal.
   *
   * \param [in] leftMessage One of the Messages to be compared.
   * \param [in] rightMessage One of the Messages to be compared.
   *****************************************************************************
   */
  bool shouldMessagesBeCombined(const Message& leftMessage,
                                const Message& rightMessage)
  {
    return (leftMessage.text().compare(rightMessage.text()) == 0 &&
            leftMessage.tag().compare(rightMessage.tag()) == 0);
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
   *
   * \pre shouldMessagesBeCombined(combined, combinee) must be true
   *****************************************************************************
   */
  void combine(Message& combined, const Message& combinee, const int ranksLimit)
  {
    combined.addRanks(combinee.ranks(), combinee.count(), ranksLimit);
  }

private:
  const std::string m_id = "TextTagCombiner";
};

}  // end namespace lumberjack
}  // end namespace axom

#endif
