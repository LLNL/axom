/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

/*!
 *******************************************************************************
 * \file MessageEqualityCombiner.hpp
 * \author Chris White (white238@llnl.gov)
 *
 * \brief This file contains the class implementation of the MessageEqualityCombiner.
 *******************************************************************************
 */

#ifndef MESSAGEEQUALITYCOMBINER_HPP
#define MESSAGEEQUALITYCOMBINER_HPP

#include "lumberjack/Combiner.hpp"
#include "lumberjack/MessageInfo.hpp"

#include <string>

namespace asctoolkit {
namespace lumberjack {

/*!
 *******************************************************************************
 * \class MessageEqualityCombiner
 *
 * \brief Combines MessageInfo classes if their (<MessageInfo>"::")<message>"()" are equal.
 *
 *  This class instance is automatically added to lumberjack's Logger for you. If you want it
 *  removed call (<Logger>"::")<removeMessageCombiner>"()" with "MessageEqualityCombiner" as it's
 *  parameter.
 *
 * \see Combiner Logger
 *******************************************************************************
 */
class MessageEqualityCombiner: public Combiner {
    public:
        MessageEqualityCombiner(): m_id("MessageEqualityCombiner") {}

        /*!
         *****************************************************************************
         * \brief Returns the unique string identifier for this combiner. Used by Logger
         *  to differentiate between other combiners.
         *****************************************************************************
         */
        const std::string id()
        {
            return m_id;
        }

        /*!
         *****************************************************************************
         * \brief Function used by Logger to indicate whether two messages should be
         * combined.  They are not actually combined by this function. MessageInfos are 
         * triggered for combination if both (<MessageInfo>"::")<message>"()" are equal.
         *
         * \param [in] leftMessageInfo One of the MessageInfos to be compared.
         * \param [in] rightMessageInfo One of the MessageInfos to be compared.
         *****************************************************************************
         */
        bool shouldMessageInfosBeCombined(const MessageInfo& leftMessageInfo,
                                          const MessageInfo& rightMessageInfo)
        {
            if (leftMessageInfo.message().compare(rightMessageInfo.message()) == 0){
                return true;
            }
            return false;
        }

        /*!
         *****************************************************************************
         * \brief Combines the combinee into the combined MessageInfo.
         *
         * The only thing truly combined in this Combiner is the ranks from combinee to
         * combined, since message is already equal.
         *
         * \param [in,out] combined the MessageInfo that will be modified.
         * \param [in] combinee the MessageInfo that is combined into the other.
         * \param [in] ranksLimit The limit on how many individual ranks are tracked in
         * the combined MessageInfo. MessageInfo::rankCount is always incremented.
         *****************************************************************************
         */
        void combine(MessageInfo& combined, const MessageInfo& combinee, const int ranksLimit)
        {
            combined.addRanks(combinee.ranks(), ranksLimit);
        }
    private:
        std::string m_id;
};

} // end namespace lumberjack
} // end namespace asctoolkit

#endif
