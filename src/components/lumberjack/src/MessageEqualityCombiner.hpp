/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#ifndef MESSAGEEQUALITYCOMBINER_HPP
#define MESSAGEEQUALITYCOMBINER_HPP

#include "lumberjack/Combiner.hpp"
#include "lumberjack/MessageInfo.hpp"

#include <string>

namespace asctoolkit {
namespace lumberjack {

class MessageEqualityCombiner: public Combiner {
    public:
        MessageEqualityCombiner(): m_id("MessageEqualityCombiner") {}

        const std::string id()
        {
            return m_id;
        }

        bool shouldMessageInfosBeCombined(const MessageInfo& leftMessage,
                                          const MessageInfo& rightMessage)
        {
            if (leftMessage.message().compare(rightMessage.message()) == 0){
                return true;
            }
            return false;
        }

        void combine(MessageInfo& combined, const MessageInfo& combinee, const int ranksLimit)
        {
            combined.addRanks(combinee.ranks(), ranksLimit);
        }
    private:
        std::string m_id;
};

}
}

#endif