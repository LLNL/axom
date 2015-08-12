/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#ifndef COMBINER_HPP
#define COMBINER_HPP

#include "lumberjack/MessageInfo.hpp"

namespace asctoolkit {
namespace lumberjack {

class Combiner {
    public:
        virtual const std::string id() = 0;
        virtual bool shouldMessageInfosBeCombined(const MessageInfo& leftMessage,
                                                  const MessageInfo& rightMessage) = 0;
        virtual void combine(MessageInfo& combined,
                             const MessageInfo& combinee, const int ranksLimit) = 0;
};

}
}

#endif