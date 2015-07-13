#ifndef COMBINER_HPP
#define COMBINER_HPP

#include "lumberjack/MessageInfo.hpp"

namespace asctoolkit {
namespace lumberjack {

class Combiner {
    public:
        virtual bool shouldMessageInfosBeCombined(const MessageInfo& leftMessage,
                                                  const MessageInfo& rightMessage) = 0;
        virtual void combine(MessageInfo& combined,
                             const MessageInfo& combinee, const int ranksLimit) = 0;
};

}
}

#endif