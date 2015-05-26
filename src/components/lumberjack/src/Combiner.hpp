#ifndef COMBINER_HPP
#define COMBINER_HPP

#include "lumberjack/MessageInfo.hpp"

namespace asctoolkit {
namespace lumberjack {

class Combiner {
	public:
		virtual bool areMessageInfosEqual(const MessageInfo& leftMessage, const MessageInfo& rightMessage) = 0;
		virtual void combine(MessageInfo& combined, const MessageInfo& combinee, int ranksLimit) = 0;
};

}
}

#endif