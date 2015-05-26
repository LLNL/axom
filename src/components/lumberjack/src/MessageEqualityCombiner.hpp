#ifndef MESSAGEEQUALITYCOMBINER_HPP
#define MESSAGEEQUALITYCOMBINER_HPP

#include "lumberjack/Combiner.hpp"
#include "lumberjack/MessageInfo.hpp"

namespace asctoolkit {
namespace lumberjack {

class MessageEqualityCombiner: public Combiner {
	public:
		bool areMessageInfosEqual(const MessageInfo& leftMessage, const MessageInfo& rightMessage)
		{
			if (leftMessage.message().compare(rightMessage.message()) == 0){
				return true;
			}
			return false;
		}

		void combine(MessageInfo& combined, const MessageInfo& combinee, int ranksLimit)
		{
			combined.addRanks(combinee.ranks(), ranksLimit);
		}
};

}
}

#endif