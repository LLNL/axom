#include "lumberjack/Logger.hpp"

#include "common/CommonTypes.hpp"

namespace asctoolkit {
namespace lumberjack {

void Logger::initialize(Communicator* communicator, int ranksLimit)
{
    m_communicator = communicator;
    m_ranksLimit = ranksLimit;
    m_combiners.push_back(new MessageEqualityCombiner);
}

void Logger::finalize()
{
    m_communicator = ATK_NULLPTR;
    clearMessageCombiners();
}

void Logger::addMessageCombiner(Combiner* combiner)
{
    bool identifierFound = false;
    for (int i=0; i<(int)m_combiners.size(); ++i){
        if (m_combiners[i]->id() == combiner->id()){
            identifierFound = true;
            break;
        }
    }
    if (!identifierFound){
        m_combiners.push_back(combiner);
    }
}

void Logger::removeMessageCombiner(const std::string& combinerIdentifier)
{
    int combinerToBeRemoved = -1;
    for (int i=0; i<(int)m_combiners.size(); ++i){
        if (m_combiners[i]->id() == combinerIdentifier){
            delete m_combiners[i];
            combinerToBeRemoved = i;
            break;
        }
    }
    if (combinerToBeRemoved != -1){
        m_combiners.erase(m_combiners.begin()+combinerToBeRemoved);
    }
}

void Logger::clearMessageCombiners()
{
    for (int i=0; i<(int)m_combiners.size(); ++i){
        delete m_combiners[i];
    }
    m_combiners.clear();
}

void Logger::getMessageInfos(std::vector<MessageInfo*>& filledVector)
{
    if (m_communicator->shouldMessagesBeOutputted()){
        filledVector.swap(m_messages);
    }
}

void Logger::ranksLimit(int ranksLimit)
{
    m_ranksLimit = ranksLimit;
}

int Logger::ranksLimit()
{
    return m_ranksLimit;
}

void Logger::pushMessagesOnce()
{
    m_communicator->pushMessagesOnce(m_messages);
    combineMessages();
}

void Logger::pushMessagesFully()
{
    m_communicator->pushMessagesFully(m_messages);
    combineMessages();
}

void Logger::queueMessage(const std::string& message)
{
    queueMessage(message, "", -1);
}

void Logger::queueMessage(const std::string& message, const std::string& fileName, const int lineNumber)
{
    MessageInfo* mi = new MessageInfo(message, m_communicator->rank(), fileName, lineNumber);
    m_messages.push_back(mi);
}

void Logger::combineMessages()
{
    int messagesSize = (int)m_messages.size();
    if (messagesSize < 2){
        return;
    }

    std::vector<MessageInfo*> finalMessageInfos;
    std::vector<int> indexesToBeDeleted;
    int combinersSize = (int)m_combiners.size();
    bool combinedMessageInfo = false;
    finalMessageInfos.push_back(m_messages[0]);
    for (int allIndex=1; allIndex<messagesSize; ++allIndex){
        combinedMessageInfo = false;
        for (int finalIndex=0; finalIndex<(int)finalMessageInfos.size(); ++finalIndex){
            for (int combinerIndex=0; combinerIndex<combinersSize; ++combinerIndex){
                if (m_combiners[combinerIndex]->shouldMessageInfosBeCombined(*finalMessageInfos[finalIndex],
                                                                             *m_messages[allIndex])){
                    m_combiners[combinerIndex]->combine(*finalMessageInfos[finalIndex],
                                                        *m_messages[allIndex], m_ranksLimit);
                    indexesToBeDeleted.push_back(allIndex);
                    combinedMessageInfo = true;
                    break;
                }
            }
            if (combinedMessageInfo){
                break;
            }
        }
        if (!combinedMessageInfo){
            finalMessageInfos.push_back(m_messages[allIndex]);
        }
    }

    for (int i=0; i<(int)indexesToBeDeleted.size(); ++i){
        delete m_messages[indexesToBeDeleted[i]];
    }
    m_messages.swap(finalMessageInfos);
}

}
}
