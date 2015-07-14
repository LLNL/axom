#include "MessageInfo.hpp"
#include "Utility.hpp"

namespace asctoolkit {
namespace lumberjack {

std::string MessageInfo::message() const
{
    return m_message;
}

void MessageInfo::message(const std::string& newMessage)
{
    m_message = newMessage;
}

std::vector<int> MessageInfo::ranks() const
{
    return m_ranks;
}

int MessageInfo::rankCount() const
{
    return m_rankCount;
}

std::string MessageInfo::fileName() const
{
    return m_fileName;
}

std::string MessageInfo::stringOfRanks(std::string delimiter) const
{
    std::string returnString = "";
    int ranksSize = m_ranks.size();
    for(int i=0; i<ranksSize;++i){
        returnString += intToString(m_ranks[i]);
        if (i < (ranksSize-1)) {
            returnString += delimiter;
        }
    }
    return returnString;
}

void MessageInfo::fileName(const std::string& newFileName)
{
    m_fileName = newFileName;
}

int MessageInfo::lineNumber() const
{
    return m_lineNumber;
}

void MessageInfo::lineNumber(int newLineNumber)
{
    m_lineNumber = newLineNumber;
}

void MessageInfo::addRank(int newRank, int ranksLimit)
{
    if (m_ranks.size() < (std::vector<int>::size_type)ranksLimit){
        m_ranks.push_back(newRank);
    }
    m_rankCount++;
}

void MessageInfo::addRanks(const std::vector<int>& newRanks, int ranksLimit)
{
    int newRanksSize = newRanks.size();
    for(int i=0; i<newRanksSize; ++i){
        if (m_ranks.size() >= (std::vector<int>::size_type)ranksLimit){
            break;
        }
        m_ranks.push_back(newRanks[i]);
    }
    m_rankCount += newRanksSize;
}

const char memberDelimiter = '*';
const char rankDelimiter = ',';

std::string MessageInfo::pack()
{
    std::string packedMessageInfo;
    int ranksSize = (int)m_ranks.size();
    for (int i=0; i<ranksSize; ++i){
        packedMessageInfo += intToString(m_ranks[i]);
        if (i < (ranksSize-1)) {
            packedMessageInfo += rankDelimiter;
        }
    }
    packedMessageInfo += memberDelimiter + intToString(m_rankCount);
    packedMessageInfo += memberDelimiter + m_fileName + memberDelimiter;

    if (m_lineNumber != -1){
        packedMessageInfo += intToString(m_lineNumber);
    }
    packedMessageInfo += memberDelimiter + m_message;

    return packedMessageInfo;
}

void MessageInfo::unpack(const std::string& packedMessage, int ranksLimit)
{
    int messageLength = (int)packedMessage.length();
    std::string currString;
    int i = 0;

    // Grab ranks
    m_ranks.clear();
    int currRank;
    for (i=0; i<messageLength; ++i) {
        if ((packedMessage[i] == memberDelimiter) ||
            (packedMessage[i] == rankDelimiter)) {
            currRank = stringToInt(currString);
            currString = "";
            addRank(currRank, ranksLimit);
            if (packedMessage[i] == memberDelimiter) {
                ++i;
                break;
            }
            else {
                continue;
            }
        }
        currString += packedMessage[i];
    }

    //Grab rank count since it can differ from list that is sent
    currString = "";
    for (; i<messageLength; ++i) {
        if (packedMessage[i] == memberDelimiter) {
            m_rankCount = stringToInt(currString);
            ++i;
            break;
        }
        currString += packedMessage[i];
    }

    //Grab file name
    currString = "";
    for (; i<messageLength; ++i) {
        if (packedMessage[i] == memberDelimiter) {
            m_fileName = currString;
            ++i;
            break;
        }
        currString += packedMessage[i];
    }

    //Grab line number
    currString = "";
    for (; i<messageLength; ++i) {
        if (packedMessage[i] == memberDelimiter) {
            m_lineNumber = stringToInt(currString);
            ++i;
            break;
        }
        currString += packedMessage[i];
    }

    //Grab message
    currString = "";
    for (; i<messageLength; ++i) {
        currString += packedMessage[i];
    }
    m_message = currString;
}

}
}
