#include "MessageInfo.hpp"

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
	if (m_ranks.size() < ranksLimit){
		m_ranks.push_back(newRank);
	}
	m_rankCount++;
}

void MessageInfo::addRanks(const std::vector<int>& newRanks, int ranksLimit)
{
	int newRanksSize = newRanks.size();
	for(int i=0; i<newRanksSize; ++i){
		if (m_ranks.size() >= ranksLimit){
			break;
		}
		m_ranks.push_back(newRanks[i]);
	}
	m_rankCount += newRanksSize;
}

}
}
