#include "Message.hpp"

namespace asctoolkit {
namespace lumberjack {

std::string Message::message() const
{
	return m_message;
}

void Message::message(const std::string& newMessage)
{
	m_message = newMessage;
}

std::vector<int> Message::ranks() const
{
	return m_ranks;
}

int Message::rankCount() const
{
	return m_rankCount;
}

std::string Message::fileName() const
{
	return m_fileName;
}

void Message::fileName(const std::string& newFileName)
{
	m_fileName = newFileName;
}

int Message::lineNumber() const
{
	return m_lineNumber;
}

void Message::lineNumber(int newLineNumber)
{
	m_lineNumber = newLineNumber;
}

void Message::addRank(int newRank, int rankLimit)
{
	if (m_ranks.size() < rankLimit){
		m_ranks.push_back(newRank);
	}
	m_rankCount++;
}

void Message::addRanks(const std::vector<int>& newRanks, int rankLimit)
{
	int newRanksSize = newRanks.size();
	for(int i=0; i<newRanksSize; ++i){
		if (m_ranks.size() >= rankLimit){
			break;
		}
		m_ranks.push_back(newRanks[i]);
	}
	m_rankCount += newRanksSize;
}

}
}
