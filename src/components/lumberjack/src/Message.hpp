#ifndef MESSAGE_HPP
#define MESSAGE_HPP

#include <string>
#include <vector>

namespace asctoolkit {
namespace lumberjack {

class Message {
	public:
		// Constructors
		Message()
		: m_message("")
		, m_rankCount(0)
		, m_fileName("")
		, m_lineNumber(0) {}

		Message(const std::string& message, int rank, const std::string& fileName, int lineNumber)
		: m_message(message)
		, m_rankCount(1)
		, m_fileName(fileName)
		, m_lineNumber(lineNumber)
		{
			m_ranks.push_back(rank);
		}

		Message(const std::string& message, std::vector<int> ranks, const std::string& fileName, int lineNumber)
		: m_message(message)
		, m_ranks(ranks)
		, m_rankCount(ranks.size())
		, m_fileName(fileName)
		, m_lineNumber(lineNumber) {}

		// Getters
		std::string message() const;
		std::vector<int> ranks() const;
		int rankCount() const;
		std::string fileName() const;
		int lineNumber() const;

		// Setters/Utilities
		void message(const std::string& newMessage);
		void fileName(const std::string& newFileName);
		void lineNumber(int newLineNumber);
		void addRank(int newRank, int rankLimit);
		void addRanks(const std::vector<int>& newRanks, int rankLimit);
	private:
		std::string m_message;
		std::vector<int> m_ranks;
		int m_rankCount;
		std::string m_fileName;
		int m_lineNumber;
};

}
}

#endif