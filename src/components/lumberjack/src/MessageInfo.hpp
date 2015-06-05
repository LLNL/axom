#ifndef MESSAGEINFO_HPP
#define MESSAGEINFO_HPP

#include <string>
#include <vector>

namespace asctoolkit {
namespace lumberjack {

class MessageInfo {
	public:
		// Constructors
		MessageInfo()
		: m_message("")
		, m_rankCount(0)
		, m_fileName("")
		, m_lineNumber(0) {}

		MessageInfo(const std::string& message, int rank,
			        const std::string& fileName, int lineNumber)
		: m_message(message)
		, m_rankCount(1)
		, m_fileName(fileName)
		, m_lineNumber(lineNumber)
		{
			m_ranks.push_back(rank);
		}

		MessageInfo(const std::string& message, std::vector<int> ranks,
			        const std::string& fileName, int lineNumber)
		: m_message(message)
		, m_ranks(ranks)
		, m_rankCount(ranks.size())
		, m_fileName(fileName)
		, m_lineNumber(lineNumber) {}

		// Getters
		std::string message() const;
		std::vector<int> ranks() const;
		std::vector<int>::size_type rankCount() const;
		std::string fileName() const;
		int lineNumber() const;

		// Setters/Utilities
		void message(const std::string& newMessage);
		void fileName(const std::string& newFileName);
		void lineNumber(int newLineNumber);
		void addRank(int newRank, std::vector<int>::size_type ranksLimit);
		void addRanks(const std::vector<int>& newRanks, std::vector<int>::size_type ranksLimit);
	private:
		std::string m_message;
		std::vector<int> m_ranks;
		std::vector<int>::size_type m_rankCount;
		std::string m_fileName;
		int m_lineNumber;
};

}
}

#endif