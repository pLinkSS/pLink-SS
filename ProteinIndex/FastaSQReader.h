#ifndef FASTAREADER_H_
#define FASTAREADER_H_
#include <string>
#include "FastaParser.h"
using namespace std;
using namespace proteomics_sdk;
namespace ProteinIndex
{

class CFastaSQReader : public CProteinReader
{
public:
	CFastaSQReader();
	~CFastaSQReader();
	void Open(const std::string &strWorkDir, const std::string & strDBName);
	void Close();
	bool GetNext(string &strSQ);
	size_t GetCurrentID() const;

private:
	ifstream m_ifIn ;
	CFastaParser m_Parser;
	size_t m_tProteinID;
};

}

#endif /*FASTAREADER_H_*/
