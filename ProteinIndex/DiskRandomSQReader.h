#ifndef DISKRANDOMSQREADER_H_
#define DISKRANDOMSQREADER_H_

#include "ProteinReader.h"
#include "ProteinHandler.h"

using namespace std;

namespace ProteinIndex {

class CDiskRandomSQReader : public CRandomReader {
public:
	CDiskRandomSQReader();
	virtual ~CDiskRandomSQReader();

	void Open(const std::string &strWorkDir, const std::string & strDBName);
	void Close();

	bool GetNext(string& strProtein);
	size_t GetByID(string& strDAT, size_t proID) ;
	
	size_t GetCurrentID() const;
	void SetCurrentID(size_t tCurrentID);
	
	size_t GetProNum(void) const;
	
	void ReadPepSQ(string& strDesDAT, size_t tProID, size_t tStartPos, unsigned char cLen) ;
	//void ReadPepSQ(string& strDesDAT, unsigned char cFileID, size_t tStartPos, unsigned char cLen) ;

	void ReadPepSQ_EX(string& strDesDAT, size_t tProID, size_t tStartPos, unsigned char cLen) ;
	//void ReadPepSQ_EX(string& strDesDAT, unsigned char cFileID, size_t tStartPos, unsigned char cLen);

private:
	string m_strIDXType ;
	string m_strDATType ;

	CProteinHandler m_contr;
 
	unsigned char * m_pIDXBuf;
	unsigned char * m_pDATBuf;
	
	void _ReadPepSQ(string& strDesDAT, unsigned char cLen);
};

}

#endif /*DISKRANDOMSQREADER_H_*/
