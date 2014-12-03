#ifndef DISKRANDOMACREADER_H_
#define DISKRANDOMACREADER_H_
#include <string>
#include <cstdio>
#include <iostream>
#include "ProteinReader.h"
#include "ProteinHandler.h"

using namespace std;

class CProteinHandler;

namespace ProteinIndex
{

class CDiskRandomACReader: public CRandomReader
{
public:
	CDiskRandomACReader();
	virtual ~CDiskRandomACReader();
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
//	FILE* m_fpProIDX;
//	FILE* m_fpProDAT;

	CProteinHandler m_contr;
 
	unsigned char * m_pIDXBuf;
	unsigned char * m_pDATBuf;
	
	void _ReadPepSQ(string& strDesDAT, unsigned char cLen);
};

}

#endif /*DISKRANDOMACREADER_H_*/
