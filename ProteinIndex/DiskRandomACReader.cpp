#include <string>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include "../include/predefine.h"
#include "ProteinIndex.h"
#include "ProteinReader.h"
#include "DiskRandomACReader.h"

using namespace std;

namespace ProteinIndex {

CDiskRandomACReader::CDiskRandomACReader() 
{
	m_tProteinID = (size_t)(-1);

	m_strIDXType = PRO_IDX_EXT;
	m_strDATType = PRO_AC_EXT;

	m_pIDXBuf = NULL;
	m_pDATBuf = NULL;
}

CDiskRandomACReader::~CDiskRandomACReader() 
{	
	Close();
}



void CDiskRandomACReader::Open(const std::string &strWorkDir,
		const std::string & strDBName) {
	try
	{
		m_contr.Init(strWorkDir, strDBName, "");
		m_contr.ReadMeta();
		m_contr.InitAllFile(m_strDATType);
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CDiskRandomACReader", "Open");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CDiskRandomACReader", "Open","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
}

void CDiskRandomACReader::Close() 
{
	m_tProteinID = (size_t)(-1);
	m_contr.Close();
}

bool CDiskRandomACReader::GetNext(string& strDAT) {
	if (m_tProteinID+1 >= m_tProNum) 
	{
		strDAT = "";
		return false;
	}
	++m_tProteinID;
	GetByID(strDAT, m_tProteinID);

	return true;
}

void CDiskRandomACReader::SetCurrentID(size_t tCurrentID)
{
	m_tProteinID = tCurrentID;
}

size_t CDiskRandomACReader::GetCurrentID() const
{
	return m_tProteinID;
}

size_t CDiskRandomACReader::GetProNum(void) const
{
	return m_contr.GetProNum();
}

size_t CDiskRandomACReader::GetByID(string& strDAT, size_t tProID) 
{
	try
	{
		m_contr.InitDiskDatFileByProID(tProID, 0, m_strIDXType, m_strDATType);

		size_t tDATCurPos = ftell(m_contr.m_fpProDAT);
		
		size_t tIdxNextPos = m_contr.GetEndPosInFile(tProID + 1) + sizeof(long) * 0;//czhou
		fseek(m_contr.m_fpProIDX,tIdxNextPos , SEEK_SET);
		size_t tDATNextPos = 0;
		fread((char *)&tDATNextPos,sizeof(size_t), 1, m_contr.m_fpProIDX);
		
		std::streamsize cLen = tDATNextPos - tDATCurPos;
		char* pszBuf1 = new char[cLen + 1];
		fread(pszBuf1, sizeof(char), cLen,  m_contr.m_fpProDAT);
		pszBuf1[cLen] = 0;
		strDAT = pszBuf1;
	
		delete []pszBuf1;

	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CDiskRandomACReader", "GetByID");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CDiskRandomACReader", "GetByID","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
	
	return true;
}

void CDiskRandomACReader::_ReadPepSQ(string& strDesDAT, unsigned char cLen)
{
	try
	{
		char* pszBuf = new char[cLen + 1];
		fread(pszBuf, sizeof(char), cLen, m_contr.m_fpProDAT);
		pszBuf[cLen] = '\0';
		strDesDAT = pszBuf;
		delete []pszBuf;	
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CDiskRandomACReader", "_ReadPepSQ");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CDiskRandomACReader", "_ReadPepSQ","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
}

void CDiskRandomACReader::ReadPepSQ(string& strDesSQ, size_t tProID, size_t tStartPos, unsigned char cLen) 
{
	m_contr.InitDiskDatFileByProID(tProID, tStartPos,  m_strIDXType, m_strDATType);
	_ReadPepSQ(strDesSQ, cLen);
}

//void CDiskRandomACReader::ReadPepSQ(string& strDesSQ, unsigned char cFileID, size_t tStartPos, unsigned char cLen)
//{
//	m_contr.InitDiskDatFileByFileID(cFileID, tStartPos, m_strIDXType, m_strDATType);
//	_ReadPepSQ(strDesSQ, cLen);
//}


void CDiskRandomACReader::ReadPepSQ_EX(string& strDesSQ, size_t tProID, size_t tStartPos, unsigned char cLen) 
{
	if(tStartPos==0)
	{
		ReadPepSQ(strDesSQ, tProID, tStartPos, cLen+1);
		strDesSQ.insert(0, "-.");
		strDesSQ.insert(strDesSQ.length()-1, ".");
		if(strDesSQ[strDesSQ.length()-1]==0)
			strDesSQ[strDesSQ.length()-1] = '-';
	}
	else
	{
		ReadPepSQ(strDesSQ, tProID, tStartPos-1, cLen+2);
		strDesSQ.insert(1,".");
		strDesSQ.insert(strDesSQ.length()-1, ".");
		if(strDesSQ[0]==0)
			strDesSQ[0] = '-';
		if(strDesSQ[strDesSQ.length()-1]==0)
			strDesSQ[strDesSQ.length()-1] = '-';
	}
}	

//void CDiskRandomACReader::ReadPepSQ_EX(string& strDesSQ, unsigned char cFileID, size_t tStartPos, unsigned char cLen)
//{
//	if(tStartPos==0)
//	{
//		ReadPepSQ(strDesSQ, cFileID, tStartPos, cLen+1);
//		strDesSQ.insert(0, "-.");
//		strDesSQ.insert(strDesSQ.length()-1, ".");
//		if(strDesSQ[strDesSQ.length()-1]==0)
//			strDesSQ[strDesSQ.length()-1] = '-';
//	}
//	else
//	{
//		ReadPepSQ(strDesSQ, cFileID, tStartPos-1, cLen+2);
//		strDesSQ.insert(1,".");
//		strDesSQ.insert(strDesSQ.length()-1, ".");
//		if(strDesSQ[0]==0)
//			strDesSQ[0] = '-';
//		if(strDesSQ[strDesSQ.length()-1]==0)
//			strDesSQ[strDesSQ.length()-1] = '-';
//	}
//}

}

