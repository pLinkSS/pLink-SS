#include <string>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include "../include/sdk.h"
#include "DiskRandomSQReader.h"

using namespace std;

namespace ProteinIndex {

CDiskRandomSQReader::CDiskRandomSQReader() 
{
	m_tProteinID = (size_t)(-1);

	m_strIDXType = PRO_IDX_EXT;
	m_strDATType = PRO_SQ_EXT;
	
	m_pIDXBuf = NULL;
	m_pDATBuf = NULL;
}

CDiskRandomSQReader::~CDiskRandomSQReader() 
{	
	Close();
}



void CDiskRandomSQReader::Open(const std::string &strWorkDir,
		const std::string & strDBName) {

	try
	{
		CTrace::GetInstance()->Debug("Protein Disk tag", MODULE_PROIDX);
		m_contr.Init(strWorkDir, strDBName, "");
		CTrace::GetInstance()->Debug("Read meta", MODULE_PROIDX);
		m_contr.ReadMeta();
		CTrace::GetInstance()->Debug("Init All File", MODULE_PROIDX);
		m_contr.InitAllFile(m_strDATType);
		CTrace::GetInstance()->Debug("InitDiskDatFile", MODULE_PROIDX);
		m_contr.InitDiskDatFileByProID(0, 0, m_strIDXType, m_strDATType);
		CTrace::GetInstance()->Debug("Open SQ File Completed", MODULE_PROIDX);
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CDiskRandomSQReader", "Open");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CDiskRandomSQReader", "Open","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
}

void CDiskRandomSQReader::Close() 
{
	m_tProteinID = (size_t)(-1);
	m_contr.Close();

}

bool CDiskRandomSQReader::GetNext(string& strDAT) 
{
	if (m_tProteinID+1 >= m_tProNum)
	{
		strDAT = "";
		return false;
	}
	++m_tProteinID;
	GetByID(strDAT, m_tProteinID);
	if("" == strDAT)
	{
		CTrace::GetInstance()->Alert("There are empty proteins", MODULE_MASS2PEP);
		return GetNext(strDAT);
	}
	return true;
}

void CDiskRandomSQReader::SetCurrentID(size_t tCurrentID)
{
	m_tProteinID = tCurrentID;
}

size_t CDiskRandomSQReader::GetCurrentID() const
{
	return m_tProteinID;
}

size_t CDiskRandomSQReader::GetProNum(void) const
{
	return m_contr.GetProNum();
}

size_t CDiskRandomSQReader::GetByID(string& strDAT, size_t tProID) 
{
	try
	{
		m_contr.InitDiskDatFileByProID(tProID, 0, m_strIDXType, m_strDATType);

		size_t tDATCurPos = ftell(m_contr.m_fpProDAT);
		
		size_t tIdxNextPos = m_contr.GetEndPosInFile(tProID + 1) + sizeof(long) * 2;
		fseek(m_contr.m_fpProIDX,tIdxNextPos , SEEK_SET);
		size_t tDATNextPos = 0;
		fread((char *)&tDATNextPos,sizeof(size_t), 1, m_contr.m_fpProIDX);
		
		std::streamsize cLen = tDATNextPos - tDATCurPos;
		char* pszBuf1 = new char[cLen + 1];
		fread(pszBuf1, sizeof(char), cLen,  m_contr.m_fpProDAT);
		pszBuf1[cLen] = 0;
		//czhou to deal with DiskProteinWriter, every string end with a '@',that useful for suffix
		if(pszBuf1[cLen - 1] <'A' || pszBuf1[cLen - 1] > 'Z') pszBuf1[cLen - 1] = 0;
		strDAT = pszBuf1;
	
		delete []pszBuf1;
		
		if(strDAT[0] && (strDAT[0] <'A' || strDAT[0] > 'Z'))
		{ 
			strDAT = strDAT.substr(1, strDAT.size() -1);
			return tDATCurPos + 1;
		}
	
		return tDATCurPos;	

	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CDiskRandomSQReader", "GetByID");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CDiskRandomSQReader", "GetByID","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
	return 0;
}

void CDiskRandomSQReader::_ReadPepSQ(string& strDesDAT, unsigned char cLen)
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
		proteomics_sdk::CErrInfo err_info("CDiskRandomSQReader", "_ReadPepSQ");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CDiskRandomSQReader", "_ReadPepSQ","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
}

void CDiskRandomSQReader::ReadPepSQ(string& strDesSQ, size_t tProID, size_t tStartPos, unsigned char cLen) 
{
	m_contr.InitDiskFileByPepPos(tProID, tStartPos);
//	m_contr.InitDiskFileByPepPos(tProID, tStartPos, m_strIDXType, m_strDATType);
	_ReadPepSQ(strDesSQ, cLen);
}

//void CDiskRandomSQReader::ReadPepSQ(string& strDesSQ, unsigned char cFileID, size_t tStartPos, unsigned char cLen)
//{
//	m_contr.InitDiskDatFileByFileID(cFileID, tStartPos,  m_strIDXType, m_strDATType);
//	_ReadPepSQ(strDesSQ, cLen);
//}


void CDiskRandomSQReader::ReadPepSQ_EX(string& strDesSQ, size_t tProID, size_t tStartPos, unsigned char cLen) 
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

//void CDiskRandomSQReader::ReadPepSQ_EX(string& strDesSQ, unsigned char cFileID, size_t tStartPos, unsigned char cLen)
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
