#include "../include/sdk.h"
#include "../ProteinIndex/DiskRandomSQReader.h"
#include "../ProteinIndex/ProteinReaderFactory.h"

#include "DiskMass2PepIndexReader.h"

using namespace std;
using namespace ProteinIndex;
using namespace proteomics_sdk;

namespace Mass2PepIndex
{

CDiskMass2PepIndexReader::CDiskMass2PepIndexReader()
:m_tPepID(0), m_proLoader(NULL)
{
}

CDiskMass2PepIndexReader::~CDiskMass2PepIndexReader()
{
	Close();
}

void CDiskMass2PepIndexReader::Open(const std::string strWorkPath, const std::string strMetaName, ProteinRandomReaderType eType, bool bSQEX)
{
	if("" == strWorkPath || "" == strMetaName)
	{
		CErrInfo info("CDiskMass2PepIndexReader", "Open", "Strings should not be null!");
		
		info.Append("strWorkPath = " + strWorkPath);
		info.Append("strMetaName = " + strMetaName);
		
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
		
	try
	{
		CTrace::GetInstance()->Debug("Peptide Disk tag", MODULE_MASS2PEP);
		m_bSQEX = bSQEX;
		m_tPepID = 0;
		strIDX = PEP_INFO_IDX;
		strDAT = PEP_INFO_DAT;
		
		m_pepHandler.Init(strWorkPath, strMetaName);

		m_proLoader = m_ProRandomFactory.GetRandomSQReader(eType);
		m_proLoader->Open(strWorkPath, m_pepHandler.GetDBName());		

		m_tDatNum = 0;
	}
	catch(runtime_error &e)
	{
		CErrInfo info("CDiskMass2PepIndexReader", "Open");
		
		info.Append("strWorkPath = " + strWorkPath);
		info.Append("strMetaName = " + strMetaName);
		
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDiskMass2PepIndexReader", "Open", "Caught an unkown exception!");
		
		info.Append("strWorkPath = " + strWorkPath);
		info.Append("strMetaName = " + strMetaName);
		
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}

}

void CDiskMass2PepIndexReader::Close()
{
	try
	{
		m_tPepID = 0;
		if(m_proLoader != NULL)
		{
			delete m_proLoader;
			m_proLoader = NULL;
		}
		m_pepHandler.Close();	
		
		if(m_tDatNum)
		{	
			m_pepHandler.CloseInvertedFiles(aPos, PEP_INFO_POS, m_tDatNum);
			delete[] aPos;aPos = NULL;
			m_pepHandler.CloseInvertedFiles(aPid, PEP_INFO_PID, m_tDatNum);
			delete[] aPid;aPid = NULL;	
			
			m_tDatNum = 0;
		}
	}
	

	catch(runtime_error &e)
	{
		CErrInfo info("CDiskMass2PepIndexReader", "Close");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDiskMass2PepIndexReader", "Close", "Caught an unkown exception!");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}

}

bool CDiskMass2PepIndexReader::GetNext(PEP_SQ & stPepSQ)
{
	if(m_tPepID >= m_pepHandler.GetMaxPepNum())
	{
		char temp[200] = {0};
		sprintf(temp, "CDiskMass2PepIndexReader::GetNext()\n\tm_tPepID=%d\n\tm_pepHandler.GetMaxPepNum()=%d",m_tPepID, m_pepHandler.GetMaxPepNum());
		CTrace::GetInstance()->Debug(temp, MODULE_MASS2PEP);		
		
		return false;
	}
	
	try
	{
		GetByID(stPepSQ, m_tPepID);
		++m_tPepID;
	}
	catch(runtime_error &e)
	{
		char str[20] = {0};
		sprintf(str,"%d",m_tPepID);
		CErrInfo info("CDiskMass2PepIndexReader", "GetNext");
		info.Append("m_tPepID = " + (string)(str));
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		char str[20] = {0};
		sprintf(str,"%d",m_tPepID);
		CErrInfo info("CDiskMass2PepIndexReader", "GetNext", "Caught an unkown exception!");
		info.Append("m_tPepID = " + (string)(str));
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}

	return true;
}

bool CDiskMass2PepIndexReader::GetByID(PEP_SQ & stPepSQ, size_t tPepID)
{
	if(tPepID >= m_pepHandler.GetMaxPepNum())
	{
		char temp[200] = {0};
		sprintf(temp, "CDiskMass2PepIndexReader::GetByID()\n\tm_tPepID=%d\n\tm_pepHandler.GetMaxPepNum()=%d",m_tPepID, m_pepHandler.GetMaxPepNum());
		CTrace::GetInstance()->Debug(temp, MODULE_MASS2PEP);		
		
		return false;
	}
	
	try
	{
		PEP_INFO_EX_TAG stPepInfo;
		
		GetPepInfoByID(stPepInfo, tPepID);
		
		stPepSQ.dfMass = stPepInfo.tMass*1.0 / m_pepHandler.GetMultiplier();
		stPepSQ.cMiss = stPepInfo.cMiss;
		stPepSQ.cEnd = stPepInfo.cEnd;
		if(m_bSQEX)// true means *.sq.*, fasle means sq.
		{
			m_proLoader->ReadPepSQ_EX(stPepSQ.strSQ, stPepInfo.tProID, stPepInfo.tPos, stPepInfo.cLen);	
		}
		else
		{
			m_proLoader->ReadPepSQ(stPepSQ.strSQ, stPepInfo.tProID, stPepInfo.tPos, stPepInfo.cLen);//*****
		}
		//todo czhou
		m_pepHandler.GetCurFileID(tPepID, stPepSQ.cDatNum, stPepSQ.tInvertedFilesPos);
		stPepSQ.tInvertedFilesPos *= sizeof(size_t);	
	}
	catch(runtime_error &e)
	{
		char str[20];
		sprintf(str,"%d",tPepID);
		CErrInfo info("CDiskMass2PepIndexReader", "GetPepInfoByID", "PEP_SQ");
		info.Append("tPepID = " + (string)(str));
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		char str[20];
		sprintf(str,"%d",tPepID);
		CErrInfo info("CDiskMass2PepIndexReader", "GetPepInfoByID", "Caught an unkown exception(PEP_SQ)!");
		info.Append("tPepID = " + (string)(str));
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}

	return true;
}

bool CDiskMass2PepIndexReader::GetProIDByPos(uchar cDatNum, size_t tPos, vector<size_t> & vtProID)
{
	//todo vector need to be clear?
	try
	{
		if(false == m_pepHandler.m_metaReader.GetMetaHead().bPep2Pro) return false;//there is no pep2pro
		
		if(0 == m_tDatNum)
		{
			m_tDatNum = m_pepHandler.m_metaReader.GetMetaHead().tIdxNum;
			aPos = new FILE *[m_tDatNum];
			aPid = new FILE *[m_tDatNum];
			m_pepHandler.OpenInvertedFiles(aPos, PEP_INFO_POS, m_tDatNum);
			m_pepHandler.OpenInvertedFiles(aPid, PEP_INFO_PID, m_tDatNum);		
		}
		size_t tPoint = 0, tNum = 0;
		fseek(aPos[cDatNum],tPos, SEEK_SET);
		fread((char*)&tPoint, sizeof(size_t),1,aPos[cDatNum]);
		
		fseek(aPid[cDatNum], tPoint, SEEK_SET);
		fread((char*)&tNum, sizeof(size_t),1,aPid[cDatNum]);
//		ACE_OS::lseek(aPos[cDatNum].handle(),tPos, SEEK_SET);
//		ACE_OS::read( aPos[cDatNum].handle(), (char*)&tPoint, sizeof(size_t));
//		
//		ACE_OS::lseek(aPid[cDatNum].handle(),tPoint, SEEK_SET);
//		ACE_OS::read( aPid[cDatNum].handle(), (char*)&tNum, sizeof(size_t));
		
		for(size_t t = 0; t < tNum; ++t)
		{
			size_t tmp;
//			ACE_OS::read(aPid[cDatNum].handle(), (char *)&tmp, sizeof(size_t));
			fread((char *)&tmp, sizeof(size_t),1,aPid[cDatNum]);
			vtProID.push_back(tmp);
		}
	}
	catch(runtime_error &e)
	{
		CErrInfo info("CDiskMass2PepIndexReader", "GetProIDByPos");
		char str[20];
		sprintf(str,"%d",cDatNum);
		info.Append("cDatNum = " + (string)(str));
		sprintf(str,"%d",tPos);
		info.Append("tPos = " + (string)(str));
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDiskMass2PepIndexReader", "GetProIDByPos", "Caught an unkown exception(PEP_INFO_EX)!");
		char str[20];
		sprintf(str,"%d",cDatNum);
		info.Append("cDatNum = " + (string)(str));
		sprintf(str,"%d",tPos);
		info.Append("tPos = " + (string)(str));
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}	
	return true;
}

bool CDiskMass2PepIndexReader::GetPepInfoByID(PEP_INFO_EX_TAG & stPepInfo, size_t tPepID)
{
	if(tPepID >= m_pepHandler.GetMaxPepNum())
	{
		char temp[200] = {0};
		sprintf(temp, "CDiskMass2PepIndexReader::GetPepInfoByID()\n\tm_tPepID=%d\n\tm_pepHandler.GetMaxPepNum()=%d",m_tPepID, m_pepHandler.GetMaxPepNum());
		CTrace::GetInstance()->Debug(temp, MODULE_MASS2PEP);		
		return false;
	}
	
	try
	{
//		if(tPepID >= m_pepHandler.GetCurrentPepFileEnd() || m_tPepID != tPepID )
//		{
//			m_pepHandler.InitDiskDatFileByPepID(tPepID, 0, strIDX, strDAT);
//		}
		if(tPepID >= m_pepHandler.GetCurrentPepFileEnd()  || tPepID < m_pepHandler.GetCurrentPepFileHead() ) // 2014.3.18 Change file and jump
		{
			m_pepHandler.InitDiskDatFileByPepID(tPepID, 0, strIDX, strDAT);
		}
		else
		{
			//todo: 对于GetNext调用，可以不执行else之中定位语句，以提高速度。但目前接口中未留判断是GetNext还是GetByID的参数。可考虑利用tPepID的最高位信息。
			m_pepHandler.SeekToPlaceByPepID(tPepID);
		}

		_ReadPepFromFile(stPepInfo);
	}
	catch(runtime_error &e)
	{
		char str[20];
		sprintf(str,"%d",tPepID);
		CErrInfo info("CDiskMass2PepIndexReader", "GetPepInfoByID", "PEP_INFO");
		info.Append("tPepID = " + (string)(str));
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		char str[20];
		sprintf(str,"%d",tPepID);
		CErrInfo info("CDiskMass2PepIndexReader", "GetPepInfoByID", "Caught an unkown exception(PEP_INFO)!");
		info.Append("tPepID = " + (string)(str));
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
	
	return true;
}



void CDiskMass2PepIndexReader::MoveTo(double dMass)
{	
	m_tPepID = GetIDByMass(dMass);
	m_pepHandler.InitDiskDatFileByPepID(m_tPepID, 0, strIDX, strDAT);
}

size_t CDiskMass2PepIndexReader::GetIDByMass(double dMass)
{
	if(dMass < 0) dMass = 0;
	size_t tMass = (size_t)(dMass * m_pepHandler.GetMultiplier());
	
	long ts = 0;//(long)m_tPepID;
	long te = 0;
	long tPass = 0;
	uchar uFileID = m_pepHandler.GetFileIDByMass(tMass,tPass,te);
	if(m_pepHandler.m_metaReader.GetMetaHead().tIdxNum == uFileID)
	{
		return m_pepHandler.GetMaxPepNum();
	}
	
	m_pepHandler.InitDiskDatFileByFileID(uFileID, 0, strIDX, strDAT);
	
	long tm = (ts + te) >> 1;
	long tp = -1;
	
	size_t tCurMass = 0;
	
	while(ts <= te)
	{
		tm = (ts + te) >> 1;
		
		fseek(m_pepHandler.m_fpDAT, tm * m_pepHandler.CountPepStructSize(), SEEK_SET);
		fread((char *)&tCurMass, sizeof(size_t), 1, m_pepHandler.m_fpDAT);

		if(tCurMass < tMass)
		{
			tp = tm;
			ts = tm + 1;
		}
		else te= tm - 1;
	}
	
	return ++tp + tPass;	
}

//size_t CDiskMass2PepIndexReader::GetIDByMass(double dMass)
//{
//	if(dMass < 0) dMass = 0;
//	size_t tMass = (size_t)(dMass * m_pepHandler.GetMultiplier());
//	return m_pepHandler.GetIDByMass(tMass);
//}

inline size_t gsize(unsigned char *ch)
{
	return (((size_t)ch[3]) << 24) + (((size_t)ch[2]) << 16) + (((size_t)ch[1]) << 8) + (size_t)ch[0];
}

void CDiskMass2PepIndexReader::_ReadPepFromFile(PEP_INFO_EX_TAG & sDatExTagBlock)
{
	try
	{		
//		size_t cLen = m_pepHandler.CountPepStructSize();
//		unsigned char chVal[cLen];
//		fread(chVal, sizeof(char), cLen, m_pepHandler.m_fpDAT );
//		
//		sDatExTagBlock.tMass = gsize(chVal);
//		sDatExTagBlock.tProID = gsize(chVal+4);
//		sDatExTagBlock.tPos = gsize(chVal+8);
//		sDatExTagBlock.cLen = (unsigned char)chVal[12];
//		sDatExTagBlock.cMiss = (unsigned char)chVal[13];
		
		fread((char *)&sDatExTagBlock.tMass, sizeof(size_t), 1, m_pepHandler.m_fpDAT);
		fread((char *)&sDatExTagBlock.tProID, sizeof(size_t), 1, m_pepHandler.m_fpDAT);
		fread((char *)&sDatExTagBlock.tPos, sizeof(size_t), 1, m_pepHandler.m_fpDAT);
		fread((char *)&sDatExTagBlock.cLen, sizeof(unsigned char), 1, m_pepHandler.m_fpDAT);
		fread((char *)&sDatExTagBlock.cMiss, sizeof(unsigned char), 1, m_pepHandler.m_fpDAT);
		fread((char *)&sDatExTagBlock.cEnd, sizeof(unsigned char), 1, m_pepHandler.m_fpDAT);
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CDiskMass2PepIndexReader", "_ReadPepFromFile");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDiskMass2PepIndexReader", "_ReadPepFromFile");		
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

}
