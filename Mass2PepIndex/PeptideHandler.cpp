
#include "PeptideHandler.h"

using namespace proteomics_sdk;

namespace Mass2PepIndex
{

CPeptideHandler::CPeptideHandler()
{
	m_fpIDX = NULL;
	m_fpDAT = NULL;
	m_uCurFileID = -1;
}

CPeptideHandler::~CPeptideHandler()
{
	Close();
}

void CPeptideHandler::Init(const string & strWorkDir, const string &strMetaName)
{
	try
	{
		m_metaReader.Load(strWorkDir, strMetaName);
		
		m_fpIDX = NULL;
		m_fpDAT = NULL;
		m_uCurFileID = -1;
		m_tPriorPepNum = 0;
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CPeptideHandler", "Init");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CPeptideHandler", "Init");		
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

void CPeptideHandler::Close()
{	
	try
	{	
//		if(MAP_FAILED != &m_MapIDX)
//		{
//			ACE_OS::close(m_MapIDX.handle());
//			m_MapIDX.unmap();
//			m_MapIDX.close();
//		}
//		if(MAP_FAILED != &m_MapDAT)
//		{
//			ACE_OS::close(m_MapDAT.handle());
//			m_MapDAT.unmap();
//			m_MapDAT.close();
//		}
		
		if(m_fpIDX)
		{
			fclose(m_fpIDX);
			m_fpIDX = NULL;
		}
		
		if(m_fpDAT)
		{
			fclose(m_fpDAT);
			m_fpDAT = NULL;
		}
		
		m_uCurFileID = -1;
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CPeptideHandler", "Close");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CPeptideHandler", "Close");		
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

size_t CPeptideHandler::GetMaxPepNum()
{
	return m_metaReader.GetMetaHead().tUniquePepSQNum;
}

string CPeptideHandler::GetDBName()
{
	return m_metaReader.GetMetaHead().cDBName;
}

size_t CPeptideHandler::GetMultiplier()
{
	return m_metaReader.GetMetaHead().tMultiplier;
}

string CPeptideHandler::GenerFileName(uchar uFileID,string strApp)
{
	try
	{
		string str = _GenerFileName();
		char sz[10];
		sprintf(sz,".%u", uFileID);
		return str + string(sz) + strApp;
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CPeptideHandler", "GenerFileName");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CPeptideHandler", "GenerFileName");		
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

string CPeptideHandler::_GenerFileName()
{
	return m_metaReader.GetWorkDir() + m_metaReader.GetMetaName().substr(0, m_metaReader.GetMetaName().size()- strlen(PEP_META ));
}

void CPeptideHandler::GetCurFileID(size_t tPepID, uchar &cDatNum, size_t &tCurID)
{
	try
	{
		size_t tPriorPepNum = 0;
		for(uchar  uFileID = 0; uFileID < m_metaReader.GetMetaHead().tIdxNum; ++uFileID)	
		{
			size_t tCurPepNum = m_metaReader.GetMetaItems()[uFileID].tUniquePepSQNum;
			tPriorPepNum += tCurPepNum;
			if(tPepID < tPriorPepNum)
			{
				cDatNum = uFileID;
				tCurID = tPepID - (tPriorPepNum - tCurPepNum);
				return;
			}
		}
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CPeptideHandler", "GetCurFileID");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CPeptideHandler", "GetCurFileID");		
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}	
}

size_t CPeptideHandler::GetCurrentPepFileEnd()
{
	return m_tPriorPepNum;
}

size_t CPeptideHandler::GetCurrentPepFileHead()
{
	return m_tHeadPepNum;
}
uchar CPeptideHandler::GetFileIDByMass(size_t tMass, long &tPass,long &te)//czhou
{
	try
	{
		for(uchar  uFileID = 0; uFileID < m_metaReader.GetMetaHead().tIdxNum; ++uFileID)	
		{
			if(tMass < m_metaReader.GetMetaItems()[uFileID].tMaxMass) 
			{
				te = m_metaReader.GetMetaItems()[uFileID].tUniquePepSQNum -1;
				return uFileID;
			}
			else tPass += m_metaReader.GetMetaItems()[uFileID].tUniquePepSQNum;
				
		}
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CPeptideHandler", "GetFileIDByMass");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CPeptideHandler", "GetFileIDByMass");		
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}		
	return m_metaReader.GetMetaHead().tIdxNum ;
}


size_t CPeptideHandler::CountIdxStructSize()
{
	return sizeof(size_t) + sizeof(long) ;
}

void CPeptideHandler::ReadIdxFromFile(size_t &tMass, long &lPos)
{	
	try
	{
		fread((char *)&tMass, sizeof(size_t), 1, m_fpIDX);
		fread((char *)&lPos, sizeof(long), 1, m_fpIDX);	
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CPeptideHandler", "ReadIdxFromFile");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CPeptideHandler", "ReadIdxFromFile");		
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

size_t CPeptideHandler::CountPepStructSize()// if there is a pattern designer to return this function automatic
{
	return sizeof(size_t) * 3 + sizeof(unsigned char) * 3;
}

//void CPeptideHandler::OpenAllFILE(ACE_Mem_Map * aDat, string strDat, size_t m_tDatNum)
//{
//	try
//	{
//		for(uchar t = 0; (size_t)t < m_tDatNum; ++t)
//		{
//			ACE_HANDLE handle = ACE_OS::open (GenerFileName(t,strDat).c_str(), O_RDONLY);
//			aDat[t].map(handle, (size_t)(-1),PROT_READ, ACE_MAP_SHARED);
//		}
//	}
//	catch(runtime_error & e)
//	{
//		CErrInfo info("CPeptideHandler", "OpenAllFILE");
//		cerr << info.Get(e) << endl;
//		throw runtime_error(info.Get(e).c_str());
//	}
//	catch(...)
//	{
//		CErrInfo info("CPeptideHandler", "OpenAllFILE");		
//		cerr << info.Get() << endl;
//		throw runtime_error(info.Get().c_str());
//	}	
//}
//
//void CPeptideHandler::CloseAllFILE(ACE_Mem_Map * aDat, string strDat, size_t m_tDatNum)
//{
//	try
//	{
//		for(uchar t = 0; (size_t)t < m_tDatNum; ++t)
//		{
//			ACE_OS::close(aDat[t].handle());
//			aDat[t].unmap();
//			aDat[t].close();
//		}
//	}
//	catch(runtime_error & e)
//	{
//		CErrInfo info("CPeptideHandler", "CloseAllFILE");
//		cerr << info.Get(e) << endl;
//		throw runtime_error(info.Get(e).c_str());
//	}
//	catch(...)
//	{
//		CErrInfo info("CPeptideHandler", "CloseAllFILE");		
//		cerr << info.Get() << endl;
//		throw runtime_error(info.Get().c_str());
//	}		
//}

void CPeptideHandler::OpenInvertedFiles(FILE ** aDat, string strDat, size_t m_tDatNum)
{
	try
	{
		for(uchar t = 0; (size_t)t < m_tDatNum; ++t)
		{
//			ACE_HANDLE handle = ACE_OS::open (GenerFileName(t,strDat).c_str(), O_RDONLY);
//			aDat[t].map(handle, (size_t)(-1),PROT_READ, ACE_MAP_SHARED);
			aDat[t] = fopen(GenerFileName(t,strDat).c_str(), "rb");
		}
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CPeptideHandler", "OpenAllFILE");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CPeptideHandler", "OpenAllFILE");		
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}	
}

void CPeptideHandler::CloseInvertedFiles(FILE ** aDat, string strDat, size_t m_tDatNum)
{
	try
	{
		for(uchar t = 0; (size_t)t < m_tDatNum; ++t)
		{
//			ACE_OS::close(aDat[t].handle());
//			aDat[t].unmap();
//			aDat[t].close();
			fclose(aDat[t]);
		}
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CPeptideHandler", "CloseAllFILE");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CPeptideHandler", "CloseAllFILE");		
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}		
}


bool CPeptideHandler::InitDiskDatFileByPepID(size_t tPepID, size_t tStartPos, string &strIDX, string &strDAT)
{
	try
	{
		m_tPriorPepNum = 0;
		m_tHeadPepNum = 0;
		for(uchar  uFileID = 0; uFileID < m_metaReader.GetMetaHead().tIdxNum; ++uFileID)	
		{
			size_t tCurPepNum = m_metaReader.GetMetaItems()[uFileID].tUniquePepSQNum;
			m_tHeadPepNum = m_tPriorPepNum;
			m_tPriorPepNum += tCurPepNum;
			if(tPepID < m_tPriorPepNum)
			{
				if(m_uCurFileID != uFileID)
				{	
					if(-1 != m_uCurFileID)
					{
						fclose(m_fpIDX);
						fclose(m_fpDAT);
					}

					m_fpIDX = fopen(GenerFileName(uFileID,strIDX).c_str(), "rb");
					m_fpDAT = fopen(GenerFileName(uFileID,strDAT).c_str(), "rb");
					
					if(!m_fpIDX || !m_fpDAT)
					{
						proteomics_sdk::CErrInfo err_info("CPeptideHandler", "InitDiskDatFileByPepID","The file don't exist!");
						throw runtime_error(err_info.Get().c_str());
					}
					
					m_uCurFileID = uFileID;
				}		
				
//				fseek(m_fpDAT,(tPepID - (m_tPriorPepNum -tCurPepNum)) * CountPepStructSize(), SEEK_SET);
				SeekToPlaceByPepID(tPepID); //added at 2014.3.18
				return true;
			}
		}
		return false;
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CPeptideHandler", "InitDiskDatFileByPepID");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CPeptideHandler", "InitDiskDatFileByPepID","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
	return true;
}
void CPeptideHandler::SeekToPlaceByPepID(size_t tPepID)
{
	fseek(m_fpDAT,(tPepID - (m_tPriorPepNum - m_metaReader.GetMetaItems()[m_uCurFileID].tUniquePepSQNum)) * CountPepStructSize(), SEEK_SET);
}

void CPeptideHandler::InitDiskDatFileByFileID(uchar uFileID, size_t tStartPos, string &strIDX, string &strDAT)
{
	try
	{
		if(m_uCurFileID != uFileID)
		{	
			if(-1 != m_uCurFileID)
			{
				fclose(m_fpIDX);
				fclose(m_fpDAT);
			}
			m_fpIDX = fopen(GenerFileName(uFileID,strIDX).c_str(), "rb");
			fseek(m_fpIDX, 0, SEEK_SET);
			
			m_fpDAT = fopen(GenerFileName(uFileID,strDAT).c_str(), "rb");
			fseek(m_fpDAT, 0, SEEK_SET);
		}				
		fseek(m_fpDAT,tStartPos,SEEK_SET);
		
		m_uCurFileID = uFileID;
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CPeptideHandler", "InitDiskDatFileByFileID");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CPeptideHandler", "InitDiskDatFileByFileID","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
}

//bool CPeptideHandler::InitMMapDatFileByPepID(size_t tPepID, size_t tStartPos, string &strIDX, string &strDAT)
//{
//	try
//	{
//		m_tPriorPepNum = 0;
//		for(uchar  uFileID = 0; uFileID < m_metaReader.GetMetaHead().tIdxNum; ++uFileID)	
//		{
//			size_t tCurPepNum = m_metaReader.GetMetaItems()[uFileID].tUniquePepSQNum;
//			m_tPriorPepNum += tCurPepNum;
//			if(tPepID < m_tPriorPepNum)
//			{
//				if(m_uCurFileID != uFileID)
//				{	
//					if(-1 != m_uCurFileID)
//					{
//						ACE_OS::close(m_MapIDX.handle());
//						m_MapIDX.unmap();
//						m_MapIDX.close();
//						
//						ACE_OS::close(m_MapDAT.handle());
//						m_MapDAT.unmap();
//						m_MapDAT.close();
//					}
//
//					
//					ACE_HANDLE handle = ACE_OS::open (GenerFileName(uFileID,strIDX).c_str(), O_RDONLY);
//					m_MapIDX.map(handle, (size_t)(-1),PROT_READ, ACE_MAP_SHARED);
//					handle = ACE_OS::open (GenerFileName(uFileID,strDAT).c_str(), O_RDONLY);
//					m_MapDAT.map(handle, (size_t)(-1),PROT_READ, ACE_MAP_SHARED);
//					
//					if(handle == ACE_INVALID_HANDLE || MAP_FAILED == &m_MapIDX || MAP_FAILED == &m_MapDAT )
//					{
//						proteomics_sdk::CErrInfo err_info("CPeptideHandler", "InitMMapDatFileByPepID","The file don't exist!");
//						throw runtime_error(err_info.Get().c_str());
//					}
//					
//					
//					m_uCurFileID = uFileID;
//				}		
//				ACE_OS::lseek(m_MapIDX.handle(), (tPepID - (m_tPriorPepNum -tCurPepNum)) * CountIdxStructSize(), SEEK_SET);
//				ACE_OS::lseek(m_MapDAT.handle(), (tPepID - (m_tPriorPepNum -tCurPepNum)) * CountPepStructSize(), SEEK_SET);
//				return true;
//			}
//		}
//		return false;
//	}
//	catch(runtime_error &e)
//	{
//		proteomics_sdk::CErrInfo err_info("CPeptideHandler", "InitMMapDatFileByProID");
//		throw runtime_error(err_info.Get(e).c_str());
//	}
//	catch(...)
//	{
//		proteomics_sdk::CErrInfo err_info("CPeptideHandler", "InitMMapDatFileByProID","Unknown Error!");
//		throw runtime_error(err_info.Get().c_str());
//	}
//	return true;
//}
//
//void CPeptideHandler::InitMMapDatFileByFileID(uchar uFileID, size_t tStartPos, string &strIDX, string &strDAT)
//{
//	try
//	{
//		if(m_uCurFileID != uFileID)
//		{	
//			if(-1 != m_uCurFileID)
//			{
//				ACE_OS::close(m_MapIDX.handle());
//				m_MapIDX.unmap();
//				m_MapIDX.close();
//				
//				ACE_OS::close(m_MapDAT.handle());
//				m_MapDAT.unmap();
//				m_MapDAT.close();
//			}
//			
//			ACE_HANDLE handle = ACE_OS::open (GenerFileName(uFileID,strIDX).c_str(), O_RDONLY);
//			m_MapIDX.map(handle, (size_t)(-1),PROT_READ, ACE_MAP_SHARED);
//			handle = ACE_OS::open (GenerFileName(uFileID,strDAT).c_str(), O_RDONLY);
//			m_MapDAT.map(handle, (size_t)(-1),PROT_READ, ACE_MAP_SHARED);
//			
//		}				
//
//		ACE_OS::lseek(m_MapDAT.handle(), tStartPos, SEEK_SET);
//		ACE_OS::lseek(m_MapIDX.handle(), 0, SEEK_SET);
//		
//		m_uCurFileID = uFileID;
//	}
//	catch(runtime_error &e)
//	{
//		proteomics_sdk::CErrInfo err_info("CPeptideHandler", "InitMMapDatFileByFileID");
//		throw runtime_error(err_info.Get(e).c_str());
//	}
//	catch(...)
//	{
//		proteomics_sdk::CErrInfo err_info("CPeptideHandler", "InitMMapDatFileByFileID","Unknown Error!");
//		throw runtime_error(err_info.Get().c_str());
//	}
//}


}
