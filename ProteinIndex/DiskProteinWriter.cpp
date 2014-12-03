#include <iostream>
#include "FastaParser.h"
#include "FastaProteinDB.h"
#include "DiskProteinWriter.h"
#include "File_ptr.h"

using namespace std;
using namespace ProteinIndex;

CDiskProteinWriter::CDiskProteinWriter(void)
:m_strOrgDBPath("")
{
}

CDiskProteinWriter::~CDiskProteinWriter(void)
{
}


// Load DB for to be indexed.
void CDiskProteinWriter::Init(const char * szDestDir, const char * szDBName,const char * szOrgDBPath, size_t &tMaxFileSize,PINDEX_HEAD &pIndexHead, PINDEX_ITEM &pIndexItem)
{
	try
	{
		_SetWorkDir(szDestDir);
		_SetDBName( szDBName );
		_SetOrgDBPath( szOrgDBPath );
		_SetFileSize(tMaxFileSize);
	
		m_contr.Init(m_strWorkDir, m_strDBName, m_strOrgDBPath);
		
		m_proDB.SetDBName(GetDBName());
		m_proDB.SetDBPath(GetOrgDBPath());
		
		m_pIndexHead = pIndexHead;
		m_pIndexItem = pIndexItem;
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CDiskProteinWriter", "Init");
		throw runtime_error(err_info.Get(e).c_str());	
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CDiskProteinWriter", "Init","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());			
	}
}

void CDiskProteinWriter::WriteProteinSuffix(size_t tSufType)
{
	m_contr.Init(m_strWorkDir, m_strDBName, m_strOrgDBPath);
	m_contr.InitPara(m_pIndexHead, m_pIndexItem);
	m_contr.SuffixSort(tSufType);
	m_contr.Close();
}

void CDiskProteinWriter::WriteProteinIndex()
{
	if (GetDBName().empty())
	{
		proteomics_sdk::CErrInfo err_info("CDiskProteinWriter", "WriteProteinIndex","DBName empty!");
		throw runtime_error(err_info.Get().c_str());
	}		

	if (GetOrgDBPath().empty() && GetDBName().find( "\\" ) == string::npos && GetDBName().find( "/" ) == string::npos) 
	{
		proteomics_sdk::CErrInfo err_info("CDiskProteinWriter", "WriteProteinIndex","DBPath error!");
		err_info.Append("OrgDBPath="+GetOrgDBPath());
		err_info.Append("DBName="+GetDBName());
		throw runtime_error(err_info.Get().c_str());
	}
	
	//check if proDB is open
	try
	{
		m_proDB.OpenFile();
		
		fstream fIdx, fAC, fDE, fSQ;
		uchar uFileID = 0;

		size_t tProID = 0, tAllProID = 0;	// 0 ~ ProNum - 1
		
		
		m_contr.ObtOStream(fIdx, fAC, fDE, fSQ, uFileID++);
		
		BLOCK_ITEM mItem;
		mItem.tBeg = 0;
		while ( !m_proDB.IsEof() ) 
		{
			m_proDB.ReadOnePrteinEntry();
			
			size_t tlen = m_proDB.GetProEntry().m_strSQ.length();
			size_t tlen2 = fSQ.tellp();

			if(m_tMaxFileSize < tlen + tlen2 )
			{
				mItem.tEnd = tAllProID;
				m_contr.InsertItem(mItem);
				
				m_stProteinPosition.lPosAC = fAC.tellp();
				m_stProteinPosition.lPosDE = fDE.tellp();
				m_stProteinPosition.lPosSQ = fSQ.tellp();
				_FillPositionBlock( fIdx,tProID );					
				
				m_contr.ObtOStream(fIdx, fAC, fDE, fSQ, uFileID++);				
				mItem.tBeg = tAllProID;
				tProID = 0;
			}
			
			_AppendEntryBlock( fIdx, fAC, fDE, fSQ);
			_FillPositionBlock( fIdx,tProID );

			++tProID;
			++tAllProID;
		}
		
		mItem.tEnd = tAllProID;
		m_contr.InsertItem(mItem);
		
		// Important: write The No.tNumProEntry+1 items in POINTER BLOCK. The last item stores the end position of file.
		m_stProteinPosition.lPosAC = fAC.tellp();
		m_stProteinPosition.lPosDE = fDE.tellp();
		m_stProteinPosition.lPosSQ = fSQ.tellp();
		_FillPositionBlock( fIdx,tProID );	

		m_proDB.CloseFile();

		fIdx.close();
		fAC.close();
		fDE.close();
		fSQ.close();
		
		m_contr.SetProInfor(tAllProID, uFileID);
		m_contr.WriteMeta();	
		m_contr.Close();
			
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CDiskProteinWriter", "WriteProteinIndex");
		err_info.Append("m_proDB.GetDBName()="+m_proDB.GetDBName());
		throw runtime_error(err_info.Get(e).c_str());	
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CDiskProteinWriter", "WriteProteinIndex","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());			
	}
}

void CDiskProteinWriter::_WriteHeadBlock(fstream& fIdx)
{
	try
	{
		fIdx.seekp( 0, ios::beg);
	
		m_stHead.lProIdxHeadOffset = 0;
		m_stHead.tProteinNum = m_proDB.GetEntryNum();
		strncpy(m_stHead.szOrgDBPath,GetOrgDBPath().c_str(),PATH_MAX);
		
		//The first section of .pro.IDX is the offset of the head
		fIdx.write((char * )&m_stHead.lProIdxHeadOffset , sizeof(long));
		// The second section of .pro.IDX is the number of protein.
		fIdx.write((char * )&m_stHead.tProteinNum , sizeof(size_t));
		// The third section of .pro.IDX is the original protein database path( MAX_PATH 260 chars ).
		//char[MAX_PATH] szPath = {'\0'};
		fIdx.write((char * )&m_stHead.szOrgDBPath , sizeof(char) * PATH_MAX);
	
		//rewrite the offset
		m_stHead.lProIdxHeadOffset = fIdx.tellp();
		fIdx.seekp( 0, ios::beg);
		fIdx.write((char * )&m_stHead.lProIdxHeadOffset , sizeof(long));
	
		fIdx.seekp( 0, ios::end);
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CDiskProteinWriter", "_WriteHeadBlock");
		throw runtime_error(err_info.Get(e).c_str());	
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CDiskProteinWriter", "_WriteHeadBlock","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());			
	}
}

// Fill one item in the POSITION BLOCK with the starting position of AC,DE,SQ .
void CDiskProteinWriter::_FillPositionBlock(fstream& fIdx, size_t tProID)
{
	try
	{
		fIdx.write((char * )&m_stProteinPosition.lPosAC , sizeof(long));
		fIdx.write((char * )&m_stProteinPosition.lPosDE , sizeof(long));	
		fIdx.write((char * )&m_stProteinPosition.lPosSQ , sizeof(long));
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CDiskProteinWriter", "_FillPositionBlock");
		throw runtime_error(err_info.Get(e).c_str());	
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CDiskProteinWriter", "_FillPositionBlock","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());			
	}	
}

// Append the protein's information( AC,DE,SQ,etc.) in ENTRY BLOCK.
void CDiskProteinWriter::_AppendEntryBlock(fstream& fIdx, fstream& fAC, fstream& fDE, fstream& fSQ)
{
	try
	{
		m_stProteinPosition.lPosAC = fAC.tellp();
		fAC.write((char * )m_proDB.GetProEntry().m_strAC.c_str(), sizeof(char) * long(m_proDB.GetProEntry().m_strAC.length()) + 1);
		
		m_stProteinPosition.lPosDE = fDE.tellp();
	//	if(m_strDE[0])
		string strDE = m_proDB.GetProEntry().m_strDE;
		if("" == strDE) strDE = m_proDB.GetProEntry().m_strAC;
		if('>' == strDE[0]) strDE = strDE.substr(1, strDE.size() - 1);
		fDE.write(strDE.c_str(), sizeof(char) * long(strDE.length()) + 1);
		
		m_stProteinPosition.lPosSQ = fSQ.tellp();
		if(0 == m_stProteinPosition.lPosSQ) 
			m_proDB.GetProEntry().m_strSQ = '?' + m_proDB.GetProEntry().m_strSQ;
		m_proDB.GetProEntry().m_strSQ += '?';
		fSQ.write((char * )m_proDB.GetProEntry().m_strSQ.c_str(), sizeof(char) * long(m_proDB.GetProEntry().m_strSQ.length()) );
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CDiskProteinWriter", "_AppendEntryBlock");
		throw runtime_error(err_info.Get(e).c_str());	
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CDiskProteinWriter", "_AppendEntryBlock","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());			
	}		
}

string& CDiskProteinWriter::GetOrgDBPath(void)
{
	return m_strOrgDBPath;
}

string CDiskProteinWriter::GetDBName(void) const
{
	return m_strDBName;
}

string CDiskProteinWriter::GetWorkDir(void) const 
{
	return m_strWorkDir;
}

void CDiskProteinWriter::_SetWorkDir(const char * szDestDir)
{
	string strDestDir( szDestDir );
	string::size_type idx = 0;
	idx = strDestDir.find('\\', idx);
	while (idx != string::npos ) {
		strDestDir[idx] = '/';
		idx = strDestDir.find('\\', idx);
	}

	idx = strDestDir.find_last_not_of(" ");
	if (strDestDir[idx] != '/' ) 
	{
		m_strWorkDir = strDestDir.substr(0,idx+1) + "/";
	}
	else
		m_strWorkDir = strDestDir;
		return;
}

void CDiskProteinWriter::_SetDBName(const char * szDBName)
{
	m_strDBName = szDBName;
}

void CDiskProteinWriter::_SetOrgDBPath(const char * szPath)
{
	m_strOrgDBPath = szPath;	
}

void CDiskProteinWriter::_SetFileSize(size_t &tMaxFileSize)
{
	 m_tMaxFileSize= tMaxFileSize;	
}

