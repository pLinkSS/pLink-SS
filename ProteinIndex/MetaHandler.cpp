#include "MetaHandler.h"

namespace ProteinIndex
{

CMetaHandler::CMetaHandler()
{
}

CMetaHandler::~CMetaHandler()
{
}

string CMetaHandler::GenerFileName(uchar uFileID,string strApp)
{
	try
	{
		string str = _GenerFileName();
		char sz[10];
		sprintf(sz,".%u", uFileID);
		return str + string(sz) + strApp;
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CMetaHandler", "GenerFileName");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CMetaHandler", "GenerFileName","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
}

string CMetaHandler::_GenerFileName()
{
	return m_strWorkDir + m_strDBName;
}


void CMetaHandler::WriteMeta()
{
	_WriteHead();
	_WriteItem();
}

void CMetaHandler::ReadMeta()
{
	_ReadHead();
	_ReadItem();
}

void CMetaHandler::_WriteHead()
{
	try
	{
		string strMeta = _GenerFileName() + PRO_META_EXE;
		string strTxt = strMeta + ".txt";
		
		ofstream fMeta(strMeta.c_str(), ios::out | ios::trunc | ios::binary );
		ofstream fTxt(strTxt.c_str(), ios::out | ios::trunc );
		
		fMeta.write((char *)&m_HeadSize, sizeof(size_t));
		fMeta.write((char *)&m_uFastaNum, sizeof(uchar));
		fMeta.write((char *)&m_tProNum, sizeof(size_t));		
		
		fMeta.write(m_szOrgDBPath.c_str() , m_szOrgDBPath.size());
		char abc[PATH_MAX + 1];
		memset(abc,0, PATH_MAX + 1);
		fMeta.write(abc , PATH_MAX - m_szOrgDBPath.size());

		m_HeadSize = fMeta.tellp();
		fMeta.seekp(0, ios::beg);
		fMeta.write((char *)&m_HeadSize, sizeof(size_t));
		fMeta.close();
		
		fTxt << "[Head]" << endl;
		fTxt << "FastaNum = " << (int)(m_uFastaNum )<< endl;
		fTxt << "ProNum = " << m_tProNum << endl;
		fTxt << "OrgDBPath = " << m_szOrgDBPath << endl;
		fTxt.close();
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CMetaHandler", "_WriteHead");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CMetaHandler", "_WriteHead","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
}

void CMetaHandler::_ReadHead()
{
	try
	{
		string strMeta = _GenerFileName() + PRO_META_EXE;
	
		ifstream fMeta(strMeta.c_str(), ios::binary );
		fMeta.seekg(0, ios::beg);
		fMeta.read((char *)&m_HeadSize, sizeof(size_t));
		fMeta.read((char *)&m_uFastaNum, sizeof(uchar));
		fMeta.read((char *)&m_tProNum, sizeof(size_t));
		char *szOrgDBPath = new char[PATH_MAX + 1];
		fMeta.read(szOrgDBPath, sizeof(char) * PATH_MAX);//****
		m_szOrgDBPath = szOrgDBPath;
		fMeta.close();
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CMetaHandler", "_ReadHead");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CMetaHandler", "_ReadHead","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
}

void CMetaHandler::_WriteItem()
{
	try
	{
		string strMeta = _GenerFileName() + PRO_META_EXE;
		string strTxt = strMeta + ".txt";
		
		ofstream fMeta(strMeta.c_str(), ios::out  | ios::app |ios::binary );
		ofstream fTxt(strTxt.c_str(), ios::out | ios::app );
		
		fMeta.seekp(0, ios::end);
		fTxt.seekp(0, ios::end);
		
		size_t tDatSize = m_vmItems.size();
		for(size_t t = 0; t < tDatSize; ++t)
		{
			fMeta.write((char*)&m_vmItems[t].tBeg, sizeof(size_t));
			fMeta.write((char*)&m_vmItems[t].tEnd, sizeof(size_t));
			
			char sz[10];
			sprintf(sz,"%d",t);
			fTxt << "[dat" << string(sz) << "]" << endl;
			
			fTxt << "First Pro = " << m_vmItems[t].tBeg << endl;
			fTxt << "Last Pro = " << m_vmItems[t].tEnd << endl;
		}
		
		fMeta.close();
		fTxt.close();
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CMetaHandler", "_WriteItem");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CMetaHandler", "_WriteItem","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
}

void CMetaHandler::_ReadItem()
{
	try
	{
		string strMeta = _GenerFileName() + PRO_META_EXE;
		
		ifstream fMeta(strMeta.c_str(), ios::binary );	
		fMeta.seekg(m_HeadSize, ios::beg);
		m_vmItems.clear();
		BLOCK_ITEM mItem;
		for(size_t t = 0; t < m_uFastaNum; ++t)
		{		
			fMeta.read((char*)&mItem.tBeg, sizeof(size_t));
			fMeta.read((char*)&mItem.tEnd, sizeof(size_t));
			m_vmItems.push_back(mItem);
		}
		fMeta.close();
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CMetaHandler", "_ReadItem");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CMetaHandler", "_ReadItem","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
}


}
