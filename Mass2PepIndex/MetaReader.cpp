
#include "MetaReader.h"

using namespace std;
using namespace proteomics_sdk;

namespace Mass2PepIndex
{

CMetaReader::CMetaReader():m_fMeta(NULL)
{
}

CMetaReader::~CMetaReader()
{
	CloseFile();
}
void	CMetaReader::Load( const string strWorkPath, const string strMetaName)
{
	if("" == strWorkPath || "" == strMetaName)
	{
		CErrInfo info("CMetaReader", "Load", "Strings should not be null!");
		
		info.Append("strWorkPath = " + strWorkPath);
		info.Append("strMetaName = " + strMetaName);
		
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
	try
	{
		SetWorkDir(strWorkPath.c_str());
		SetMetaName(strMetaName); 
		
		OpenFile();
	
		SetMetaHead();
		SetMetaItems();
	}
	catch(runtime_error &e)
	{
		CErrInfo info("CMetaReader", "Load");
		
		info.Append("strWorkPath = " + strWorkPath);
		info.Append("strMetaName = " + strMetaName);
		
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CMetaReader", "Load", "Caught an unkown exception!");
		
		info.Append("strWorkPath = " + strWorkPath);
		info.Append("strMetaName = " + strMetaName);
		
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

void	CMetaReader::OpenFile(void)
{
	string strMetaName = m_strWorkDir + m_strMetaName;
	try
	{		
		m_fMeta = fopen(strMetaName.c_str(), "rb");
	}
	catch(runtime_error &e)
	{
		CErrInfo info("CMetaReader", "OpenFile");
		info.Append("strMetaName = " + strMetaName);		
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CMetaReader", "OpenFile", "Caught an unkown exception!");
		info.Append("strMetaName = " + strMetaName);	
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}
void	CMetaReader::CloseFile(void)
{
	try
	{
		if(m_fMeta)
		{
			fclose(m_fMeta);
			m_fMeta = NULL;
		}
	}
	catch(runtime_error &e)
	{
		CErrInfo info("CMetaReader", "CloseFile");		
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CMetaReader", "CloseFile", "Caught an unkown exception!");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

void	CMetaReader::SetMetaName(string strMetaName)
{
	m_strMetaName = strMetaName;
}

void	CMetaReader::SetMetaHead()
{
	try
	{
		fseek(m_fMeta, 0, SEEK_SET);
		fread((char *)&m_stMetaHead, _GetMetaHeadSize(), 1, m_fMeta);
	}
	catch(runtime_error &e)
	{
		CErrInfo info("CMetaReader", "SetMetaHead");		
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CMetaReader", "SetMetaHead", "Caught an unkown exception!");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

void	CMetaReader::SetMetaItems()
{
	try
	{
		m_vstMetaItems.clear();
		META_ITEM stMetaItem;
		fseek(m_fMeta, m_stMetaHead.tOffset, SEEK_SET);
		
		for(size_t t=0; t<m_stMetaHead.tIdxNum; ++t)
		{
			fread((char *)&stMetaItem, _GetMetaItemSize(), 1, m_fMeta);
			
//			fread((char *)&stMetaItem.tPepSQNum, sizeof(size_t), 1, m_fMeta);
//			fread((char *)&stMetaItem.tUniquePepSQNum, sizeof(size_t), 1, m_fMeta);
//			fread((char *)&stMetaItem.tUniqueMassNum, sizeof(size_t), 1, m_fMeta);
//			fread((char *)&stMetaItem.tMinMass, sizeof(size_t), 1, m_fMeta);
//			fread((char *)&stMetaItem.tMaxMass, sizeof(size_t), 1, m_fMeta);
//			
			m_vstMetaItems.push_back(stMetaItem);
		}	
	}
	catch(runtime_error &e)
	{
		CErrInfo info("CMetaReader", "SetMetaItems");		
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CMetaReader", "SetMetaItems", "Caught an unkown exception!");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

void 	CMetaReader::SetWorkDir(const char * szDestDir)
{
	try
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
			m_strWorkDir = strDestDir.substr(0,idx + 1) + "/";
		}
		else
			m_strWorkDir = strDestDir;
	}
	catch(runtime_error &e)
	{
		CErrInfo info("CMetaReader", "SetWorkDir");		
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CMetaReader", "SetWorkDir", "Caught an unkown exception!");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

size_t	CMetaReader::_GetMetaHeadSize(void)const
{
	//read the offset part
	size_t tOffset;
	try
	{
		fseek(m_fMeta, 0, SEEK_SET);
		fread((char *)&tOffset, sizeof(size_t), 1, m_fMeta);
		fseek(m_fMeta, 0, SEEK_SET);
	}
	catch(runtime_error &e)
	{
		CErrInfo info("CMetaReader", "_GetMetaHeadSize");		
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CMetaReader", "_GetMetaHeadSize", "Caught an unkown exception!");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
	return tOffset;
}

size_t	CMetaReader::_GetMetaItemSize(void)const
{
	return sizeof(META_ITEM);
}

}
