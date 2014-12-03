#include <iostream>
#include "FastaParser.h"
#include "FastaProteinDB.h"

using namespace ProteinIndex;

CFastaProteinDB::CFastaProteinDB(void)
{
	m_strDBName = "";
	m_strPath = "";
}

CFastaProteinDB::~CFastaProteinDB(void)
{
	CloseFile();
}

void CFastaProteinDB::SetDBPath(const string& filepath)
{
	m_strPath = filepath;
}

void CFastaProteinDB::SetDBName(const string& dbname)
{
	m_strDBName = dbname;
}

string& CFastaProteinDB::GetDBPath(void)
{
	return m_strPath;
}

string& CFastaProteinDB::GetDBName(void)
{
	return m_strDBName;
}

// Get current protein entry that read from sequence database.
CProtein& CFastaProteinDB::GetProEntry(void)
{
	return m_proEntry;
}

CProtein& CFastaProteinDB::ReadOnePrteinEntry()
{
	m_proEntry.m_strAC.clear();
	m_proEntry.m_strDE.clear();
	m_proEntry.m_strSQ.clear();

	m_Parser.ReadOnePrteinEntry(m_proEntry, m_ifIn);

	return m_proEntry;
}

// Note: 
//		Before using this function, must call OpenFile() earlier,
//		and after this function must call CloseFile().
//	 like:
//		if ( !OpenFile(strPath) ) 
//			return false;
//		....
//		CProtein protein;
//      if(!ReadOnePrteinEntry(protein))
//        ShowError();
//      else strAC = protein.m_strAC;
//		....
//		CloseFile();
bool CFastaProteinDB::ReadOnePrteinEntry(proteomics_sdk::CProtein& rPro)
{
	try
	{
		if (!m_Parser.ReadOnePrteinEntry(rPro, m_ifIn))
			return false;
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CFastaProteinDB", "ReadOnePrteinEntry");
		throw runtime_error(err_info.Get(e).c_str());
				
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CFastaProteinDB", "ReadOnePrteinEntry","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
	return true;
}


// Open file stream
void CFastaProteinDB::OpenFile(const char* szFilePath)
{
	m_strPath = szFilePath;
	try
	{
		OpenFile();
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CFastaProteinDB", "OpenFile");
		err_info.Append("m_strPath="+m_strPath);
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CFastaProteinDB", "OpenFile","Unknown Error!");
		err_info.Append("m_strPath="+m_strPath);
		throw runtime_error(err_info.Get().c_str());
	}
}

void CFastaProteinDB::OpenFile()
{
	CloseFile();

	if (GetDBPath().empty())
	{
		proteomics_sdk::CErrInfo err_info("CFastaProteinDB", "OpenFile","DBPath empty!");
		throw runtime_error(err_info.Get().c_str());
	}
	m_ifIn.open(GetDBPath().c_str());

	if(m_ifIn.fail())
	{
		proteomics_sdk::CErrInfo err_info("CFastaProteinDB", "OpenFile","Can't find fasta file!");
		throw runtime_error(err_info.Get().c_str());
	}
}
// Close file stream
// Close file stream
void CFastaProteinDB::CloseFile(void)
{
	if(m_ifIn.is_open())
	{
		m_ifIn.clear();
		m_ifIn.close();
	}
}

// get the number of all the protein entries in the sequence database
size_t CFastaProteinDB::GetEntryNum(void)
{
//	m_ifIn.clear();
	return m_Parser.GetEntryNum(m_ifIn);

}

void CFastaProteinDB::SeekBegin(void)
{
	m_ifIn.clear(); 
	m_ifIn.seekg(0,std::ios::beg);
}
