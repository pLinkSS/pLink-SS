#include <string>
#include <stdexcept>
#include <iostream>
#include "../include/predefine.h"
#include "ProteinIndex.h"
#include "ProteinReader.h"
#include "FastaSQReader.h"

using namespace std;
using namespace proteomics_sdk;
namespace ProteinIndex
{

CFastaSQReader::CFastaSQReader()
{
	m_tProteinID = (size_t)(-1);
}

CFastaSQReader::~CFastaSQReader()
{
	Close();
}
void CFastaSQReader::Open(const std::string &strWorkDir,const std::string & strDBName) 
{
	Close();
	string strFastaName = strWorkDir + strDBName + FILE_FASTA;

	m_ifIn.open(strFastaName.c_str());
	if(m_ifIn.fail())
	{
		proteomics_sdk::CErrInfo err_info("CFastaSQReader", "Open","Open fasta file error!");
		throw runtime_error(err_info.Get().c_str());
	}
	return;
}

void CFastaSQReader::Close()
{
	if(m_ifIn.is_open())
	{
		m_ifIn.clear();
		m_ifIn.close();
	}
}
bool CFastaSQReader::GetNext(string &strSQ)
{
	CProtein protein;
	
	try
	{
		if(!m_Parser.ReadOnePrteinEntry(protein, m_ifIn))
		{
			return false;
		}
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CFastaSQReader", "GetNext");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CFastaSQReader", "GetNext","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}

	++m_tProteinID;
	strSQ = protein.m_strSQ;

	return true;
}
size_t CFastaSQReader::GetCurrentID() const
{
	return m_tProteinID;
}
}
