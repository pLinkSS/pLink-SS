#pragma once
#include <string>
#include "ProteinIndex.h"

using namespace std;
namespace ProteinIndex{
class CProteinWriter
{
public:

	CProteinWriter(void)
		: m_strWorkDir("./")
	{
	}

	virtual ~CProteinWriter(void)
	{
	}

	virtual string GetDBName(void)const =0;
	virtual string GetWorkDir(void)const = 0;

	// Load the original DB for to be indexed.
//	virtual void LoadOrgDB(const char * szDBName,const char * szOrgDBPath) = 0;
	virtual void WriteProteinIndex() = 0;
	virtual void WriteProteinSuffix(size_t tSufType) = 0;
	virtual void Init(const char * szDestDir, const char * szDBName,const char * szOrgDBPath, size_t &tMaxFileSize, PINDEX_HEAD &PIndexHead, PINDEX_ITEM&PIndexItem) = 0;
protected:

	// The default work directory for store the index files.
	string	m_strWorkDir;

	// The name (or path) for protein sequence database.
	string	m_strDBName;

};
}
