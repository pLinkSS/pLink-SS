#pragma once
#include <string>
#include <fstream>
#include "ProteinIndex.h"
#include "ProteinWriter.h"
#include "FastaProteinDB.h"
#include "ProteinHandler.h"

using namespace std;

namespace ProteinIndex{

	
class CDiskProteinWriter :public CProteinWriter
{

public:

	CDiskProteinWriter(void);
	~CDiskProteinWriter(void);
	
	virtual string GetWorkDir(void)const ;
	virtual string GetDBName(void)const;
	
	virtual void Init(const char * szDestDir, const char * szDBName,const char * szOrgDBPath, size_t &tMaxFileSize,PINDEX_HEAD &PIndexHead, PINDEX_ITEM&PIndexItem);
	virtual void WriteProteinIndex();	
	virtual void WriteProteinSuffix(size_t tSufType);
	
protected:

	// The original protein sequence database.
	string	m_strOrgDBPath;
	string	m_strWorkDir;
	string	m_strDBName;
	
	size_t	m_tMaxFileSize;
	CProteinHandler m_contr;
	CFastaProteinDB m_proDB;
	
	void _SetWorkDir(const char * szDestDir);
	void _SetOrgDBPath(const char * szPath);
	void _SetDBName(const char * szDBName);
	void _SetFileSize(size_t &tMaxFileSize);
	
	// Write the head of .pro.IDX ( offset of the head, protein number, original database file path )
	void _WriteHeadBlock(fstream& fIdx);

	void _FillPositionBlock(fstream& fIdx, size_t tProID);
	
	void _AppendEntryBlock(fstream& fIdx, fstream& fAC, fstream& fDE, fstream& fSQ);

	string& GetOrgDBPath(void);

private:
	PRO_IDX_HEAD m_stHead;
	PRO_IDX_POS_POINTER m_stProteinPosition;

public:
	PINDEX_HEAD m_pIndexHead;
	PINDEX_ITEM m_pIndexItem;
};
}
