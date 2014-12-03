#ifndef METAHANDLER_H_
#define METAHANDLER_H_

typedef unsigned char uchar;


#include "../include/sdk.h"
#include "ProteinIndex.h"
#include "PepCalcFunc.h"

using namespace proteomics_sdk;

namespace ProteinIndex
{

class CMetaHandler
{
public:
	CMetaHandler();
	virtual ~CMetaHandler();

	void WriteMeta();
	void ReadMeta();	
	
	string GenerFileName(uchar uFileIdx,string strApp);
	
	
	
	size_t	m_tProNum;
	uchar	m_uFastaNum;
	vector <BLOCK_ITEM> m_vmItems;
	string m_strWorkDir, m_strDBName, m_szOrgDBPath;
//	short	m_uCurFileID;
	
//protected:	
public:
	
	
	size_t	m_HeadSize;	
	
	string _GenerFileName();
	
	void _WriteHead();
	void _ReadHead();
	
	void _WriteItem();
	void _ReadItem();	
};

}

#endif /*METAHANDLER_H_*/
