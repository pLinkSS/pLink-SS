#ifndef PROTEININDEXWRITER_H_
#define PROTEININDEXWRITER_H_

#include "ACEConfigTool.h"
#include "DBConf.h"

#include "../include/sdk.h"
#include "../include/predefine.h"
#include "../include/option.h"

#include "../ProteinIndex/ProteinWriterFactory.h"
#include "../ProteinIndex/DiskRandomSQReader.h"
#include "../ProteinIndex/ProteinReaderFactory.h"

using namespace std;
using namespace ProteinIndex;
using namespace proteomics_sdk;


class CProteinIndexWriter
{
public:
	CProteinIndexWriter();
	~CProteinIndexWriter();
	
	void Init(const PINDEX_HEAD& pIndexHead, const PINDEX_ITEM& pIndexItem);
	void Write(CProteinWriter* m_pProIdxWriter);
	
protected:
	PINDEX_HEAD m_PIndexHead ;
	PINDEX_ITEM m_PIndexItem ;
	
	size_t m_tDepth;
	
	void	 _WriteTargetDecoyData(string strSourcePath, string strNewPath);	
		
	CProtein _ReadNext(ifstream &m_ifIn);
	string	_ReadEntrySQ(ifstream &m_ifIn, char*  pszBuf);
	void	_WriteProtein(ofstream &ofFile, const CProtein &protein, const size_t nSQSize, const size_t line);

	void	_ChangeInvalid(string& strContent);
	void 	_ChangeChar(string& strContent, string chChanged);

	void	_ShowFBeginInfo(string strFName);

	void	_ShowFEndInfo(string strFName, double dfTime);
};

#endif /*PROTEININDEXWRITER_H_*/
