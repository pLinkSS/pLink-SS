#ifndef INDEXER_H_
#define INDEXER_H_

#include "../Mass2PepIndex/Mass2PepIndex.h"
#include "../Mass2PepIndex/Mass2PepIndexCreatorFactory.h"

#include "../ProteinIndex/ProteinWriterFactory.h"
#include "../ProteinIndex/DiskRandomSQReader.h"
#include "../ProteinIndex/ProteinReaderFactory.h"

#include "../include/sdk.h"
#include "../include/predefine.h"
#include "../include/option.h"

#include "ACEConfigTool.h"
#include "DBConf.h"
#include "Configure.h"
#include "ProteinIndexWriter.h"

using namespace ProteinIndex;
using namespace Mass2PepIndex;
using namespace std;


class CExcuteIndexer
{
public:
	CExcuteIndexer();
	virtual ~CExcuteIndexer();
	
	void Init(ProteinWriterType eProteinWriterType, IndexCreatorType ePeptideCreatorType, string &strPIndexName);
	
	void Run(void);
	
	void Tester(void);
	
	void Reader(void);
	
	void Timer(void);
protected:
	CProteinWriter* m_proteinIndexerCreator;
	CProteinIndexWriter m_proteinWriter;
	
	CMass2PepIndexCreator * m_peptideIndexerCreator;
	
	CConfigure m_configure;
	
	void _WriteConf(PINDEX_HEAD &pIndexHead, PINDEX_ITEM& pIndexItem, size_t, bool);
	void _WriteMeta(PINDEX_HEAD &pIndexHead, PINDEX_ITEM& pIndexItem);
};

#endif /*INDEXER_H_*/
