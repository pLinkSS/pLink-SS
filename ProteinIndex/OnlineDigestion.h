#ifndef ONLINEDIGESTION_H_
#define ONLINEDIGESTION_H_

#define DIVPEP '?'
typedef unsigned char uchar;

//#include "ProteinIndex.h"
#include "MetaHandler.h"

namespace ProteinIndex
{

class COnlineDigestion
{
public:
	virtual ~COnlineDigestion(){};
	virtual bool InitPara(PINDEX_HEAD &pIndexHead,PINDEX_ITEM &pIndexItem, CMetaHandler &pMetaHandler) = 0;
	virtual bool Close() = 0;
	virtual bool GetNextPep(PEP_SQ &strPepSQ) = 0;
	virtual bool GetNextBlock() = 0;
	virtual bool SuffixSort(size_t) = 0;
	virtual bool GetPep2Pro(size_t tp, size_t tLen, vector<size_t> & vtProID) = 0;
};

}

#endif /*ONLINEDIGESTION_H_*/
