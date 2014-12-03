#ifndef NTERMONLINEDIGESTION_H_
#define NTERMONLINEDIGESTION_H_

#include "../include/sdk.h"
#include "PepCalcFunc.h"
#include "OnlineDigestion.h"
#include "MetaHandler.h"


namespace ProteinIndex
{

class CNTermOnlineDigestion : public ProteinIndex::COnlineDigestion
{
public:
	CNTermOnlineDigestion();
	virtual ~CNTermOnlineDigestion();
	
	virtual bool InitPara(PINDEX_HEAD &pIndexHead,PINDEX_ITEM &pIndexItem, CMetaHandler &pMetaHandler) ;
	virtual bool Close() ;
	virtual bool GetNextPep(PEP_SQ &strPepSQ) ;
	virtual bool GetNextBlock() ;
	virtual bool SuffixSort(size_t) ;
	virtual bool GetPep2Pro(size_t tp, size_t tLen, vector<size_t> & vtProID);
	
	inline bool _IsDiv(size_t p);
	inline bool _IsEnzyme(size_t p);
	inline bool _IsLeftTerm(size_t p);
	inline bool _IsRightTerm(size_t p);
	inline bool _IsLeft(size_t p);
	inline bool _IsRight(size_t p);
	inline bool _IsPep(size_t s, size_t e);


	void _DeleteAll();
	void _InitDigest();
	
	bool GetNextAllPro(size_t tSufType);
	bool GetNextBeg();
	bool GetNextLCP();
	bool GetNextSA();

	void _Print(size_t s, size_t e, uchar t,PEP_SQ &strPepSQ);	
	bool _GetSemiSpecificPep1(PEP_SQ &strPepSQ);
	bool _GetSemiSpecificPep2(PEP_SQ &strPepSQ);
	
	bool GetSpecificPep(PEP_SQ &strPepSQ);
	bool GetSemiSpecificPep(PEP_SQ &strPepSQ);
	bool GetNonSpecificPep(PEP_SQ &strPepSQ);
	
	size_t GetProID(size_t tp);
	long long StaticNonSpecificPep();
	
	void _GetLCP();
	void _DealSpecial();	
public:
	short m_uCurFileID;
	CMetaHandler m_pMetaHandler;
	
public:
	PINDEX_HEAD m_pIndexHead;
	PINDEX_ITEM m_pIndexItem;
	
	size_t m_tLen, m_tNum; 
	
	uchar  *m_chDat, *LCP;
	size_t *SA;
	size_t *nClvy, ncp;
	size_t *m_tMass;
	size_t *beg;
	size_t m_b, m_e, m_c, tbg, tMaxClvy;
	
	CEnzyme	enzyme;
	size_t *m_tAAMass;
	
	string m_strCleave;
	string m_strNotCleave;
};
}

#endif /*NTERMONLINEDIGESTION_H_*/
