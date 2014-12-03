#ifndef PROTEINHANDLER_H_
#define PROTEINHANDLER_H_



#include "../include/sdk.h"

#include "PepCalcFunc.h"
#include "MetaHandler.h"
#include "ProteinIndex.h"
#include "OnlineDigestion.h"

#define DIVPEP '?'
typedef unsigned char uchar;

using namespace std;


namespace ProteinIndex
{




class CProteinHandler
{
public:
	CProteinHandler();
	virtual ~CProteinHandler();
	
	void Init(const string & strWorkDir, const string &strDBName,const string &szOrgDBPath = "");
	void Close();
	void _DeleteAll();

	void ObtOStream(fstream &fIdx, fstream &fAC, fstream &fDE, fstream &fSQ, uchar uFileID);

	void InitPara(	PINDEX_HEAD m_pIndexHead,PINDEX_ITEM m_pIndexItem );

	inline bool IsDiv(size_t p);
	inline bool IsEnzyme(size_t p);
	inline bool IsLeftTerm(size_t p);
	inline bool IsRightTerm(size_t p);
	inline bool IsLeft(size_t p);
	inline bool IsRight(size_t p);
	inline bool IsPep(size_t s, size_t e);
	
	bool GetNextBlock();
	bool GetNextAllPro(size_t tSufType);
	bool GetNextBeg();
	bool GetNextLCP();
	bool GetNextSA();
	void InitDigest();
	
	void Print(size_t s, size_t e, uchar t,PEP_SQ &strPepSQ);
	long long StaticNonSpecificPep();
	
	size_t GetProID(size_t tp);
	bool GetPep2Pro(size_t tp, size_t tLen, vector<size_t> & vtProID);
	
	bool GetSpecificPep(PEP_SQ &strPepSQ);
	bool GetSemiSpecificPep(PEP_SQ &strPepSQ);
	bool GetNonSpecificPep(PEP_SQ &strPepSQ);
	bool GetNextPep(PEP_SQ &strPepSQ);
	
	bool _GetSemiSpecificPep1(PEP_SQ &strPepSQ);
	bool _GetSemiSpecificPep2(PEP_SQ &strPepSQ);
	
	bool SuffixSort(size_t tSufType);
	void _GetLCP();
	void _DealSpecial();
	
	void WriteMeta();
	void ReadMeta();	
	
	size_t GetProNum() const;
	size_t GetPosInFile(size_t tProID);
	size_t GetEndPosInFile(size_t tProID);
	string GenerFileName(uchar uFileIdx,string strApp);
	void InsertItem(BLOCK_ITEM);
	void SetProInfor(size_t &tProNum, uchar & uFastaNum);
	
	//bool InitDiskFileByPepPos(size_t tProID, size_t tStartPos, string &strIDX, string &strDAT);
	bool InitDiskDatFileByProID(size_t tProID, size_t tStartPos, string &strIDX, string &strDat);
	//void InitDiskDatFileByFileID(uchar uFileID, size_t tStartPos, string &strIDX, string &strDat);

//	bool InitMMapFileByPepPos(size_t tProID, size_t tStartPos, string &strIDX, string &strDAT);
//	bool InitMMapDatFileByProID(size_t tProID, size_t tStartPos,string &strIDX, string &strDat);
//	void InitMMapDatFileByFileID(uchar uFileID, size_t tStartPos,string &strIDX, string &strDat);
	
	void InitAllFile(const string & strDAT)
	{
		for(uchar uFileID = 0; uFileID < m_vmItems.size(); ++uFileID)		
		{
			m_DAT[uFileID] = fopen(GenerFileName(uFileID,strDAT).c_str(), "rb");
		}
	}
	
	bool InitDiskFileByPepPos(size_t tProID, size_t tStartPos)
	{
		for(uchar  uFileID = 0; uFileID < m_vmItems.size(); ++uFileID)		
		{
			if( tProID < m_vmItems[uFileID].tEnd)
			{
				m_uCurFileID = uFileID;
				fseek(m_DAT[uFileID], tStartPos,SEEK_SET);		
				m_fpProDAT = m_DAT[uFileID];
				return true;
			}
		}
		return false;
		
	}
	
	FILE * m_fpProIDX;
	FILE * m_fpProDAT;
	
//	ACE_Mem_Map m_MapIDX;
//	ACE_Mem_Map m_MapDAT;
	
	FILE * m_IDX[50];
	FILE * m_DAT[50];
	
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
	
public:
	size_t	m_tProNum;
	uchar	m_uFastaNum;
	size_t	m_HeadSize;
	vector <BLOCK_ITEM> m_vmItems;

	short	m_uCurFileID;
	
	string m_strWorkDir, m_strDBName, m_szOrgDBPath;
	
	string _GenerFileName();
	
	void _WriteHead();
	void _ReadHead();
	
	void _WriteItem();
	void _ReadItem();	
	
public:
	COnlineDigestion *m_pDigester; 
	
};

}

#endif /*PROTEINHANDLER_H_*/
