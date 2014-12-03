#ifndef XLINKOPENSCORER_H_
#define XLINKOPENSCORER_H_

#include "CommonClass.h"

typedef unsigned char UCHAR;

namespace proteomics_sdk
{
	
class CXLinkOpenScorer
{
public:
	CXLinkOpenScorer();
	virtual ~CXLinkOpenScorer();
	virtual void Initialize(const CCondition & condition);
	virtual void Close(void);
	virtual void SetSpectrum(CSpectrum & spec);
	virtual void SetPeptide(CXLinkOpenPepResult & pep);
	virtual double Score();
	void _Match();
	void _ComputeMZ();
	inline void _GetMassBorder(int nMz, int nChg, int &nMin, int &nMax);
	void _SetCondition(const CCondition & cond);
	void _CalMatchInfo();
	void _SetPeakRank();
	
	void SetPeptide(CXLinkOpenMiddlePepResult & pep); //added at 2014.4.28

	size_t GetMatchNo(); //added at 2014.5.13. Return how many theory peaks are matched

protected:
	static void _InitLogNPerm();
	static void _InitLogTagLenCDF();
	double _CalProb(int N, int n, int n0, int l, int x, int x1);
	double _CalProb(int N, int n, int l, int x);
	double _C(int D,int U);
	int _ShiftMass(size_t tIonId);
	double GetMaxTagLength();
	
protected:	
	
	//members:
	CCondition m_Condition;
	CXLinkOpenPepResult * m_pPeptide;
	CSpectrum * m_pSpectrum;
	size_t m_tionmzSize;
	size_t m_tIonNum;
	bool m_iontag_pep[MAX_IONTYPE_NUM][MAX_PEPTIDE_LENGTH];
	
	int * m_pnMatched;
	size_t * m_pPeakRank;
	size_t m_nMaxPeaksNum;
	CMzTriple * m_ionmz;
	
	//aamass precomputing
	int m_nAAMass[26];
	bool m_bComputeMz;
	//0 %;1 mmu;2 ppm;3 da
	char m_cTolType;
	//to compute the multiplier according to the m_cTolType
	double m_lfTolMultiplier;
	// modify by emily
	char m_cTolBaseType;
	double m_lfTolBaseMultiplier;
	
	static double m_lfLogNPerm[100000];
	static double m_lfLogStdNormCDF[101];
	static double m_lfLogGamaCDF[100];
	static double m_lfLogTagLenCDF[120][21];
	//FILE * m_fp;

	CXLinkOpenMiddlePepResult m_pMidPeptide; //注意，此类型目前并非指针，p的命名待改
	bool m_bMiddlePepUsed; //是否使用中间肽段计算理论碎片离子，默认为false

	int m_nPepLen; //肽段长度，新加字段。设置肽段时重置此值，计算时使用此值作为肽段长度。初始化为0。
	size_t m_nMatchedNo; //2014.5.13 新加字段，记录匹配到谱峰数目
};

}

#endif /*XLINKOPENSCORER_H_*/
