#ifndef CXLinkRefineScorer_H_
#define CXLinkRefineScorer_H_

#include "CommonClass.h"

typedef unsigned char UCHAR;
namespace proteomics_sdk
{

class CXLinkRefineScorer
{
public:
	CXLinkRefineScorer():m_pPeptide(NULL),m_pSpectrum(NULL),m_pbMatched(NULL),/*m_pPeakWeight(NULL),*/m_nMaxPeaksNum(5000), m_ionmz(NULL){};
	virtual ~CXLinkRefineScorer(){Close();};
	
public:
	virtual void Initialize(const CCondition & condition);
	virtual void Close(void);

	virtual void SetSpectrum(CSpectrum & spec);
	virtual void SetPeptide(CXLinkPepResult & pep);
	
	virtual double Score();
	
	void _Match();
	
	void _ComputeMZ();
	
	void _GetMassBorder(int nMz, int nChg, int &nMin, int &nMax);
	
	void _SetCondition(const CCondition & cond);
	
	void _CalMatchInfo();
	double i2Score();
	double i2HyperGeoMetric();
	double CalcuCombine(int a, int b);

	void _PrintMZ();//added at 2013.11.21

protected:

	//members:
	CCondition m_Condition;
	CXLinkPepResult * m_pPeptide;
	CSpectrum * m_pSpectrum;
	size_t m_tionmzSize;
	
	// for inter-continous window
	int m_InterContiWnd[MAX_IONTYPE_NUM][2*MAX_PEPTIDE_LENGTH];
	
	int m_FullInterContiWnd[MAX_IONTYPE_NUM];
	
	// for intra-continous window
	int m_iontag_pep[MAX_IONTYPE_NUM][2*MAX_PEPTIDE_LENGTH];
	
	//double m_lfMatchedIons1[MAX_IONTYPE_NUM][2*MAX_PEPTIDE_LENGTH];
	//double m_lfMatchedIons2[MAX_IONTYPE_NUM][2*MAX_PEPTIDE_LENGTH];
	
	//double * m_pPeakWeight;
	
	char * m_pbMatched;
	size_t m_nMaxPeaksNum;
	vector<CMzTriple> m_ionmz;
	
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
	
	//some constants
	double alfa;
	double gamma;
	double exp_table[2 * 2 * MAX_PEPTIDE_LENGTH];
	int l_win;
	// left and right of correlative window
	int l1,l2; 
	
	//FILE * m_fp;
};
}

#endif /*CXLinkRefineScorer_H_*/
