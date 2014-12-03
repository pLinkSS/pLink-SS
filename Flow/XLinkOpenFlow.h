#ifndef XLINKOPENFLOW_H_
#define XLINKOPENFLOW_H_

#include "../Mass2PepIndex/Mass2PepIndex.h"
#include "../Mass2PepIndex/PeptideReaderFactory.h"
#include "../XLinkScore/XLinkOpenScorer.h"
#include "../XLinkScore/XLinkRefineScorer.h"
#include "common.h"
using namespace proteomics_search;
using namespace proteomics_sdk;
using namespace Mass2PepIndex;

#define LINK_SITE_NUM 64

/* different kinds of linker sites, all amino acid is * */
#define PROTEIN_N_CHAR '[' /* N-terminal of protein symbol */
#define PROTEIN_N_SITE 'A' + 26 /* N-terminal of protein index */
#define PROTEIN_C_CHAR ']' /* C-terminal of protein symbol */
#define PROTEIN_C_SITE 'A' + 27 /* C-terminal of protein index */
#define PEPTIDE_N_CHAR '(' /* N-terminal of peptide symbol */
#define PEPTIDE_N_SITE 'A' + 28 /* N-terminal of peptide index */
#define PEPTIDE_C_CHAR ')' /* C-terminal of peptide symbol */
#define PEPTIDE_C_SITE 'A' + 29 /* C-terminal of peptide index */
#define MIN_PEPTIDE_MASS 300

namespace proteomics_search{

class CXLinkOpenFlow:
	public CFlow
{
public:
	CXLinkOpenFlow(vector<CSpectrum> & vSpec);
	virtual ~CXLinkOpenFlow();
	
	virtual void Init(const CCondition &condition, CSearchState * pState, size_t tID);

	virtual void Run(const string & strTmpFile);

	virtual void Close(void);

	virtual string GetVersion(void) const {return "pfind-xlink-openflow";};

	virtual CSearchState * GetSearchState(void){return m_pState;}
	
protected:
	void _InitMatchResult();
	
	void _CreateSpecIndex();
	
	void _CreateLinkerMap();
	
	bool _IsXLinked(const CSimplePeptide &peptide);
	
	void _Preproc();
	
	double _Score();
	
//	int _GetLinkSiteIndex(CXLinkOpenPepResult & pep_res1);

//	bool _CheckLinkSiteLegal(CXLinkOpenPepResult & pep_res1, CXLinkOpenPepResult & pep_res2);

	void _GotoRefineStage(size_t idx);
	
	void _RefineDispatch(size_t idx);
	
	void _RefineTrack(size_t tIdx);
	
	void _GenerateRedundantPepCandidates(const CXLinkPepResult & modelPepCandidate,vector<CXLinkPepResult> & vPepCandidates, int nMinNum);
	
	void _Evaluate(size_t idx);
	
	void _GenMapAA2Mod();
	
	void _OpenDatabase();
	
	void _CloseDatabase();

	void _GetAllPeptide();
	
	void _AttatchFixedMod();

	void _AttatchVarMod(int nCurrLen);
	
	char _GetCandLinkSiteIndex(const CSimplePeptide &peptide,
			int nIdx, int nLinkerId);

	bool IsModLinkSiteConflict(const CSimplePeptide &peptide,
			const int nIdx);
	bool IsLinkSiteCleave(const CSimplePeptide &peptide, const int ); //added at 2013.10.30

	void _Dispatch();

	void _track(size_t tIdx,bool bLowOrHigh);
	
	double _Calc_Theoretical_MH();
	
	double _Calc_Theoretical_MH(const CSimplePeptide & pep_res);
	
	double _GetXLinkerMass(bool bMono,int nLinkerId);
	
	double _GetWaterMass();
	
	double _GetProtonMass();
	
	void _ResultCheck();
	void _WritePFDFile(const string & strOutputPath);
	
	void _CheckResult(const string & strOutputPath);
	

public:
	vector<CSpectrum> & m_vSpectra; /* holding spectra */
	CCondition m_Condition; /* holding parameters */
	CSearchState * m_pState;
	size_t m_tID; /* seems no use, always be zero, refer to DynamicSlave */
	
	/* different stages objects */
	CPreProcess * m_pPreProc;
	CPreProcess * m_pRefinePreProc;
	CXLinkOpenScorer *  m_pScorer;
	CXLinkRefineScorer *  m_pRefineScorer;
	CXLinkPeptideEvaluater * m_pEV;
	CPeptideGenerator* m_pGenerator;

	/* modifications */
	map<char, vector<CSimpleMod> > m_mapVarMod;
	map<char, vector<CSimpleMod> > m_mapFixMod;
	
	/* linkers */
	bool m_bMonoLinked[MAX_LINKER_NUM][LINK_SITE_NUM];
	bool m_bLinked[MAX_LINKER_NUM][LINK_SITE_NUM][LINK_SITE_NUM];
	
	/* spectra index information */
	vector<SPEC_WND_INFO> m_vSpectraIndex; /* windows boundary */
	vector<int> m_vSpecWndCnt; /* windows number */

	/* holding open results */
	vector<vector<CXLinkOpenPepResult> > m_vLowMassResults;
	vector<vector<CXLinkOpenPepResult> > m_vHighMassResults;

	/* holding refinement results */
	vector<CXLinkMatchResult> m_vResults;
	vector<size_t> m_vResultsCandidateNum;

	CXLinkOpenPepResult m_AssignedPep;
	CXLinkPepResult m_AssignedXLinkPep;

	size_t m_tSalvoLine;
	size_t m_tCurrentSpec;
	
	/* for database access */
	CMass2PepIndexReader *m_pReader;

	/* for log and debug */
	CTrace *m_pTrace;
	string m_strInfo;
};

}

#endif /* XLINKOPENFLOW_H_ */
