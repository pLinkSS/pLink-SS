#ifndef XLINK_PEP_SALVO_FLOW_H_
#define XLINK_PEP_SALVO_FLOW_H_

//#define _USE_PEP_LIST

#include "../Mass2PepIndex/Mass2PepIndex.h"
#include "../Mass2PepIndex/PeptideReaderFactory.h"
//#include "../XLinkScore/XLinkKSDPScorer.h"
#include "../XLinkScore/XLinkRefineScorer.h"
#include "../XLinkEvaluate/XLinkPeptideEvaluater.h"
#include "PepPairList.h"
#include "../PepPairContainer/common.h"
#include "../PepPairContainer/PepPairContainer.h"
#include "SpectraSearcher.h"
using namespace proteomics_search;
using namespace proteomics_sdk;
using namespace Mass2PepIndex;

#define LINK_SITE_NUM 64

// N-terminal of protein is '[' - 'A' + 26
// C-terminal of protein is ']' - 'A' + 27
// N-terminal of peptide is '(' - 'A' + 28
// C-terminal of peptide is ')' - 'A' + 29

#define PROTEIN_N_CHAR '['
#define PROTEIN_N_SITE 'A' + 26
#define PROTEIN_C_CHAR ']'
#define PROTEIN_C_SITE 'A' + 27
#define PEPTIDE_N_CHAR '('
#define PEPTIDE_N_SITE 'A' + 28
#define PEPTIDE_C_CHAR ')'
#define PEPTIDE_C_SITE 'A' + 29

//class proteomics_sdk::CTrace;

namespace proteomics_search{

class CXLinkPepSalvoFlow:
	public CFlow
{
public:

	CXLinkPepSalvoFlow(vector<CSpectrum> & vSpec);
	virtual ~CXLinkPepSalvoFlow(void);

	virtual void Init(const CCondition &condition, CSearchState * pState, size_t tID);

	virtual void Run(const string & strTmpFile);

	virtual void Close(void);

	virtual string GetVersion(void) const {return "pfind-xlink";};

	virtual CSearchState * GetSearchState(void){return m_pState;}

public:
	void CreatePepPairList();
	
protected:
	/*
	 * generate at least nMinNum of redundant peptide candidates based on a model peptide candidate
	 */
	void _GenerateRedundantPepCandidates(const CXLinkPepResult & modelPepCandidate,vector<CXLinkPepResult> & vPepCandidates, int nMinNum);
	/*
	 * creat peptide pair list for xlink search
	 */
	void _CreatePepPairList();
	
	/*
	 * for debug
	 */
	void _PrintPepPairList();
	/*
	 * creat matrix to record linker sites
	 */
	void _CreateLinkerMap();
	
	void _InitMatchResult();
	/*
	 * create an index for spectra to yield a O(1) query time complexity
	 */
	void _CreateSpecIndex();
	
	/*
	 * preprocess: algorithms for filter invalid spectrum peaks 
	 * or re-calibur precursor masses
	 * score: peptide-spectrum similarity measurement
	 * evaluate: for the match of one spectrum and its top-k scored peptides,
	 * calculate e-value according to all the scored peptides
	 */
	void _Preproc();
	double _Score();
	void _Evaluate(size_t idx);
	
	/*
	 * write temporary file named num.XXXXX.pfd
	 * num: m_tID
	 * XXXXX: an identifier generated randomly to distinguish different experiments
	 */
	void _WritePFDFile(const string & strOutputPath);
	void _ResultCheck();
	
	/*
	 * modifications are stored by their names,
	 * but when generating candidate modified peptides,
	 * it is better to store all modifications according to their sites,
	 * m_mapVarMod and m_mapFixMod are filled by this function,
	 * the key of the two maps is amino acids and 
	 * the value of them are potential modifications on them
	 */
	void _GenMapAA2Mod();
	
	/*
	 * get the linker mass
	 */
	double _GetXLinkerMass(size_t nLinkerId);
	
	/*
	 * check whether two peptides could be x-linked 
	 */
	bool _IsXLinked(const CSimplePeptide &peptide);
	
	void _OpenDatabase();
	
	void _CloseDatabase();
	/*
	 * for debug
	 */
	void _printAssignedPep();

	/*
	 * the following functions is responsible for generating candidate peptides(pairs)
	 * considering modification(s) and cross linker(s) 
	 * for each assigned peptide (pair) result 
	 * and then scoring each one against the corresponding spectra   
	 */
	
	/*
	 * generate all peptides , take into consideration of :
	 * fixed modifications  
	 * variable modifications
	 * cross links , including mono-link loop-link and x-link. 
	 */
	void _GetAllPeptide();
	
	/*
	 * called by _GetAllPeptide to generate each peptide (of each peptide pair) with fixed modification
	 * nPeptideId specify which peptide of a peptide pair is subjected 
	 */
	
	void _AttatchFixedMod(CSimplePeptide *const pPeptide);

	/*
	 * first called by _GetAllPeptide to begin the process of generatiing each peptide (of each peptide pair) with var mod
	 * then called by itself recursively to permute  
	 * nPeptideId specify which peptide of a peptide pair is subjected 
	 * nCurrLen specify the current position (mod site) of the peptide 
	 */
	void _AttatchVarMod(int nPeptideId , int nCurrLen);
	
	/*
	 * called by _AttatchVarMod to attatch the selected linker to the peptide (pair) on each possible site(s)  
	 * when all the modification been attatched  
	 */
	void _AttatchXLink(size_t nLinkerId);
	char _GetCandLinkSiteIndex(const CSimplePeptide &peptide, int nIdx, int nLinkerId);
	bool IsModLinkSiteConflict(const CSimplePeptide &peptide, const int );
	bool IsLinkSiteCleave(const CSimplePeptide &peptide, const int ); //added at 2013.10.30
	
	
	/*
	 * called by _AttatchXLink when a candidate peptide (pair) is complited with mod and link
	 * lookup the spectra index to find proper 
	 * spectra and score each peptide(pair)-spectrum match
	 */
	void _Dispatch(size_t nLinkerId);

	/*
	 * called by _Dispatch
	 * after scoring a peptide(pair)-spectrum match, store the scored peptides into the result
	 * record if its score is high enough
	 *  (usually top-10 peptides are stored for each spectrum)
	 */
	void _track(size_t tIdx, size_t nLinkerId);
	
public:
	
	//members:
	
	/*
	 * static data member to specify whether the peptide pair list has once been created
	 */
	static bool m_bPepPairListCreated;
	
	/*
	 * static data member to store the peptide pairs
	 * the list should be created only once
	 */

#ifdef _USE_PEP_LIST
	static CPepPairList m_CandidatePepPairs;
#else
	static CPepPairContainer m_CandidatePepPairs;
#endif
	/*
	 * assigned peptide (pair)
	 * read in from the peptide pair list 
	 */
	CXLinkPepResult m_AssignedPep;
	
	/*
	 * preprocessing procedure
	 */
	CPreProcess * m_pPreProc;
	
	/*
	 * brand new scorer for xlink
	 */
	//CXLinkKSDPScorer *  m_pScorer;
	CXLinkRefineScorer *  m_pScorer;
	/*
	 * evaluater for xlink 
	 */
	CXLinkPeptideEvaluater * m_pEV;
	
	/*
	 * should be used differently here for xlink
	 */
	CPeptideGenerator* m_pGenerator;
	/*
	 * spectra index
	 */
	vector<SPEC_SORTER_INDEX_INFO> m_vInfo;
	vector<int> m_vSpecWndCnt; /* windows number */
	
	/*
	 * results of each spectra
	 */
	vector<CXLinkMatchResult> m_vResults;
	
	/*
	 * thread id of the outer caller
	 */
	size_t m_tID;
	
	CCondition m_Condition;
	
	map<char, vector<CSimpleMod> > m_mapVarMod;

	map<char, vector<CSimpleMod> > m_mapFixMod;
	
	/*
	 *  save the linked sites in the matrix
	 */
	bool m_bLinked[MAX_LINKER_NUM][LINK_SITE_NUM][LINK_SITE_NUM];

	/*
	 *  save the mono linked sites in the array
	 */
	bool m_bMonoLinked[MAX_LINKER_NUM][LINK_SITE_NUM];
	
	size_t m_tSalvoLine;

	/*
	 * database reader
	 */
	CMass2PepIndexReader * m_pReader;
	
	CSearchState * m_pState;
	
	CSpectraSearcher m_SpecSearcher;
	
	vector<CSpectrum> & m_vSpectra;

	CTrace *m_pTrace;
	string m_strInfo;

	CMapAAMass mapAAMass; //added at 2013.11.21, will be mapped only 1 time in init
};
}
#endif /*PEP_SALVO_FLOW_H_*/
