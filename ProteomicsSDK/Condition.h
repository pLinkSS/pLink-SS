#ifndef CONDITON_H_
#define CONDITON_H_

#include <string>
#include <list>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <time.h>
#include <string.h>
#ifndef IN_VC_PROGRAM_
#include <dirent.h>
#endif

using namespace std;

namespace proteomics_sdk {
struct PEP_TOL_WND {
	PEP_TOL_WND() {
		m_lfPepTol = 0.0;
		m_strPepTolType = "Da";
		m_lfPepTolBase = 0.0;
		m_strPepTolBaseType = "Da";
	}

	double m_lfPepTol;
	string m_strPepTolType;
	double m_lfPepTolBase;
	string m_strPepTolBaseType;
};

class CCondition {
public:
	CCondition();
public:
	/* database */
	string m_strXLinkerPath; /* xlink_List xlink.ini */
	IndexContent m_eIndexContent; /* index_conent */
	IndexType m_eIndex; /* index_type */
	vector<string> m_vAllDataBaseName; /* all database */
	vector<string> m_vSelectedDBName;

	/* enzymes */
	vector<CEnzyme> m_vAllEnzyme;
	CEnzyme m_SelectedEnzyme;

	/* modifications */
	vector<CModification> m_vAllModification;
	vector<CModification> m_vSelectedFixMod;
	vector<CModification> m_vSelectedVarMod;
	int m_nMaxModifyNumber;

	/* xlinkers */
	vector<CXLinker> m_vAllXLinker;
	vector<CXLinker> m_vSelectedXLinker; /* light & heavy linker */

	/* ions */
	bool m_bPepMono; /* for peptide: use mono mass or average mass */
	bool m_bFragmentMono; /* for fragment: use mono mass or average mass */
	vector<PEP_TOL_WND> m_vPepTolWnds; /* multiple mass tolerance windows */
	double m_lfFragmentTol;
	string m_strFragmentTolType; /* e.g. ppm or Da */
	double m_lfFragmentTolBase;
	string m_strFragmentTolBaseType; /* e.g. ppm or Da */
	double m_lfMinPepTol;
	string m_strMinPepTolType;
	vector<CIonType> m_vIonTypes; /* ion type vector */

	/* simple score: for xlink open flow */
	PreProcType m_eSimplePreProcMethod;
	double m_lfSimpleScoreCutoff;
	size_t m_nSimpleReportPep; /* the output number of peptide */
	vector<CIonType> m_vSimpleIonTypes;

	/* flow: for refinement score */
	size_t m_tProcessorNum;
	size_t m_tMinScoreNum;
	size_t m_tMaxScoreNum;
	size_t m_tSalvoBatchSize;
	LogRankType m_eLogRank; /* log rank */
	ModuleType m_eModule;
	LogAppenderType m_eLogAppender;
	SearchFlowType m_eFlow;
	InstrumentType m_eInstrumentType;
	PreProcType m_ePreProcMethod;
	ScoreType m_eScoreMethod;
	EvaluateType m_eEV;
	double m_lfMaxEV; /* max E-value */
	size_t m_nReportPep; /* the number of peptides to report */
	size_t m_nReportPro; /* the number of proteins to report */
	ProteinReportType m_eProReport;
	double m_lfFPRThreshold; /* false positive rate */
	string m_strFPSign; /* e.g. REVERSE */
	ProteinInferType m_eInferType; /* ????? */
	size_t m_tInterContinuousWndNum; /* ????? */

	/* spectrum */
	string m_strSpectraTitle;
	/*Spec path and output path 20140525*/
	string m_strSpectraPath;
	string m_strOutputPath;

	/* cluster */
	string m_strSpectraIndexPath; /* block_path */
	SpectraIndexType m_eSpectraIndex; /* block_format_type */
	LoadBalanceType m_eLoadBalance; /* load_balance_type */
	int m_nMinBlockNum;
	int m_nMaxBlockSize;

	size_t m_nCleaveWay; /* add by czhou */
	bool m_bPep2Pro;
	size_t m_nMaxMissCleaves;
	string m_strDBConfPath;
	string m_strEnzymeListPath; /* enzyme.ini */
	string m_strModifyPath; /* modify.ini */
	string m_strAAListPath; /* aa.ini */
	string m_strIndexPath; // zhangkun
	CMapAAMass m_mapAA;

	/* for peptides generation */
	PepGeneType m_ePepGen;
	int m_nKSDP_L;
	double m_lfKSDP_G;
	double m_lfKSDP_A;

	/* refined-search: abandoned  */
	bool m_bRSSaveInfo;
	bool m_bRSLoadInfo;
	string m_strRSPath;
	double m_lfRSFDR;
	string m_strRSPindexPath;
	string m_strRSDBName;
	string m_strRSTag;

	/* for sumo */
	string m_strSumoSeq;
	char m_strSumoSite;// For the changing pep, usually K
	char m_strFixSumoSite;//For the fixed pep, usually K

	/* for triple peptide */
	bool m_bSingleCAsSitePep;
	string m_strReverseTag;

	bool ValidateElement(void) const;
	bool ValidateAll(void) const;
	void clear(void);

	void SetDefaultIons(void);

	double GetExpMass(double lfMH, int nCharge) const;
	double GetMZ(double lfMH, int nCharge) const;
	double GetNeutraMass(double lfExpMass) const;
	double GetPepTolBase(double lfMZ, int nCharge) const;
	double GetPepTolBase(int nWnd, double lfMZ, int nCharge) const;
	double GetPepTol(double lfMZ, int nCharge) const;
	double GetPepTol(int nWnd, double lfMZ, int nCharge) const;
	double GetMinPepTol(double lfMZ, int nCharge) const;
	void GetMHMassBorder(double lfMH, int nCharge, double &lfLower,
			double &lfUpper) const;
	void GetMHMassBorder(int nWnd, double lfMH, int nCharge, double &lfLower,
			double &lfUpper) const;
	void GetMassBorder(double lfMH, int nCharge, double &lfLower,
			double &lfUpper) const;
	void GetMassBorder(int nWnd, double lfMH, int nCharge, double &lfLower,
			double &lfUpper) const;
	double GetLeastNegativeMod() const;
	double GetMaxMassDif() const;
	string GetInstrumentType(void);

	void GetIonMassBorder(double lfMz, int nCharge, double &lfLower,
				double &lfUpper) const;
}; /* end of CCondition */
} /* end of namespace proteomics_sdk */

#endif /* CONDITON_H_ */
