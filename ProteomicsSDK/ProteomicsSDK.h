#ifndef PROTEOMICS_SDK_H_
#define PROTEOMICS_SDK_H_

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
#include <dirent.h>

using namespace std;


#define MAX_PEP_TOL_WND_NUM 10
#define MAX_LINKER_NUM 3

namespace proteomics_sdk {

const double IonMass_Proton = 1.00727647012;
const double IonMass_Mono_H = 1.007825035;
const double IonMass_Aver_H = 1.00794;
const double IonMass_Mono_C = 12;
const double IonMass_Aver_C = 12.0107;
const double IonMass_Mono_N = 14.003074;
const double IonMass_Aver_N = 14.0067;
const double IonMass_Mono_O = 15.99491463;
const double IonMass_Aver_O = 15.9994;
const double IonMass_Mono_Na = 22.98976967;
const double IonMass_Aver_Na = 22.98976967;
const double IonMass_Mono_K = 38.9637069;
const double IonMass_Aver_K = 39.0983;
const double IonMass_Mono_Ca = 39.9625912;
const double IonMass_Aver_Ca = 40.0784;
const double NH3Mass_Mono = (IonMass_Mono_N + IonMass_Mono_H * 3);
const double NH3Mass_Aver = (IonMass_Aver_N + IonMass_Aver_H * 3);
const double H2OMass_Mono = (IonMass_Mono_O + IonMass_Mono_H * 2);
const double H2OMass_Aver = (IonMass_Aver_O + IonMass_Aver_H * 2);


typedef struct _PEP_SQ
{
	double			dfMass;
	string			strSQ;
	unsigned char	cMiss;
	size_t			tInvertedFilesPos;
	unsigned char	cDatNum;
	unsigned char	cEnd;//0 means nothing, 1 means left, 2 means right, 3 means left&right
} PEP_SQ;


// add by emily 
// guarantee the correctness of path 
class CCheck {
public:

public:
	static void CheckPath(string & strPath);
};

class CErrInfo {
	string m_strClass;
	string m_strMethod;
	string m_strDetail;
	string m_strInfo;
	string m_strException;
public:
	CErrInfo(const string &strClass, const string &strMethod,
			const string &strDetail = "");

	CErrInfo(const string &strClass, const string &strMethod,
			const string &strDetail, const exception & e);
	void Append(const string &strInfo);

	string Get() const;

	string Get(const exception& e);

	void SetException(const exception & e);

	friend ofstream& operator<<(ofstream& os, const CErrInfo& info);
	friend ostream& operator<<(ostream& os, const CErrInfo& info);

};

//to denote an amino acid
class CAA {
public:
	//single letter abbreviation for the amino acid
	char m_cAA;
	//three letter abbreviation for the amino acid
	string m_strAA;
	//the description of the amino acid
	string m_strAADetail;
	//average mass
	double m_lfAvrgMass;
	//monoisotopic mass
	double m_lfMonoMass;
	//mass of the modification,depending on current bMono .
	double m_lfMod;
	//Average mass(integer)
	size_t m_tavrgMass;
	//Monoisotopic mass(integer)
	size_t m_tmonoMass;

	CAA(const CAA & r);

	CAA();

	CAA(char cAA, string strAA, string strAADe, double lfAvg, double lfMono,
			double lfMod = 0);

	virtual ~CAA(void);

	void InitByString(string& strInit);
};

class CMapAAMass {
public:

	CMapAAMass(const CMapAAMass & m);

	CMapAAMass();
	void SetAAMap(char cAA, CAA& aa);
	CAA& GetAA(char cAA);

	//members
	map<char, CAA> m_mapAAMass;
};

// pfind-xlink
class CXLinker {
public:
	// name for the linker
	string m_strName;

	// AA specify for site1 
	string m_strAlphaAA;

	// AA specify for site2
	string m_strBetaAA;

	// homobifunctional linker or heterobifunctional linker 
	// if m_strAlphaAA == m_strBetaAA then m_bHomo = true else m_bHomo = false
	bool m_bHomo;

	// x-linker mass shift (mono)
	double m_lfMonoMass_dif;

	// x-linker mass shift (avrg)
	double m_lfAvrgMass_dif;

	// mono-linker mass shift (mono)
	double m_lfMLMonoMass_dif;

	// mono-linker mass shift (avrg)
	double m_lfMLAvrgMass_dif;

	CXLinker();

	~CXLinker();

	void InitByString(string & strInit);

};

class CModification {
public:
	CModification();

	~CModification();

	void InitByString(string& strInit);

	ModType m_eModType;
	//  mass added to the  residue , in mono mass
	double m_lfMonoMass_dif;
	//  mass added to the  residue , in average mass
	double m_lfAvrgMass_dif;

	vector<double> m_vlfAvrgNeutralLoss_dif;

	vector<double> m_vlfMonoNeutralLoss_dif;

	// amino acid (AA) that to modify  mass
	string m_strAA;

	// the modification name
	string m_strName;

	// the number of neutral losses
	size_t m_tNLSize;
};


class CPeak {
public:
	CPeak(double mz, double intensity, int nmz);

	CPeak();

	CPeak& operator=(const CPeak& rhs);

	//members
public:
	double lfMz;
	double lfIntensity;
	int nMz;

	// add by emily
	int nCharge;
};
class CEnzyme {
public:
	CEnzyme();
	virtual ~CEnzyme();

	CEnzyme& operator=(const CEnzyme& rhs);

	bool GetIsNTerm() const;
	string GetCleaveString() const;
	string GetName() const;
	string GetNotCleave() const;

	void SetIsNTerm(bool bNTerm);
	void SetCleaveString(const char * szCleave);
	void SetName(const char * szName);
	void SetNotCleave(const char * szNoClv);

	bool empty(void) const;
	void InitByString(string& strInit);

	// NTerm or CTerm
	bool m_bNTerm;
	// Enzyme name
	string m_strName;
	// Cleave amino acids
	string m_strCleave;
	// Exception amino acids
	string m_strNotCleave;
};

class CSimplePeptide {
public:

	CSimplePeptide();
	CSimplePeptide(const CSimplePeptide & peptide);

	CSimplePeptide &operator=(const CSimplePeptide& peptide);
	void SetPeptideInfor(const char * pszSequence, size_t tLength,
			double lfMass, unsigned char cMissCleave, unsigned char cEnd,
			bool bCopyNCTerm);
	bool EqualPeptide(const CSimplePeptide & peptide);
	bool IsProNTerm();
	bool IsProNTerm() const;
	bool IsProCTerm();
	bool IsProCTerm() const;
	static double Calc_Theoretical_MH(const CSimplePeptide & pep,
			bool bPepMono);
	static double Calc_Theoretical_M(const CSimplePeptide & pep); //added at 2014.5.4, calculate the M (+H2O)
	//members:
	char m_szSequence[MAX_PEPTIDE_LENGTH + 1];

	size_t m_tLength;

	double m_lfMass;

	size_t m_tModCnt;
	size_t m_tFixedModCnt;

	unsigned char m_cMissCleave;

	size_t m_tModSites[MAX_MODIFY_NUM][2];

	char m_cPrev;
	char m_cNext;

	unsigned char m_cEnd;
};

// pfind-xlink
class CXLinkItem {
public:
	CXLinkItem() {
		m_eXLinkType = (XLinkType) 0;
		m_tAlphaSite = -1;
		m_tBetaSite = -1;
		//todo :
		m_nLinkerId = 0;
	}
	CXLinkItem & operator=(const CXLinkItem &);
	void setXLinkItemInfo(const XLinkType &eXLinkType,
			const int nAlphaSite=-1, const int nBetaSite=-1,
			const int nLinkerId=-1);
	XLinkType m_eXLinkType;
	int m_tAlphaSite;
	int m_tBetaSite;
	// various linker
	// light & heavy linker
	int m_nLinkerId;
};

struct CMatchInfo {
	int aPepConf[2 * MAX_PEPTIDE_LENGTH];
	double lfMatchedSpecInt;
	double lfUnMatchedSpecInt;
};


// pfind-xlink
class CXLinkPepResult {
public:
	CXLinkPepResult();

	CMatchInfo m_stMatchInfo;
	bool m_bPair;
	CXLinkItem m_XLink;

	CSimplePeptide m_AlphaPeptide;
	CSimplePeptide m_BetaPeptide;

	static bool Score_Greater(const CXLinkPepResult & a,
			const CXLinkPepResult & b);
	void SwapAlphaBeta();
	CXLinkPepResult & operator=(const CXLinkPepResult &);

	double m_lfScore;

	double m_lfCalc_MH;

	double m_lfEvalue;

	double m_lfAlphaEvalue; // added at 2014.9.10 E-value for single peptide
	double m_lfBetaEvalue;

	bool m_bEV;

	double m_lfFPR;

	vector<string> m_vAlphaProteinAC;
	vector<size_t> m_vAlphaProteinID;
	vector<string> m_vBetaProteinAC;
	vector<size_t> m_vBetaProteinID;

	// site of the first AA in of peptide in the protein sequence
	vector<size_t> m_vAlphaProteinSite;
	vector<size_t> m_vBetaProteinSite;
};

// add by xlink
class CXLinkOpenPepResult {
public:

	CSimplePeptide m_peptide;

	double m_lfScore;

	int m_nLinkSite;

	double m_lfOpenMass;

	int m_nLinkerId;
	double m_lfOtherMinMass;
	double m_lfOtherMaxMass;
	CXLinkOpenPepResult();
	void resetPeptideInfo(string strSq, double lfMass,
			unsigned char cMissCleave, unsigned char cEnd, bool bCopyNCTerm);
};

//add at 2014.5.11, to store link and mod info for middle peptide
class CXLinkMiddleModAndLink
{
public:
	size_t m_nSidePepSite : 7; //与第一条边肽的连接位点。肽段最长100，7位足够
	size_t m_nLinkSite : 7; //与另一条未知边肽的连接位点。

	size_t m_tModCnt : 5; //修饰数目。最大20，5位足够
	size_t m_tFixedModCnt : 5; //固定修饰数目

	vector <size_t> m_tModPos; //记录修饰位置的vector
	vector <size_t> m_tModIdx; //记录修饰Index的vector

	size_t m_nPepIdx : 5; //这个修饰对应的中间肽段序列的Index，即这张谱图的哪条中间肽段。最大20，5位足够。
	size_t m_nSidePepIdx : 5; //这个修饰对应的侧边肽段序列的Index，即这张谱图的哪条边肽段。最大20，5位足够。

	//Notice: LinkerID信息没有记录

	float m_fScore; //记录当前中间肽段的打分，由于后面暂未用，先不加入 //2014.8.26 加入，用于边肽-中间肽组合筛选，由400名改为保留更少。注意内存占用变大，后续不加粗打分筛选时要去掉此属性。
//	double m_lfOpenMass; //剩余质量，后面重算，不必加入
};

//add at 2014.4.21, for middle peptide
class CXLinkOpenMiddlePepResult {
public:

	CSimplePeptide m_peptide;

	double m_lfScore;

	int m_nLinkSite;

	double m_lfOpenMass;

	double m_lfSidePepMass; //The known side pep mass (this side pep has already been known, will find another)
	int m_nSidePepSite; //The pep site, where the known side pep is on.

	int m_nLinkerId;
	double m_lfOtherMinMass;
	double m_lfOtherMaxMass;
	CXLinkOpenMiddlePepResult();
	void resetPeptideInfo(string strSq, double lfMass,
			unsigned char cMissCleave, unsigned char cEnd, bool bCopyNCTerm);

	void SetPep(CXLinkOpenPepResult pep); //用CXLinkOpenPepResult设置本类型

	CXLinkMiddleModAndLink GetModAndLink(); //从本类型中获取CXLinkMiddleModAndLink类型
	void SetModAndLink(const CXLinkMiddleModAndLink & pep); //用CXLinkMiddleModAndLink设置本类型的修饰和交联位点

};



class CPeptideResult {
public:

	CPeptideResult();
	static bool
	Score_Greater(const CPeptideResult & a, const CPeptideResult & b);

	CSimplePeptide m_peptide;

	double m_lfScore;

	//double s, t;

	double m_lfEvalue;

	double m_lfCalc_MH;

	bool m_bEV;

	double m_lfFPR;

	vector<string> m_vProteinAC;
	vector<size_t> m_vProteinID;
};

/*
 * refined-search needed
 * store evalue coefficients
 */

struct EVCoefficient {
	double lfCoef0;
	double lfCoef1;

	EVCoefficient();

	int GetSize();

	bool WriteToFile(FILE * fp);

	bool ReadFromFile(FILE * fp);

};

/*
 enum ActivationType
 {
 AT_CID_IT = 0,
 AT_CID_FT = 1,
 AT_ETD_IT = 2,
 AT_ETD_FT = 3
 }
 */

class CMS1Info {
public:
	size_t m_tScanNo;
	double m_lfRetentionTime;
	string m_strInstrumentType;

};

class CSpectrum {
public:
	CSpectrum();
	virtual ~CSpectrum();

	void CopyBasicItems(const CSpectrum & spec);
	CSpectrum & operator =(const CSpectrum & spec);
	CSpectrum(const CSpectrum & spec);

	static bool intesity_greater(const CPeak & elem1, const CPeak & elem2);
	static bool intesity_lesser(const CPeak & elem1, const CPeak & elem2);
	static bool mz_lesser(const CPeak & elem1, const CPeak & elem2);
	void clear(void);
	void ReAlloc(size_t tNum);

	// refined-search needed
	// store Evalue coefficient of the first round of search

	struct EVCoefficient m_stEVCoef;

	int m_nEvalueNo;

	int m_nCharge;

	double m_lfIntensity;

	//the sqrt of the highest peak's intensity
	double m_lfSqtMaxInten;

	double m_lfMH;

	double m_lfMZ;

	string m_strFilePath;

	CPeak * m_pPeaks;
	size_t m_tPeaksNum;

	size_t m_tScanNo;
	double m_lfRetentionTime;
	string m_strActivationType;

	size_t m_tPrecursorScanNo;
	string m_strInstrumentType;

	//double s, t;
	vector<int> m_vHash;
	void Create_Hash_Table();

};


class CXLinkMatchResult {
public:
	size_t m_tScore;
	vector<double> m_vlfScores;
	size_t m_tRealCandidate;
	vector<CXLinkPepResult> m_vPeptideResults;

	CXLinkMatchResult();

	void remove_invalid(void);

	static double Calc_Theoretical_MH(const CXLinkPepResult & pep_res,
			bool bPepMono);

	CXLinkMatchResult & operator=(const CXLinkMatchResult &);
};

class CSimpleMod {
public:
	size_t m_nIdx;
	double m_lfMass;
	ModType m_eModType;
};

class CProtein {
public:

	size_t m_nID;

	//ACCESSION ID
	string m_strAC;

	//Description
	string m_strDE;

	//Sequence
	string m_strSQ;

	CProtein();
	virtual ~CProtein();

	CProtein(string strAC, string strDE, string strSQ);

	CProtein(const CProtein & pro);
	CProtein Reverse(void);
};

class CBaseProtein: public CProtein {
public:

	string m_strID;

	string m_strDT;

	string m_strGN;

	string m_strOS;

	string m_strOG;

	string m_strOC;

	string m_strOX;

	string m_strOH;

	string m_strRN;

	string m_strRP;

	string m_strRC;

	string m_strRX;

	string m_strRG;

	string m_strRA;

	string m_strRT;

	string m_strRL;

	string m_strCC;

	string m_strDR;

	string m_strPE;

	string m_strKW;

	string m_strFT;

	string m_strSQHeader;

	CBaseProtein();
	virtual ~CBaseProtein();
};

class CCleaver {
public:
	void Cleave_specific(const string & strSQ);
	void SetEnzyme(const CEnzyme & enzyme);

	void InitValue(double pMass[], const CEnzyme & enzyme);

public:
	vector<int> m_viCleaves;
	vector<double> m_vlfMasses;

protected:
	CEnzyme m_enzyme;

	double m_lfMass[26];
};
class CAssignedProtein {
public:

	string m_strAC;

	string m_strDE;

	size_t m_nProID;

	vector<string> m_vSameSet;

	vector<size_t> m_vSameSetID;

	vector<string> m_vSubSet;

	vector<size_t> m_vSubSetID;

	int m_nPeptide;

	vector<pair<const CPeptideResult *, size_t> > m_vpContainPep;

	void EmptyContents();

	double GetTotalScore(void) const;

	static bool Total_Score_Greater(const CAssignedProtein& elem1,
			CAssignedProtein& elem2);
	double GetTotalEValue(void) const;

	static bool Total_EValue_Lesser(const CAssignedProtein& elem1,
			const CAssignedProtein& elem2);

	double GetTotalFPR(void) const;

	double GetPepAverRank(void) const;
};

class CFileFind {
public:
	bool Open(const string & strFilePath);
	bool GetNextFile();
	void Close();
	//members:
	DIR * m_pDir;
	struct dirent * m_pent;
};

class CThreadState {
public:
	size_t m_tID;
	int m_nCurrentSpectra;
	int m_nTotalSpectra;
	// add for xlink output
	float m_lfCurrentDbRatio;
	float m_lfTotalDbRatio;
};

class CSearchState {
public:
	CSearchState();

	//add for xlink output;
	void XlinkConsole();
	void GetProgress(string & strInfo);

	void SetTotalSpectra(int nTotal);
	void SetCurrentSpectra(int nCurrent);
	void SetThreadTotal(size_t tID, int nTotal);
	void SetThreadCurrent(size_t tID, int nCurrent);
	// set for xlink output
	void XlinkSetThreadCurrent(size_t tID, float lfCur, float lfTotal);

	bool m_bSet;
	//string m_strCondition;
	double m_lfStartTime;
	time_t m_tStartTime;
	int m_nTotalSpectra;
	int m_nCurrentSpectra;
	vector<CThreadState> m_vThreadState;

};

}

#endif /*PROTEOMICS_SDK_H_*/
