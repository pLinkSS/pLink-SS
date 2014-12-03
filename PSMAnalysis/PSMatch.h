#ifndef PSMATCH_H_
#define PSMATCH_H_

struct PeakInfo
{
	int nPeakType;
	// 0: common
	// 1: isotop
	// 2: -NH3
	// 3: -H2O
	// 4: MH
	// 5: noise
	int nMainPeakId;
	double lfCalDeltaMz;
	double lfExpDeltaMz;
	int nCharge;
	
	PeakInfo()
	{
		nPeakType = 0;
		nMainPeakId = -1;
		lfCalDeltaMz = 0.0;
		lfExpDeltaMz = 0.0;
		nCharge = 0;
	}
};

struct CMzTriple
{
	int nMz;
	int nIonTypeOrder;
	// for internal ions , the order is assigned as the sequence of peptide
	int nPepPosOrder1;
	// for internal ions;
	int nPepPosOrder2;
	// for internal ion;
	bool bNTerm1 ;
	bool bNTerm2 ;
	
	bool bContainLinker;
	bool bSameSite;
	bool bAddToTag;
	
	int nAAnum;
	static bool Less(const CMzTriple & a, const CMzTriple & b)
	{
		return a.nMz < b.nMz;
	}
	
	void clear()
	{
		nMz = 0;
		nIonTypeOrder = -1;
		nPepPosOrder1 = -1;
		nPepPosOrder2 = -1;
		bNTerm1 = false;
		bNTerm2 = false;
		bContainLinker = false;
		bSameSite = false;
		bAddToTag = false;
		nAAnum = 0;
	}
};

struct IonMatchInfoEx
{
	int nIonTypeOrder;
	// index in the parameter setting
	int nIonType;
	// 0 : N-term ions
	// 1 : C-term ions
	// 2 : Internal ions
	// 3 : Precursor ions (one peptide) 
	// 4 : Precursor ions (two peptide)
	// 5 : Reporter ions
	double lfCalMz;
	double lfExpMz;
	double lfTol;
	double lfIntensity;
	int nRank;
	bool bContainLinker;
	int nCharge;
	int nLossNH3;
	int nLossH2O;
	int nPepPosOrder1;
	int nPepPosOrder2;
	bool bSameSite;
	double lfOtherLoss;
	double lfExpTotalLoss;
	int nAAnum;
	
	void clear()
	{
		nIonTypeOrder = -1;
		nIonType = -1;
		lfCalMz = 0.0;
		lfExpMz = 0.0;
		lfTol = 0.0;
		lfIntensity = 0.0;
		nRank = -1;
		bContainLinker = false;
		nCharge = -1;
		nLossNH3 = 0;
		nLossH2O = 0;
		nPepPosOrder1 = -1;
		nPepPosOrder2 = -1;
		bSameSite = false;
		lfOtherLoss = 0.0;
		lfExpTotalLoss = 0.0;
		nAAnum = 0;
	}
};


struct IonMatchInfo
{
	int nIonTypeOrder;
	int nMatchIonCount;
	int nTheoIonCount;
	// 匹配离子数目/总匹配离子数目
	double lfMatchCountRatio;
	// 匹配离子数目/总理论离子数目
	double lfMatchGainRatio;
	// 平均匹配强度
	double lfMatchAvrgItensity;
	// 平均匹配误差
	double lfMatchAvrgTolerance;
	// 谱峰连续性 = 连续离子/总离子数目
	// 目前只对cType = 0和1的离子，即肽段主干一次碎裂的离子有效
	double lfMatchContiRatio;
	// AA num: 碎片离子里所含的氨基酸数目
	double lfAvrgAAnum;
	
	
	// 重要性度量
	double lfSignificance;
	
	string strMaxItensityFileName;
	double lfMaxItensity;
	
	int nSpectraCnt ;
	void clear()
	{
		nIonTypeOrder = 0;
		nMatchIonCount = 0;
		nTheoIonCount = 0;
		lfMatchCountRatio = 0.0;
		lfMatchGainRatio = 0.0;
		lfMatchAvrgItensity = 0.0;
		lfMatchAvrgTolerance = 0.0;
		lfMatchContiRatio = 0.0;
		lfAvrgAAnum = 0.0;
		lfSignificance = 0.0;
		lfMaxItensity = 0.0;
		strMaxItensityFileName = "";
		nSpectraCnt = 0;
	}
};

const int CO_H = int((IonMass_Mono_C + IonMass_Mono_O - IonMass_Mono_H) * MZMULTIPLIER);
const int nhmass_mono_multi = (int)(IonMass_Mono_H*MZMULTIPLIER);
const int nomass_mono_multi = (int)(IonMass_Mono_O*MZMULTIPLIER);
const int nhmass_mono_multi_3 = (int)(IonMass_Mono_H*MZMULTIPLIER*3);

const int nhmass_avrg_multi = (int)(IonMass_Aver_H*MZMULTIPLIER);
const int nomass_avrg_multi = (int)(IonMass_Aver_O*MZMULTIPLIER);
const int npmass_multi = (int)(IonMass_Proton*MZMULTIPLIER);

class CPSMatch
{
public:
	CPSMatch();
	virtual ~CPSMatch();
	
public:
	void Initialize(const CPSMConf & psmconf,FILE * fp);
	void Close(void);
	void SetSpectrum(CSpectrum & spec);
	void SetPeptide(CXLinkPepResult & pep);
	
	void Preprocess();
	void ComputeMZ();
	void Match();
	void Output();
	
	void GetPSMMatchInfo();
	double GetMatchOdd(int nPepId);
	
	inline void _GetMassBorder(int nMz, int nChg, int &nMin, int &nMax);
	void SetCondition(const CPSMConf & cond);
	string GetIonDescription(int IonId);
	
	void _SetPeakRank();
	double _CalProb(int N, int n, int n0, int l, int x, int x1);
	double _CalProb(int N, int n, int l, int x);
	double _C(int D,int U);
	void _CalMatchInfo();
	bool _HasIsotope(size_t tSpecId, size_t tIonId);
	double _GetIonTol(size_t tSpecId,size_t tIonId);
	double GetMaxTagLength(int nPep);
	
	static void _InitLogNPerm();
	static void _InitLogTagLenCDF();
	
	bool _IsMatched(double lfCalMz,double lfExpMz);
protected:
	double m_lfPepMass1;
	double m_lfPepMass2;
 
protected:
	
	static double m_lfLogNPerm[100000];
	static double m_lfLogStdNormCDF[101];
	static double m_lfLogGamaCDF[100];
	
	static double m_lfLogTagLenCDF[120][21];
	
	CPSMConf m_psmconf;
	FILE * m_fp;
	
	FILE * m_debugFp;
	
	CXLinkPepResult * m_pPeptide;
	CSpectrum * m_pSpectrum;
	size_t m_tionmzSize;
	int m_iontag_pep1[10*MAX_IONTYPE_NUM][MAX_PEPTIDE_LENGTH];
	int m_iontag_pep2[10*MAX_IONTYPE_NUM][MAX_PEPTIDE_LENGTH];
	
	bool * m_pbMatched;
	vector<struct PeakInfo> m_vPeakInfo;
	vector<vector<int> > m_vMatchedIons;
	
	//CMzTriple * m_ionmz;
	vector<CMzTriple> m_ionmz;

	// for single peptide;
	CMzTriple * m_SingleIonmz;
	size_t m_tSingleIonMzSize;
	int * m_pSingleMatched;
	double * m_lfSingleMatched;
	/*
	bool m_single_iontag_pep1[MAX_IONTYPE_NUM][MAX_PEPTIDE_LENGTH];
	bool m_single_iontag_pep2[MAX_IONTYPE_NUM][MAX_PEPTIDE_LENGTH];
	*/
	
	//aamass precomputing
	int m_nAAMass[26];
	
	bool m_bComputeMz;
	//to compute the multiplier according to the m_cTolType
	double m_lfTolMultiplier;
	
	char m_cTolType;
	size_t m_tSingleIonNum1;
	size_t m_tSingleIonNum2;
	
	size_t * m_pPeakRank;
	
	vector<struct IonMatchInfoEx> m_vIonMatchInfo;
public:
	vector<struct IonMatchInfo> m_vIonTypeMatchInfo;
	
	vector<struct IonMatchInfo> m_vCommonIonTypeMatchInfo;
	vector<struct IonMatchInfo> m_vXlinkIonTypeMatchInfo;

	double m_lfTop10Ratio;
	double m_lf10PercentRatio;
	double m_lfTotalRatio;
	double m_lfNRErrorRatio;
};

#endif /*PSMATCH_H_*/
