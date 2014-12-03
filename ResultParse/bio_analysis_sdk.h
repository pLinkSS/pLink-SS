#ifndef BIO_ANALYSIS_SDK_H_
#define BIO_ANALYSIS_SDK_H_

#include <set>
#include <string>
#include <sstream>
#include <vector>
#include <string>
#include <list>
#include <vector>
#include <map>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>
#include <limits>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <functional>
#include <queue>

using namespace std;
//#include "PSM/sdk.h"
#include "option.h"



//using namespace proteomics_sdk;

//#define TESTOUTPUT
//#define TESTPARMFILE
//#define OUTPUTPAEMFILE
//#define DEBUG
//#define dataanalysisforhfchen

#ifdef WIN32
const char cSlash = '\\';
#else
const char cSlash = '/';
#endif

#ifdef WIN32
#define PRI_SIZE_T "%u"
#else
#define PRI_SIZE_T "%zu"
#endif

namespace bio_analysis
{

enum SearchEngineType
{
	ST_PFIND = 0, ST_MASCOT = 1, ST_SEQUEST = 2, ST_SQT = 3, ST_NEXT_ENGINE = 10
};

struct KeyWordSpec
{
	string m_strInPutPath;
	string m_strFirstPeptide;
	double m_lfScore;
	//vector<string> m_vProtein;
};

const double protonH = 1.00727647012;
const int SIMPLE_MAX = 32;
/* maximum length of file name */
const size_t FILE_NAME_SIZE = 256;

enum ExprotType
{
	ST_INTERFACEEXPORT = 0,
	ST_HTMLFormat = 1,
	ST_PBULIDPREFormat = 2,
	ST_PFINDFormat = 3,
	ST_PEPXMLFormat = 4,

	ST_NEXT_Format = 10
};

const int SCORE_SIZE[ST_NEXT_ENGINE] =
{ 2, 2, 2 };

class CBioMethods //relative error
{
public:
	static inline double TransDeltaMH2PPM(double lfDelta, double lfMH)
	{
		return lfDelta * 1000000 / lfMH;
	}
	static double CalcpI()
	{
		return 0.0;
	}
	static void ExtractScanNo()
	{
	}
};

class CFiltration
{
public:
	bool m_UseFilter;
	double m_lfPepMassLowerBound;
	double m_lfPepMassUpperBound;

	size_t m_tLengthLowerBound;
	size_t m_tLengthUpperBound;

	//	vector<double> m_vlfPepTolBase;
	//	string m_strPepTolBaseType;

	bool m_Sepatrate;
	vector<double> m_lfPepTolBase;
	vector<double> m_lfPepTolLowerBound;
	vector<double> m_lfPepTolUpperBound;
	string m_strPepTolType;

	vector<int> m_vChargeState;

	size_t m_tRankLimit;

	int m_bRedundant;

	bool m_bFixedDeltCn;
	double m_lfFixedDeltCn;

	bool m_bUseFDR;
	double m_lfFDR;
	int m_FDRFormula;
	double m_ScoreMax;
	double m_ScoreMin;

	size_t m_nDistinctPepLimit;
	size_t m_nDistinctSpecLimit;

	string m_strModSites;
	string m_strCTerminal;
	string m_strNTerminal;
};

struct CModificationSiteInfo
{
	string m_strModName;
	size_t m_tPos;
};

class CMatchProteinInfo //Protein information
{
public:
	string m_strAC;
	string m_strSQ;
	string m_strDE;
};

class CProteinInfo
{
public:
	double m_MW;
	double m_pI;
	double m_Coverage;
	string m_strDE;
};

class CMatchPeptideInfo
{

public:
	~CMatchPeptideInfo()
	{
		clear();
	}
	string m_strSQ;
	char m_cPrev;
	char m_cNext;//pre of next of the pep，R or K
	vector<double> m_vlfScores;
	//SearchEngineType m_eEngineType;
	double m_lfCalc_MH;
	//double m_lfCalc_M;
	vector<CModificationSiteInfo> m_vMod;
	vector<string> m_vProteinAC;
	//size_t m_tRank;
	int m_tRank;
	double m_lfDelta;
	double m_lfPPM;
	//bool m_bDecoyPep;//
	void clear()
	{
		m_vlfScores.clear();
		m_vMod.clear();
		m_vProteinAC.clear();
	}
};

class CMatchSpectraInfo
{
public:
	~CMatchSpectraInfo()
	{
		clear();
	}
	string m_strFileName;
	int m_nCharge;
	double m_lfMH;
	vector<CMatchPeptideInfo> m_vPeptides;
	//	int m_nStartScan;
	//	int m_nEndScan;
	int m_nSpectraNum;
	int m_nDataSetID;//SamplerID
	int m_nFileID;//FileID
	SearchEngineType EnginType;
	void clear()
	{
		for (size_t t = 0; t < m_vPeptides.size(); t++)
			m_vPeptides[t].clear();
		m_vPeptides.clear();
	}
};

class CConditionInfo //pFind txt //one file of the input
{
private:
	string m_strFilePath;//the path of this file
public:
	int m_nTotalSpec;
	string m_strInputPath;//
	vector<string> m_vDatabase;
	string m_Enzyme;
	vector<string> m_vFixMod;
	vector<string> m_vVarMod;
	bool m_bPepMono; /*default: true*/
	double m_lfPepBase; /*default: 0.0*/
	double m_lfPepTol; /*TOL=*/
	string m_strPepTolType; /*TOLU=*/
	bool m_bFragMono; /*MASS=*/
	double m_lfFragBase; /*default: 0.0*/
	double m_lfFragTol; /*ITOL=*/
	string m_strFragTolType; /*ITOLU=*/
	int m_nMaxMiss; /*?PFA=*/
	string m_strInstrument; /*INSTRUMENT=*/
	SearchEngineType EngineType;
};

class CParseredDataSetInfo //the data of one browse(ID)
{
public:
	int m_nDataSetID;
	SearchEngineType m_eEngineType;
	vector<CConditionInfo> m_vConditions; //the front of one file dataset
	//CConditionInfo m_vCondition;
	//first: various charge states; second: one data set
	map<int, vector<CMatchSpectraInfo> > m_mapResultData;//separate by charge
};

struct DataSetInputPathInfo //pFind, Mascot, SEQUEST result //one browse

{
	int m_nDataSetID;
	SearchEngineType m_Type;
	vector<string> m_vResultPath;
	vector<int> m_Valid;//避免使用vector<bool>
	vector<CConditionInfo> m_vConditions;//记录对应文件来自哪个谱图文件
};

class CConf //the message of the dlg
{
public:
	string m_SampleID;
	vector<DataSetInputPathInfo> m_vEngineInputInfo;
	vector<string> m_vProDBPath;
	vector<string> m_vDecoyTags;
	CFiltration m_Filter;
	string m_OutPutForder;//输出的文件夹
	string m_outPutForder_Java;
	string m_outPutForder_Charge;
	string m_OutPutForder_ErrorFigure;
	string m_outPutForder_Index;
	string m_outPutForder_pLabel;
	string m_OutPutFile;//文件文件名
	string m_OutPutName;//输出文件路径m_OutPutForderPath+m_OutPutForderPath
	vector<ExprotType> ExportFormat;
	string ExportFileType;
	string strExportPath;
	time_t time_now;

	int m_bReserveDecoy;

	int m_bCrossLink;
	vector<int> m_vLinkID;
	vector<int> m_vLinkType;

	int m_ReParse;//为了解决pFind调pBuild时必须重新解析

	string m_bDataAnalysis;

	string m_strPNovoParam;
};

typedef pair<string, long> SPECTRA_PATH;
typedef pair<SPECTRA_PATH, KeyWordSpec> SPECTRAINFO;
typedef pair<CConditionInfo, map<int, vector<SPECTRAINFO> > > OneDATASET;
//typedef map<int, vector<SPECTRAINFO> > OneDATASET;
typedef pair<string, vector<CModificationSiteInfo> > PEPTIDEINFO;
//typedef map<string, long> PROTEINIDATABASEINFO;


class IntegratedInfo
{
public:
	map<string, vector<SPECTRAINFO> > spec_peptide;
	map<string, vector<SPECTRAINFO> > peptide_spec;
	map<string, vector<SPECTRAINFO> > protein_peptide;
	map<string, set<string> > peptide_protein_set;

	set<string> m_setAllModName;
	set<string> m_setInvalidProteins;
	map<string, set<string> > m_mapSameSet;
	map<string, set<string> > m_mapSubSet;

	set<string> m_deleteProtein;
	set<string> m_deletePeptide;
	set<string> m_deleteSpectra;

	void clear()
	{
		spec_peptide.clear();
		peptide_spec.clear();
		protein_peptide.clear();
		peptide_protein_set.clear();

		m_setAllModName.clear();
		m_setInvalidProteins.clear();
		m_mapSameSet.clear();
		m_mapSubSet.clear();

		m_deleteProtein.clear();
		m_deletePeptide.clear();
		m_deleteSpectra.clear();

	}
};
///////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////

}
;

#endif /*BIO_ANALYSIS_SDK_H_*/
