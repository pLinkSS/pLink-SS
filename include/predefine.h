#ifndef PREDEFINE_H_
#define PREDEFINE_H_
#include <cstddef>

//the maximum and minimum limitation of the peptide length
const int MAX_PEPTIDE_LENGTH = 100;
const int MIN_PEPTIDE_LENGTH = 3;

//mimimum precursor mass tolerance
const double MIN_PEPTIDE_TOL = 2.0;

//the maximum limitation of modification sites
const int MAX_MODIFY_NUM = 20;

//modify by emily
//the number of score that a spectrum should save at least
//MIN_SCORE_NUM = 0 means not use pep gen as default
const size_t MIN_SCORE_NUM = 0;

//the number of score that a spectrum can save at most
const size_t MAX_SCORE_NUM = 20000;


//kinds of elements now: CHNO
const int ELEMENT_NUMBER = 4;
//multiplication coefficient
const int MZMULTIPLIER = 10000;
const int MZMULTIPLIER2 = 1000; //For ion index. Add at 2014.2.21
const int MAXMRANGE = 10000; // max m considered in ion index  //Notice: need be changed when index changed

const int nDECIMAL = 100000;
const int MMU_MULTIPLIER = 10;

//maximal limitation of iontype
const int MAX_IONTYPE_NUM = 200;

//used in SpectraIO when output the float number
const int PRECISION = 6;

//maxmal limatation of proteins
const int MAXPTN = 10005;

//to set the precision of the comparison between two float numbers
const double EPS = 0.001;

//number of spectra identified during one scanning of the database
const int SALVO_BATCH_SIZE = 5000;

//number of peaks once reserved for a spectrum
const int MAX_PEAKS_NUM = 1024;

//the separator of the operating system path,respectively Windows and linux
#ifdef WIN32
const char SLASH = '\\';
#else
const char SLASH = '/';
#endif


#ifndef PATH_MAX
#define PATH_MAX 260
#endif

#ifndef MAX_PATH
#define MAX_PATH 260
#endif

//types of instruments
enum InstrumentType
{
	INSTRUMENT_UNKNOWN = 0,
	INSTRUMENT_ESI_QUAD_TOF = 1,
	INSTRUMENT_ESI_TRAP = 2,
	INSTRUMENT_ESI_QUAD = 3,     
	INSTRUMENT_ESI_FTICR = 4,    
	INSTRUMENT_ESI_4SECTOR = 5,
	INSTRUMENT_MALDI_QUAD_TOF = 6,
	INSTRUMENT_MALDI_TOF_PSD = 7,
	INSTRUMENT_MALDI_TOF_TOF = 8,
	INSTRUMENT_MALDI_QIT_TOF = 9,
	INSTRUMENT_FTMS_ECD = 10
};

//types of the indices
enum IndexContent
{
	// prefix PEPTIDE_PAIR represent flows for xlink
	PEPTIDE_PAIR = 2,
	PEPTIDE_PAIR_OPEN = 7,
};

//types of loading indices
enum IndexType
{
	INDEX_DEFAULT = 0,
	DISK_FILE = 1,
	MEMORY_FILE = 2,
	MEMORY_MAPPED = 3
};

// used by MPI ,added by wpwang 090417
//
enum SpectraIndexType 
{
	SPECTRA_INDEX_BIN = 0,
	SPECTRA_INDEX_BIN_BALANCE = 1,
	SPECTRA_INDEX_PREPROC = 2,
	SPECTRA_INDEX_PREPROC2 = 3,
	SPECTRA_INDEX_PREPROC3 = 4,
	SPECTRA_INDEX_SIMPLE = 5,
	SPECTRA_INDEX_LARGE = 6
};
enum LoadBalanceType
{//MassRange.cpp  MassShuffle.cpp Estimate1.cpp Estimate2.cpp DynamicM3.cpp
	LOAD_BALANCE_MASS_RANGE = 0,
	LOAD_BALANCE_MASS_SHUFFLE = 1,
	LOAD_BALANCE_MASS_DYNAMIC_MUL2 = 2,
	LOAD_BALANCE_MASS_DYNAMIC_MUL3 = 3,
	LOAD_BALANCE_MASS_DYNAMIC_MUL4 = 4,
	LOAD_BALANCE_ESTIMATE1 = 5,
	LOAD_BALANCE_ESTIMATE2 = 6,
	LOAD_BALANCE_ESTIMATE3 = 7
};
//
// used by MPI ,added by wpwang 090417

//types of preprocessing algorithms
enum PreProcType
{
	PRE_PROC_XLINK_HCD = 7
};

//types of scoring algorithms
typedef enum ScoreType
{
	SCORE_KSDP = 1
}SCORE_TYPE;

//types of searching flows
enum SearchFlowType
{
	FLOW_NULL = 0,
	FLOW_PRO_IDX = 1,
	FLOW_ALL_IDX = 2,
	FLOW_XLINK = 3,
};



//types of evaluate algorithms
enum EvaluateType
{
	EV_DEFAULT = 0
};

/*
 * Spectra format conversion class
 */
enum MS1FormatType
{
	PMF_RAW = 0,
	PMF_MZML = 1,
	PMF_MS1 = 2
};

enum MS2FormatType
{
	PFF_DTA = 0,
	PFF_DTAS = 1,
	PFF_MGF = 2,
	PFF_PKL = 3,
	PFF_MZML = 4,
	PFF_MS2 = 5,
	PFF_RAW = 6,
	PFF_SDTA = 7
};

//types of the protein files
enum DataFormatType
{
        FASTA = 0,
        DAT   = 1,
        XML   = 2,
        INDEX = 3
};

//unknown, perhaps added by Li You or Liu Luoqi
enum AssessType
{
	ASSESS_None = 0,
	ASSESS_Need = 1,
	ASSESS_NoNeed = 2,
	ASSESS_Discard = 3
};

//types of reporting the identification results of proteins
enum ProteinReportType
{
	PROTEIN_REPORT_TXT  = 0,
	PROTEIN_REPORT_CANDIDATE = 1
};

//types of proteins inferring
enum ProteinInferType
{
	PROTEIN_INFER_DEFAULT = 0,
	PROTEIN_INFER_AC = 1,
	PROTEIN_INFER_NULL = 2,
	PROTEIN_INFER_SIMPLE = 3
};

//types of sorting method of the spectra
enum SpectrumSorterType
{
	SPECTRUM_MASS_GREATER = 0,
	SPECTRUM_MASS_LESS = 1
};

//types of generating peptides
enum PepGeneType
{
	PEPGENE_DEFAULT = 0
};

enum ModType
{
	MT_NORMAL = 0,
	MT_PEP_NTERM = 1,
	MT_PEP_CTERM = 2,
	MT_PRO_NTERM = 3,
	MT_PRO_CTERM = 4
};

enum LogRankType
{
	LOG_RANK_NONE=0,
	LOG_RANK_ALERT=100,
	LOG_RANK_INFO=200,
	LOG_RANK_DEBUG=300
};

enum LogAppenderType
{
	LOG_APPENDER_CMD=0,
	LOG_APPENDER_FILE=1,
	LOG_APPENDER_CMD_AND_FILE=2
	
};

enum ModuleType
{
	MODULE_ALL=0,
	MODULE_SEARCHER = 1,
	MODULE_SEARCHENGINE = 2,
	MODULE_CONDITION = 3,
	MODULE_SPECTRAIO=4,
	MODULE_FLOW = 5,
	MODULE_PREPRO = 6,
	MODULE_SCORE = 7,
	MODULE_MASS2PEP = 8,
	MODULE_PROIDX = 9,
	MODULE_GENERATOR = 10,
	MODULE_EVALUATE = 11,
	MODULE_PROINFER = 12,
	MODULE_INFEROR = 13,
	MODULE_IMPORTER = 14,
	MODULE_SPECINDEX = 15,
	MODULE_LOADBALANCE = 16,
	MODULE_MPISEARCHER = 17,
	MODULE_PFDIO = 18,
	MODULE_ITERFILTER = 19
};

// pfind-xlink
// types of xlink
enum XLinkType
{
	NONE_LINK = 0,	// no linker exist
	MONO_LINK = 1,	// hydrolyzed linker on one peptide similar to a modification
	LOOP_LINK = 2,	// loop linker on one peptide
	X_LINK = 3		// cross linker on two peptides
};

//criteriation of filterring the protein identification results
typedef struct
{
	double lfFPR;
	double lfPepMass;
	size_t nPepLength;

	int nProUniquePepNum;
	int nProSpecNum;

}FILTER_CRITERIA_INFO;

//index structure for sorting the spectra
struct SPEC_SORTER_INDEX_INFO
{
	int nIndex;
	double lfMass;
	size_t nCandidateNum;
	int nCharge;
	double lfMin;
	double lfMax;
};

//recording the ion type information
struct CIonType
{
	bool bNTerm;
	char cType;
	// 0 : N-term ions
	// 1 : C-term ions
	// 2 : Internal ions
	// 3 : Precursor ions (one peptide) 
	// 4 : Precursor ions (two peptide)
	// 5 : Reporter ions
	int nCharge;
	int nLost[ELEMENT_NUMBER];
	int nTotalLostVal;
	bool bIntraContinous;
	int nInterContinousWnd1;
	int nInterContinousWnd2;
	int nContainLinker;
	//0: both 
	//1: common ion: not contain linker
	//2: xlink ion: contain linker
	char cSymbol;//added at 2013.11.21, record b\y\a...
};



#endif /*PREDEFINE_H_*/
