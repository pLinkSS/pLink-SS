#ifndef INTERFACE_H_
#define INTERFACE_H_

#include "sdk.h"
//#include "predefine.h"
//#include <string>
using namespace proteomics_sdk;
using namespace std;
namespace proteomics_sdk
{

/*
 * preprocess algorithm class
 */
class CPreProcess
{
public:
	CPreProcess(void){};
	virtual ~CPreProcess(){};

	virtual void Init(CCondition & condition) = 0;
	virtual bool Instrument(InstrumentType eType) const = 0;
	virtual void Run(CSpectrum Input, CSpectrum &Output) = 0;
	virtual void Close(void) = 0;
};

/*
 * Score algorithm class
 */
class CPeptideScore
{
public:
	CPeptideScore(void){};
	virtual ~CPeptideScore(){};
	
	virtual int GetVersion(void) const = 0;

	virtual void Initialize(const CCondition & cond) = 0;
	virtual void Close(void) = 0;

	virtual void SetSpectrum(CSpectrum & spec) = 0;
	virtual void SetPeptide(const CSimplePeptide & pep) = 0;
	
	virtual double Score() = 0;
};

/*
 * Evaluation algorithm class
 */
class CPeptideEvaluater
{
public:
	CPeptideEvaluater(void){};
	virtual ~CPeptideEvaluater(){};
	virtual void Init(const CCondition &)=0;

	// modify by emily , remove 'const' for the convenience of refined-search
	virtual void Run(CSpectrum & spectrum, CSimpleMatchResult * pResult) = 0;

	virtual void Close(void)=0;
};

/*
 * Protein report class
 */
class CProteinReport
{
public:
	CProteinReport(void){};
	virtual ~CProteinReport(){};

	virtual void Init(const CCondition& condition,clock_t clkStartTime = 0) = 0;

	virtual void GetLines(const vector<CSimpleMatchResult> & vResults, const vector<CAssignedProtein> & vProteins, const vector<CSpectrum> & vSpectra, string strInputPath, string & strLines) = 0;

	virtual void WriteFile(const vector<CSimpleMatchResult> & vResults, const vector<CAssignedProtein> & vProteins, const vector<CSpectrum> & vSpectra, string strInputPath, string strOutputPath, string strTitle = "") = 0;

	virtual void Close(void) = 0;

	virtual ProteinReportType GetType(void) = 0;
};
/*
 * Protein inference class
 */
class CProteinInfer
{
public:
	CProteinInfer(void){};
	virtual ~CProteinInfer(){};

	virtual void SetSign(string strFPSign) = 0;
	virtual vector<CAssignedProtein> & Infer(vector<CSimpleMatchResult> & vResults, const vector<CSpectrum> & vSpectra, 
			FILTER_CRITERIA_INFO & cs) = 0;
	virtual ProteinInferType GetType(void) = 0;
	virtual void Initialize(CCondition * pCond) = 0;
};

class CMS1Input
{
public:
	CMS1Input(void) {};
	virtual ~CMS1Input() {};

	virtual MS1FormatType GetType(void) = 0;
	virtual void LoadAll(string strPath, vector<CMS1Info> & S) = 0;

	virtual void StartLoad(string strPath, int & nTotal) = 0;
	virtual void LoadNext(CMS1Info & spec, int &idx) = 0;
	virtual void LoadNextBatch(int nNumber ,vector<CMS1Info> & S, int &idx) = 0;
	virtual void EndLoad(void) = 0;
};

class CMS2Input
{
public:
	CMS2Input(void) {};
	virtual ~CMS2Input() {};

	virtual MS2FormatType GetType(void) = 0;
	
	virtual void LoadAll(string strPath, vector<CSpectrum> & S) = 0;

	virtual void StartLoad(string strPath, int & nTotal) = 0;
	virtual void LoadNext(CSpectrum & spec, int &idx) = 0;
	virtual void LoadPrev(CSpectrum & spec, int &idx) = 0;
	virtual void LoadNextBatch(int nNumber ,vector<CSpectrum> & S, int &idx) = 0;
	virtual void LoadSpec(CSpectrum & spec, string strTitle, int& nSpecNo) = 0;
/*	virtual void LoadAll(string strPath, vector<CMS2Info> & S) = 0;

	virtual void StartLoad(string strPath, int & nTotal) = 0;
	virtual void LoadNext(CMS2Info & spec, int &idx) = 0;
	virtual void LoadNextBatch(int nNumber ,vector<CMS2Info> & S, int &idx) = 0;*/
	virtual void EndLoad(void) = 0;
};

class CMS2Output
{
public:
	CMS2Output(void) {};
	virtual ~CMS2Output() {};

	virtual MS2FormatType GetType(void) = 0;
	
	virtual void WriteAll(string strPath, vector<CSpectrum> & S) = 0;

	virtual void StartWrite(string strPath, int & nTotal) = 0;
	virtual void WriteNext(CSpectrum & spec, int &idx) = 0;
	virtual void WriteNextBatch(int nNumber ,vector<CSpectrum> & S, int &idx) = 0;	
	virtual void EndWrite(void) = 0;
};

class CAAConf;
class CPeptideGenerator{
public:
    CPeptideGenerator(void) {;};
    virtual ~CPeptideGenerator() {;};    
    virtual void Init(const CAAConf & aaconf, bool bMono) = 0;
    virtual void GetPeptide(double lfMass, string & strOut) = 0;
    virtual void Batch(double lfMass, size_t tNum, vector<string> & vstrOut) = 0;
    virtual void Close(void) = 0;
};


// used by MPI ,added by wpwang 090417
//
/*
 * spectra index creator class
*/
class CSpectraIndexCreator
{
public:
	CSpectraIndexCreator(void) {;};
	virtual ~CSpectraIndexCreator() {;};
	virtual void Init(const CCondition & nCondition)=0;
	//strOutPath is a directory
	virtual void Create(string strTitle, string strSpectraPath, MS2FormatType eType, string strOutPath, int nMinBlockNum, int nMaxBlockCapacity,string & strMetaFile, int & nBlockNumber, int & nSpectraTotal)=0;
	virtual void Close(void)=0;
	
};

/*
 * spectra index loader class
*/
class CSpectraIndexLoader
{
public:
	CSpectraIndexLoader(void) {;};
	virtual ~CSpectraIndexLoader() {;};
	
	virtual void Init(string strMetaPath, int & nBlockNum, int & nTotalSpectra) = 0;
	virtual void Close(void) = 0;


	virtual void Load(int nBlockId, vector<CSpectrum> & vS) = 0;

	
	virtual size_t GetSpectraNum( void)=0;
	virtual void GetRange(double & lfMHMin,double & lfMHMax)=0;

};
//
// used by MPI ,added by wpwang 090417


class CAppender
{
public:
	CAppender(void){;};
	virtual ~CAppender(){;};
	
	virtual void Init(const CCondition &condition) = 0;
	virtual void Out(const string &strOunt) = 0;
	virtual void Close(void) = 0;
};

}
namespace proteomics_search
{

// modify 090324
class CFlow
{
public:
	CFlow(vector<CSpectrum> & vSpec){};

	virtual ~CFlow(){};

	virtual void Init(const CCondition &condition, CSearchState * pState, size_t tID) = 0;

	virtual void Run(const string & strTmpFile) = 0;
	
	virtual void Close(void) = 0;
	
	virtual CSearchState * GetSearchState(void) = 0;
	
	virtual string GetVersion(void) const = 0;
};


}
#endif /*INTERFACE_H_*/
