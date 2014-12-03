#ifndef SEARCHENGINE_H_
#define SEARCHENGINE_H_


#include "../ProteinIndex/ProteinIndex.h"
#include "../ProteinIndex/ProteinReaderFactory.h"
#include "../ProteinIndex/ProteinReader.h"
#include "../ProteinIndex/DiskRandomSQReader.h"
#include "../ProteinIndex/DiskRandomACReader.h"

using namespace ProteinIndex;
using namespace proteomics_sdk;
#include <pthread.h>

//class proteomics_sdk::CAppender;
class proteomics_sdk::CTrace;
namespace proteomics_search
{
class CSearchEngine
{
public:
	CSearchEngine();
	virtual ~CSearchEngine();
	
    // Interface
	bool Init(const string & strParamPath);
	void Search();
	void Infer();
	void GetSeachInfo(int & nTempFileNum , int & nTotalSpec , string & strIdentifier);
	void Close(void);

	// multi-thread control functions
	void _SetFreeSignal(size_t tID, bool bFree){m_vFlowEnd[tID] = bFree;}
	bool _IsFree(size_t tID){return m_vFlowEnd[tID];}
	inline void _SetLock(){m_bLocked = true;}
	inline void _FreeLock(){m_bLocked = false;}
	bool _IsLocked(){return m_bLocked;}
	void _SetCurrentID(size_t tID){m_tCurrentThreadID = tID;}
	inline size_t _GetCurrentID(){return m_tCurrentThreadID;}
	CFlow * _GetCurrentThread(){return m_vFlow[m_tCurrentThreadID];};
	void _WaitForLock()
	{
		while(_IsLocked())
		{
#ifdef WIN32
			Sleep(100);
#else
			usleep(100);
#endif
		}
	}
	void _WaitAllThreads();
	
protected:
    // work for search
    bool _ParseCmdLine(int nArgc, char ** pArgv);
	void _LoadSpectra();
    void _ReadOutputPath();
    void _InitThreadPool();
    static void * _SearchEngine_Callback(void* pArguments);

	// work for refined-search
	bool _GetFileName(string strSuffix,string & strFilePath);
	bool _LoadNextBatchOfEvcoef(string strEvcoefFile,int nBatchStart,int nLoad,vector<CSpectrum> & vSpecInfo);
    void _ReadSpectraType();
    bool _LoadSpectra(CMS2Input * pIO, 
			int nSpectraTotal, 
			int & nBatchStart, 
			int & nLoad,
			vector<CSpectrum> & vSpecInfo);
    int _PreLoadSpec(CMS2Input * pIO);
	
protected:
	// starting time 
	time_t m_tmStartTime;
	
    // parameter
	string m_strParamFile;      // the path of parameter file, for example, "D:\\a.pfind"
    CCondition m_Condition;
	string m_strWorkDir;        // the path of Searcher.exe

    // spectra input
	string m_strSpectraPath;    // input file path
    MS2FormatType m_eType;      // spectra file type
    int m_nLoad;
	int m_nTotalSpec;

    // output
	string m_strOutputPath;     // output path
	string m_strIdentifier;     // a 5-length string to indicate the current search
	
    // multi-threading
	vector<CFlow*> m_vFlow;     // search threads information
	vector<vector<CSpectrum> > m_vSpecInfos;        // the outside vector's size equals to the number of specified threads, and the inside vector store spectra corresponding to each thread
    pthread_mutex_t mutex;      // variables used in multi-threading
	pthread_cond_t cond;
	vector<bool> m_vFlowEnd;
	size_t m_tCurrentThreadID;
    bool m_bSearched;
    bool m_bLocked;
	int m_nTempFileNum;     // number of temporary files (*.pfd) already generated

    // Trace and State
	CSearchState * m_pState;
	CTrace * m_pTrace;
};
}
#endif /*Searcher_H_*/
