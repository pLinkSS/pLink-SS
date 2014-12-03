#include <fstream>
#include <time.h>
#include <dirent.h>
#include <stdio.h>
#include <string>
#include <exception>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <limits>
#include <sys/stat.h>

#ifdef WIN32
#include <direct.h>
#endif

#include "../include/predefine.h"
#include "../include/sdk.h"
#include "../include/option.h"
#include "../include/interface.h"
#include "../SpectraIO2.0/MS2Input/MS2InputFactory.h"
#include "../Flow/FlowFactory.h"
#include "../ProteinIndex/ProteinIndex.h"
#include "../ProteinIndex/ProteinReaderFactory.h"
#include "../ProteinIndex/ProteinReader.h"
#include "../ProteinIndex/DiskRandomSQReader.h"
#include "../ProteinIndex/DiskRandomACReader.h"
#include "SearchEngine.h"


#ifdef XLINK_TEST
#include "../Flow/XLinkPepFlow.h"
#endif

using namespace std;
using namespace proteomics_sdk;
using namespace proteomics_search;
using namespace ProteinIndex;

bool bAllBusy = false;

#ifdef WIN32
const char cSlash = '\\';
#else
const char cSlash = '/';
#endif

CSearchEngine::CSearchEngine() :
		m_pState(NULL), m_nTempFileNum(0), m_eType(PFF_DTA), m_nLoad(0), m_nTotalSpec(
				0), m_bLocked(0), m_tCurrentThreadID(0) {
	m_bSearched = false;
	m_tmStartTime = time(NULL);

	char szBuf[BUFFER_SIZE] = { 0 };
	getcwd(szBuf, BUFFER_SIZE);
	m_strWorkDir = szBuf;
	if (m_strWorkDir[m_strWorkDir.length() - 1] != cSlash) {
		m_strWorkDir += cSlash;
	}

	m_pTrace = CTrace::GetInstance();
}

CSearchEngine::~CSearchEngine() {
	Close();

	CTrace::DeleteInstance();
	m_pTrace = NULL;
}

bool CSearchEngine::_GetFileName(string strSuffix, string & strFilePath) {
	if (strSuffix.size() <= 0)
		return false;

	string strTempDir = m_Condition.m_strRSPath;

	if (strSuffix[0] != '.')
		strSuffix += ".";

	string strFileName = m_Condition.m_strRSTag + strSuffix;
	strFilePath = m_Condition.m_strRSPath + strFileName;

	DIR* dir = opendir(strTempDir.c_str());
	if (NULL != dir) {
		struct dirent* de;
		while ((de = readdir(dir)) != NULL) {
			if (strcmp(de->d_name, ".") != 0 && strcmp(de->d_name, "..") != 0) {
				if (strcmp(de->d_name, strFileName.c_str()) == 0) {
					return true;
				}
			}
		}
		closedir(dir);
	}

	return false;

}

void CSearchEngine::Infer() {
#ifdef WIN32
	string strCommand = "Inferor.exe ";
#else
	string strCommand = "./Inferor ";
#endif

	CCheck::CheckPath(m_strParamFile);
	strCommand += m_strParamFile;

	strCommand += " ";
	if (m_strOutputPath.find(' ') != string::npos) {
		strCommand += "\"" + m_strOutputPath;
#ifdef WIN32
		strCommand += "\\";
#endif
		strCommand += "\"";
	} else {
		strCommand += m_strOutputPath;
	}
	strCommand += " ";
	char szBuf[256] = { 0 };
	sprintf(szBuf, "%d %d %s %ld", m_nTempFileNum, m_nTotalSpec,
			m_strIdentifier.c_str(), m_tmStartTime);
	strCommand += szBuf;

	if (0 != system(strCommand.c_str())) {
		CErrInfo info("CSalvoFlow5", "Run", "in the function system");
		info.Append("strCommand=" + strCommand);
		throw invalid_argument(info.Get().c_str());
	}

	/* remove all temporary pfd files */
	CTrace::GetInstance()->Info("Removing temporary pfd files...");
	strCommand = "del ";
	strCommand += m_strOutputPath;
	if (strCommand[strCommand.length() - 1] != '\\')
		strCommand += "\\";
	strCommand += "*.pfd";
	if (0 != system(strCommand.c_str())) {
		CErrInfo info("SearchEngine", "Infer",
				"Cannot remove temporary pfd files: " + strCommand);
		CTrace::GetInstance()->Info(info.Get().c_str());
	}
	CTrace::GetInstance()->Info("The search is completed.", MODULE_MPISEARCHER);
	m_pTrace->Info("pFind Searching completed successfully!");
}

bool CSearchEngine::Init(const string & strParamFile) {

	m_bSearched = false;
	Close();

	m_tmStartTime = time(NULL);

	m_pState = new CSearchState();

	m_pTrace->Info("SearchEngine initialize...");
	m_pTrace->Info("Initializing pFind...");

	m_strParamFile = strParamFile;
	m_pTrace->Info("Reading Parameter File : " + strParamFile);

	CConditionReader reader(m_strParamFile, m_strWorkDir);
	try {
		reader.Read();
	} catch (exception & e) {
		CErrInfo info("CSearchEngine", "Init",
				"in the function CConditionReader::Read");
		info.Append("m_strParamFile=" + m_strParamFile);
		info.Append("m_strWorkDir=" + m_strWorkDir);
		throw runtime_error(info.Get(e).c_str());
	} catch (...) {
		CErrInfo info("CSearchEngine", "Init",
				"caught an unknown exception in the function CConditionReader::Read");
		info.Append("m_strParamFile=" + m_strParamFile);
		info.Append("m_strWorkDir=" + m_strWorkDir);
		throw runtime_error(info.Get().c_str());
	}

	if (!reader.m_Condition.ValidateAll())
		return false;
	m_pTrace->Debug("Return from Condition Reader");
	m_Condition = reader.m_Condition;

	m_pTrace = CTrace::GetInstance(m_Condition); // take the same effect as m_pTrace.Set(m_Condition)
	m_pTrace->Info("Read Parameter File Successful");

	m_pTrace->Info("Create Thread Pool...");
	try {
		_InitThreadPool();
	} catch (exception & e) {
		CErrInfo info("CSearchEngine", "Init",
				"in the function InitThreadPool");
		throw runtime_error(info.Get(e).c_str());
	} catch (...) {
		CErrInfo info("CSearchEngine", "Init",
				"caught an unknown exception in the function InitThreadPool.");
		throw runtime_error(info.Get().c_str());
	}

	m_pTrace->Info("Initialization Finished!");
	return true;
}

void CSearchEngine::GetSeachInfo(int & nTempFileNum, int & nTotalSpec,
		string & strIdentifier) {
	if (!m_bSearched) {
		CErrInfo info("CSearchEngine", "GetSeachInfo", "bSearched == false");
		throw runtime_error(info.Get().c_str());
	}

	nTempFileNum = m_nTempFileNum;
	nTotalSpec = m_nTotalSpec;
	strIdentifier = m_strIdentifier;

}
void CSearchEngine::Search() {
	if (m_vFlow.empty() || !m_Condition.ValidateAll()) {
		CErrInfo info("CSearchEngine", "Search",
				"m_vFlow.empty() || !m_Condition.ValidateAll()");
		throw runtime_error(info.Get().c_str());
	}

	m_pTrace->Info("Running...");

	_ReadSpectraType();
	_ReadOutputPath();

	try {

		_LoadSpectra();

		m_bSearched = true;

	} catch (exception & e) {
		CErrInfo info("CSearchEngine", "Search",
				"in the function _LoadSpectra.");
		char temp[200] = { 0 };
		sprintf(temp, "FormatType=%d", m_eType);
		info.Append(temp);
		info.Append(m_strSpectraPath.c_str());
		info.Append(m_strOutputPath.c_str());
		throw runtime_error(info.Get(e).c_str());
	} catch (...) {
		CErrInfo info("CSearchEngine", "Search",
				"caught an unknown exception in the function _LoadSpectra.");
		char temp[200] = { 0 };
		sprintf(temp, "FormatType=%d", m_eType);
		info.Append(temp);
		info.Append(m_strSpectraPath.c_str());
		info.Append(m_strOutputPath.c_str());
		throw runtime_error(info.Get().c_str());
	}

	char szTemp[8192] = { 0 };
	sprintf(szTemp, "Total time consumed: %.3lfs",
			(clock() - m_pState->m_lfStartTime) / CLOCKS_PER_SEC);
	m_pTrace->Info(szTemp);
}

void CSearchEngine::Close(void) {
	for (size_t i = 0; i < m_vFlow.size(); ++i) {
		if (m_vFlow[i]) {
			delete m_vFlow[i];
			m_vFlow[i] = NULL;
		}
	}

	if (m_pState) {
		delete m_pState;
		m_pState = NULL;
	}

	m_bSearched = false;
}

void CSearchEngine::_InitThreadPool() {
//	m_pTrace->Info("Start _InitThreadPool()");

	m_vSpecInfos.clear();
	m_vFlow.reserve(m_Condition.m_tProcessorNum);

	m_pState->m_vThreadState.clear();

	CFlowFactory factory;
	for (size_t i = 0; i < m_Condition.m_tProcessorNum; ++i) {
		m_vSpecInfos.push_back(vector<CSpectrum>());
	}
	m_vFlow.clear();

//	if (m_Condition.m_eIndexContent == PEPTIDE_TRI_ION)
//	{
//
//	} // If want to generate ion index out. May add to ccondition to transfer to pFlow.

	for (size_t i = 0; i < m_Condition.m_tProcessorNum; ++i) {
		CFlow * pFlow = factory.GetFlow(m_Condition.m_eIndexContent,
				m_vSpecInfos[i]);

		m_vFlow.push_back(pFlow);

		try {
			pFlow->Init(m_Condition, m_pState, i);
		} catch (exception & e) {
			if (pFlow)
				delete pFlow;
			CErrInfo info("CSearchEngine", "_InitThreadPool",
					"in the function CFlow::Init");
			throw runtime_error(info.Get(e).c_str());
		} catch (...) {
			if (pFlow)
				delete pFlow;
			CErrInfo info("CSearchEngine", "_InitThreadPool",
					"caught an unknown exception in the function CFlow::Init.");
			throw runtime_error(info.Get().c_str());
		}
		m_vFlowEnd.push_back(true);

		CThreadState ts;
		ts.m_tID = i;
		ts.m_nCurrentSpectra = 0;
		ts.m_nTotalSpectra = 0;
		m_pState->m_vThreadState.push_back(ts);
	}

//	m_pTrace->Debug("End _InitThreadPool()");
}

void CSearchEngine::_ReadSpectraType() {
	COptionTool * pOption = new COptionTool("spectrum", m_strParamFile.c_str());
	m_strSpectraPath = pOption->GetString("spec_path", "null");

	if ("null" == m_strSpectraPath) {
		delete pOption;
		CErrInfo info("CSearchEngine", "_ReadSpectraType",
				"null == m_strSpectraPath");
		throw runtime_error(info.Get().c_str());
	}

	if (access(m_strSpectraPath.c_str(), 0) == -1) {
		CTrace * pTrace = CTrace::GetInstance();
		pTrace->Alert(
				"the path of the MS/MS data: " + m_strSpectraPath
						+ " is not existed");

		delete pOption;
		CErrInfo info("CSearchEngine", "_ReadSpectraType",
				m_strSpectraPath + "is not existed");
		throw runtime_error(info.Get().c_str());
	}

	string strType = pOption->GetString("spec_type", "PKL");
	m_eType = PFF_PKL;
	if ("DTA" == strType) {
		m_eType = PFF_DTA;
		if (m_strSpectraPath[m_strSpectraPath.length() - 1] != cSlash)
			m_strSpectraPath += cSlash;
	} else if ("DTAS" == strType)
		m_eType = PFF_DTAS;
	else if ("MGF" == strType)
		m_eType = PFF_MGF;
	else if ("RAW" == strType)
		m_eType = PFF_RAW;
	else if ("MS2" == strType)
		m_eType = PFF_MS2;

	CTrace * pTrace = CTrace::GetInstance();
	pTrace->Info("the path of the MS/MS data: " + m_strSpectraPath);
	pTrace->Info("MS/MS data file format: " + strType);

	delete pOption;
	pOption = NULL;
}
void CSearchEngine::_ReadOutputPath() {
	COptionTool * pOption = new COptionTool("flow", m_strParamFile.c_str());
	m_strOutputPath = pOption->GetString("output_path", "null");
	if ("null" == m_strOutputPath) {
		delete pOption;
		CErrInfo info("CSearchEngine", "_ReadOutputPath",
				"null == m_strOutputPath");
		throw runtime_error(info.Get().c_str());
	}
	delete pOption;
	pOption = NULL;

	if (m_strOutputPath[m_strOutputPath.length() - 1] != cSlash) {
		m_strOutputPath += cSlash;
	}

	//20100309
	CTrace * pTrace = CTrace::GetInstance();
	pTrace->Info("output path: " + m_strOutputPath);

	if (access(m_strOutputPath.c_str(), 0) == -1) {
		pTrace->Alert(m_strOutputPath + "is not existed");
#ifdef WIN32
		if (-1 == mkdir(m_strOutputPath.c_str()))
#else
				if(-1 == mkdir(m_strOutputPath.c_str(), 777))
#endif
				{
			delete pOption;
			CErrInfo info("CSearchEngine", "_ReadOutputPath",
					"Cannot create the folder: " + m_strOutputPath);
			throw runtime_error(info.Get().c_str());
		} else {
			pTrace->Info(m_strOutputPath + "is created.");
		}
	}

	pTrace->Break();
}
void* CSearchEngine::_SearchEngine_Callback(void* pArguments) {

	CSearchEngine * pSearchEngine = (CSearchEngine *) pArguments;
	string strOutputPath = pSearchEngine->m_strOutputPath;
	char szBuf[10] = { 0 };
	sprintf(szBuf, "%d", pSearchEngine->m_nTempFileNum);
	strOutputPath += szBuf;
	size_t tID = pSearchEngine->_GetCurrentID();
	CFlow * pFlow = NULL;
	pFlow = pSearchEngine->_GetCurrentThread();
	int nLoad = pSearchEngine->m_nLoad;

	pSearchEngine->m_pState->SetThreadTotal(tID, nLoad);

	while (pSearchEngine->_IsLocked()) {
		pSearchEngine->_FreeLock();
	}
	strOutputPath += "." + pSearchEngine->m_strIdentifier + ".pfd";

	try {
		pFlow->Run(strOutputPath);
	} catch (exception & e) {
		if (pFlow)
			delete pFlow;
		CErrInfo info("SearchEngine.cpp", "_SearchEngine_Callback",
				"in the function CFlow::Run");
		throw runtime_error(info.Get(e).c_str());
	} catch (...) {
		if (pFlow)
			delete pFlow;
		CErrInfo info("SearchEngine.cpp", "_SearchEngine_Callback",
				"caught an unknown exception in the function CFlow::Run.");
		throw runtime_error(info.Get().c_str());
	}

	pSearchEngine->_SetFreeSignal(tID, true);
	pthread_mutex_lock(&pSearchEngine->mutex);
	pthread_cond_signal(&pSearchEngine->cond);
	bAllBusy = false;

	pSearchEngine->m_pState->SetCurrentSpectra(
			pSearchEngine->m_pState->m_nCurrentSpectra += nLoad);
	pSearchEngine->m_pState->SetThreadTotal(tID, 0);
	pSearchEngine->m_pState->SetThreadCurrent(tID, 0);
	pthread_mutex_unlock(&pSearchEngine->mutex);

	return 0;
}

// true : not complete
// false : complete

bool CSearchEngine::_LoadNextBatchOfEvcoef(string strEvcoefFile,
		int nBatchStart, int nLoad, vector<CSpectrum> & vSpecInfo) {

	FILE * fp;
	if (NULL == (fp = fopen(strEvcoefFile.c_str(), "rb"))) {
		CErrInfo info("CSearchEngine", "_LoadNextBatchOfEvcoef",
				"in the function fopen");
		throw invalid_argument(info.Get().c_str());
	}

	if (vSpecInfo.size() <= 0)
		return true;

	size_t tSizePerUnit = vSpecInfo[0].m_stEVCoef.GetSize();
	size_t tTotalUnit;
	if (0 != fseek(fp, 0, SEEK_END)) {
		fclose(fp);
		CErrInfo info("CSearchEngine", "_LoadNextBatchOfEvcoef",
				"in the function fseek");
		throw invalid_argument(info.Get().c_str());
	}
	tTotalUnit = ftell(fp) / tSizePerUnit;

	if (nBatchStart >= (int) tTotalUnit) {
		fclose(fp);
		return false;
	}

	if (0 != fseek(fp, tSizePerUnit * nBatchStart, SEEK_SET)) {
		fclose(fp);
		CErrInfo info("CSearchEngine", "_LoadNextBatchOfEvcoef",
				"in the function fseek");
		throw invalid_argument(info.Get().c_str());
	}

	int nLeft = tTotalUnit - nBatchStart;

	if (nLeft < nLoad)
		nLoad = nLeft;

	cout << " nBatchStart = " << nBatchStart << endl << " nLoad = " << nLoad
			<< endl << " tTotalUnit = " << tTotalUnit << endl << " nLeft = "
			<< nLeft << endl;

	for (int i = 0; i < nLoad; ++i) {
		if (false == vSpecInfo[i].m_stEVCoef.ReadFromFile(fp)) {
			fclose(fp);
			CErrInfo info("CSearchEngine", "_LoadNextBatchOfEvcoef",
					"in the function fread");
			throw invalid_argument(info.Get().c_str());
		}

		//cout << "spectra " << i << "	: [" << vSpecInfo[i].m_stEVCoef.lfCoef0 << "," << vSpecInfo[i].m_stEVCoef.lfCoef1 << "]" << endl;
	}

	fclose(fp);

	if (nBatchStart + nLoad >= (int) tTotalUnit)
		return false;
	else
		return true;

}

bool CSearchEngine::_LoadSpectra(CMS2Input * pIO, int nSpectraTotal,
		int & nBatchStart, int & nLoad, vector<CSpectrum> & vSpecInfo) {
	nLoad = (int) m_Condition.m_tSalvoBatchSize;

	/*
	 * only a few spectra, for example, 2000, and 2 threads,
	 * batchsize=5000,
	 * then every thread is assigned with 1000 spectra rather than 
	 * assign all spectra to one thread
	 */
	if (nSpectraTotal
			< (int) (m_Condition.m_tSalvoBatchSize * m_Condition.m_tProcessorNum)) {
		nLoad = nSpectraTotal / (int) m_Condition.m_tProcessorNum;
		if (nLoad * (int) m_Condition.m_tProcessorNum < nSpectraTotal) {
			++nLoad;
		}
	}
	if (nSpectraTotal - nBatchStart < nLoad)
		nLoad = nSpectraTotal - nBatchStart;
	int nIdx = nBatchStart;

	try {
		int nTmpLoad = nLoad;
		int nTmpIdx = nIdx;
		pIO->LoadNextBatch(nLoad, vSpecInfo, nIdx);

		/*
		 * refined-search
		 * if loadinfo = 1 then load coef from *.coef
		 * seems no use
		 */
		/* if (m_Condition.m_bRSLoadInfo) {
			string strEvcoefFile;
			if (_GetFileName(".coef", strEvcoefFile)) {
				_LoadNextBatchOfEvcoef(strEvcoefFile, nTmpIdx, nTmpLoad,
						vSpecInfo);
			}
		} */
	} catch (runtime_error & e) {
		delete pIO;
		CErrInfo info("CSearchEngine", "_LoadSpectra",
				"CMS2Input::LoadNextBatch failed.");
		throw runtime_error(info.Get(e).c_str());
	} catch (...) {
		delete pIO;
		CErrInfo info("CSearchEngine", "_LoadSpectra",
				"caught an unknonw exception from CMS2Input::LoadNextBatch.");
		throw runtime_error(info.Get().c_str());
	}

	nBatchStart += nLoad;
	return (nBatchStart < nSpectraTotal);
}

void CSearchEngine::_LoadSpectra() {
	/* generate a 5-length string to indicate the current search */
	m_strIdentifier = "";
	srand(time(0));
	for (size_t i = 0; i < 5; ++i) {
		m_strIdentifier += (char) ('A' + rand() % 26);
	}

	CMS2InputFactory IOFactory;
	CMS2Input * pIO = IOFactory.GetImporter(m_eType);

	m_nTotalSpec = _PreLoadSpec(pIO);
	if (m_nTotalSpec <= 0) {
		delete pIO;
		return;
	}

	m_pTrace->Info("Count spectra completed.", MODULE_SEARCHENGINE);
	m_pState->SetTotalSpectra(m_nTotalSpec);

	int nBatchStart = 0;
	m_nLoad = 0;
	m_nTempFileNum = 0;
	bool bLoadingSpectra = true;
	pthread_mutex_init(&mutex, 0);
	pthread_cond_init(&cond, 0);

	while (bLoadingSpectra) {
		size_t nCnt = 0;
		for (size_t i = 0; i < m_Condition.m_tProcessorNum && bLoadingSpectra;
				++i) {
			_SetCurrentID(i);
			if (!_IsFree(i)) {
				++nCnt;
				continue;
			}
			_SetFreeSignal(i, false);
			try {

				bLoadingSpectra = _LoadSpectra(pIO, m_nTotalSpec, nBatchStart,
						m_nLoad, m_vSpecInfos[i]);

			} catch (exception & e) {
				CErrInfo info("CSearchEngine", "_LoadSpectra",
						"in function CSearchEngine::_LoadSpectra.");
				info.Append("m_strSpectraPath=" + m_strSpectraPath);
				info.Append("m_strOutputPath=" + m_strOutputPath);
				throw runtime_error(info.Get(e).c_str());
			} catch (...) {
				CErrInfo info("CSearchEngine", "_LoadSpectra",
						"caught an unknown exception.");
				info.Append("m_strSpectraPath=" + m_strSpectraPath);
				info.Append("m_strOutputPath=" + m_strOutputPath);
				throw runtime_error(info.Get().c_str());
			}

			pthread_t tid;
			_SetLock();
			pthread_create(&tid, 0, &CSearchEngine::_SearchEngine_Callback,
					this);
			_WaitForLock();
			++m_nTempFileNum;
		}

		/*
		 * when all threads are working, hold on the main thread
		 */
		if (bLoadingSpectra && (nCnt == m_Condition.m_tProcessorNum)) {
			pthread_mutex_lock(&mutex);
			bAllBusy = true;
			while (bAllBusy) {
				pthread_cond_wait(&cond, &mutex);
			}
			pthread_mutex_unlock(&mutex);

		}
	}

	_WaitAllThreads();

	pthread_mutex_destroy(&mutex);
	pthread_cond_destroy(&cond);

	pIO->EndLoad();

	delete pIO;
}

int CSearchEngine::_PreLoadSpec(CMS2Input * pIO) {
	if (!pIO)
		return 0;

	int nSpectraTotal = 0;
	m_pTrace->Info("Start loading Spectra...");
	try {
		pIO->StartLoad(m_strSpectraPath, nSpectraTotal);
	} catch (exception & e) {
		CErrInfo info("CSearchEngine", "_PreLoadSpec",
				"caught this runtime_error in CMS2Input::StartLoad Function.");
		info.Append("m_strSpectraPath=" + m_strSpectraPath);
		char temp[200] = { 0 };
		sprintf(temp, "nSpectraTotal=%d", nSpectraTotal);
		info.Append(temp);
		throw runtime_error(info.Get(e).c_str());
	} catch (...) {
		CErrInfo info("CSearchEngine", "_PreLoadSpec",
				"caught an unknown exception.");
		info.Append("m_strSpectraPath=" + m_strSpectraPath);
		char temp[200] = { 0 };
		sprintf(temp, "nSpectraTotal=%d", nSpectraTotal);
		info.Append(temp);
		throw runtime_error(info.Get().c_str());
	}
	return nSpectraTotal;
}

void CSearchEngine::_WaitAllThreads() {
	bool bContinue = true;
	while (bContinue) {
		bContinue = false;
		for (size_t i = 0; i < m_Condition.m_tProcessorNum && !bContinue; ++i) {
			if (!_IsFree(i)) {
				pthread_mutex_lock(&mutex);
				bAllBusy = true;
				while (bAllBusy)
					pthread_cond_wait(&cond, &mutex);
				bContinue = true;
				pthread_mutex_unlock(&mutex);
			}
		}
	}
}

