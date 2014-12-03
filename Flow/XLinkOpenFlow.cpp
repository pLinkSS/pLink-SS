#include <set>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include "../include/sdk.h"
#include "../include/interface.h"
#include "../include/option.h"
#include "../include/predefine.h"
#include "common.h"
#include "../PreProcess/PreProcessFactory.h"


#include "../Mass2PepIndex/Mass2PepIndex.h"
#include "../Mass2PepIndex/MetaReader.h" 
#include "../Mass2PepIndex/Mass2PepIndexReader.h"
#include "../Mass2PepIndex/DiskMass2PepIndexReader.h"
#include "../Mass2PepIndex/PeptideReaderFactory.h"
#include "../XLinkEvaluate/XLinkPeptideEvaluater.h"
#include "../XLinkScore/XLinkOpenScorer.h"
#include "../XLinkPFDIO/XLinkPFDFileIO.h"
#include "../XLinkResultReport/XLinkResultReport.h"
#include "../XLinkResultReport/XLinkResultReportFactory.h"
#include "../XLinkEvaluate/XLinkPeptideEvaluater.h"
#include "../XLinkScore/XLinkRefineScorer.h"
#include "../XLinkPFDIO/XLinkPFDFileIO.h"
#include "../PeptideGenerator/PepGeneratorFactory.h"
#include "XLinkOpenFlow.h"

using namespace std;
using namespace proteomics_sdk;
using namespace proteomics_search;
using namespace ProteinIndex;
using namespace Mass2PepIndex;




#define CANBELINKMODNAME "C-1"

CXLinkOpenFlow::CXLinkOpenFlow(vector<CSpectrum> & vSpec) :
		m_pPreProc(NULL), m_pRefinePreProc(NULL), m_pScorer(NULL), m_pRefineScorer(
		NULL), m_pEV(NULL), m_pGenerator(NULL), m_tID(0), m_tSalvoLine(0), m_tCurrentSpec(
				0), m_pReader(NULL), m_pState(NULL), m_vSpectra(vSpec), CFlow(
				vSpec)
{

	m_pTrace = NULL;
}

CXLinkOpenFlow::~CXLinkOpenFlow()
{
	Close();
}

void CXLinkOpenFlow::_InitMatchResult()
{
	try {
		CXLinkMatchResult mr;
		m_vLowMassResults.clear();
		m_vHighMassResults.clear();
		m_vResults.clear();
		size_t tTolWndNum = m_Condition.m_vPepTolWnds.size();
		m_vLowMassResults.reserve(tTolWndNum * m_vSpectra.size());
		m_vHighMassResults.reserve(tTolWndNum * m_vSpectra.size());

		m_vResultsCandidateNum.clear();
		m_vResultsCandidateNum.reserve(m_vSpectra.size());

		for (size_t i = 0; i < m_vSpectra.size(); ++i)
		{
			m_vResults.push_back(mr);
			m_vResults.back().m_vPeptideResults.reserve(m_Condition.m_nReportPep);
			m_vResultsCandidateNum.push_back(0);
			for (size_t j = 0; j < tTolWndNum; ++j)
			{
				m_vLowMassResults.push_back(vector<CXLinkOpenPepResult>());
				m_vLowMassResults.back().reserve(m_Condition.m_nSimpleReportPep);
				m_vHighMassResults.push_back(vector<CXLinkOpenPepResult>());
				m_vHighMassResults.back().reserve(m_Condition.m_nSimpleReportPep);
			}
		}
	} catch(exception &e) {
		CErrInfo info("CXLinkOpenFlow", "_InitMatchResult", "failed to reserve storage space.");
		throw runtime_error(info.Get(e));
	} catch(...) {
		CErrInfo info("CXLinkOpenFlow", "_InitMatchResult", "caught an unknown exception.");
		throw runtime_error(info.Get());
	}
}

void CXLinkOpenFlow::_CreateSpecIndex()
{
	m_vSpectraIndex.clear();
	SPEC_WND_INFO specIndex;
	m_vSpecWndCnt.clear();

	for (size_t i = 0; i < m_vSpectra.size(); ++i)
	{
		double lfMHMin, lfMHMax;
		for (size_t j = 0; j < m_Condition.m_vPepTolWnds.size(); ++j)
		{
			m_Condition.GetMHMassBorder(j, m_vSpectra[i].m_lfMH,
					m_vSpectra[i].m_nCharge, lfMHMin, lfMHMax);
			specIndex.lfMassMin = lfMHMin;
			specIndex.lfMassMax = lfMHMax;
			specIndex.nIndex = i;
			m_vSpectraIndex.push_back(specIndex);
		}
		m_vSpecWndCnt.push_back(m_Condition.m_vPepTolWnds.size());
	}

	sort(m_vSpectraIndex.begin(), m_vSpectraIndex.end(), XLINK_MIN_MASS_LESS);
}

void CXLinkOpenFlow::_CheckResult(const string & strOutputPath)
{

	CXLinkResultReportFactory reportFactory;
	CXLinkResultReport * report = reportFactory.GetReport(Xlink_STANDARD);
	vector<CXLinkMatchResult> vResults;
	vector<CSpectrum> vSpectra;
	report->Init("F:\\xlink_exp\\answer\\new_utp_b_etd1_0+-50ppm.pfind", 0);
	report->LoadFile(
			"F:\\xlink_exp\\answer\\new_utp_b_etd1_0+-50ppm_positive_qry.proteins.txt",
			vResults, vSpectra);
	report->Close();

	if (m_vSpectra.size() != vSpectra.size())
	{
		cout << "spectra size not match ! " << endl;
		return;
	}

	string strFile = strOutputPath + ".answer";

	FILE * fp;
	fp = fopen(strFile.c_str(), "w");
	if (fp == NULL)
	{
		cout << "can't open file :" << strFile << endl;
		return;
	}

	for (size_t i = 0; i < m_vSpectra.size(); ++i)
	{
		if (vResults[i].m_vPeptideResults.size() <= 0)
			continue;

		bool bFail = false;
		fprintf(fp, "name=%s\n", m_vSpectra[i].m_strFilePath.c_str());

		size_t tRank = 0;
		for (size_t j = 0; j < m_vHighMassResults[i].size(); ++j)
		{
			if (m_vHighMassResults[i][j].m_peptide.EqualPeptide(
					vResults[i].m_vPeptideResults[0].m_AlphaPeptide)
					&& m_vHighMassResults[i][j].m_nLinkSite
							== vResults[i].m_vPeptideResults[0].m_XLink.m_tAlphaSite)
			{

				tRank = j + 1;
				break;
			}
		}
		if (tRank == 0)
			bFail = true;

		fprintf(fp, "rank_alpha=%d\n", tRank);

		tRank = 0;
		if (vResults[i].m_vPeptideResults[0].m_bPair)
		{
			for (size_t j = 0; j < m_vLowMassResults[i].size(); ++j)
			{
				if (m_vLowMassResults[i][j].m_peptide.EqualPeptide(
						vResults[i].m_vPeptideResults[0].m_BetaPeptide)
						&& m_vLowMassResults[i][j].m_nLinkSite
								== vResults[i].m_vPeptideResults[0].m_XLink.m_tBetaSite)
				{
					tRank = j + 1;
					break;
				}
			}
			fprintf(fp, "rank_beta=%d\n", tRank);
		}

		if (tRank == 0)
			bFail = true;

		if (bFail)
		{
			fprintf(fp, "fail=true\n");
		}
	}

	fclose(fp);

}

void CXLinkOpenFlow::_ResultCheck()
{
	for (size_t tIdx = 0; tIdx < m_vResults.size(); ++tIdx)
	{
		for (size_t tOrder = 0;
				tOrder < m_vResults[tIdx].m_vPeptideResults.size(); ++tOrder)
		{
			if (m_vResults[tIdx].m_vPeptideResults[tOrder].m_XLink.m_eXLinkType != X_LINK)
					continue;
			CXLinkPepResult &pep = m_vResults[tIdx].m_vPeptideResults[tOrder];
			if (pep.m_AlphaPeptide.m_lfMass < pep.m_BetaPeptide.m_lfMass)
			{
				pep.SwapAlphaBeta();
			}
		}
	}
}
void CXLinkOpenFlow::_WritePFDFile(const string & strOutputPath)
{
	try
	{
		CXLinkPFDFileIO pfdIO;
		_ResultCheck();
		pfdIO.WriteFile(strOutputPath, m_vResults, m_vSpectra);
	} catch (exception & e)
	{
		CErrInfo info("_WritePFDFile", "_WritePFDFile",
				"in the function CPFDIO::WriteFile.");
		info.Append("strOutputPath=" + strOutputPath);
		throw runtime_error(info.Get(e).c_str());
	} catch (...)
	{
		CErrInfo info("_WritePFDFile", "_WritePFDFile",
				"caught an unknown exception in the function CPFDIO::WriteFile.");
		info.Append("strOutputPath=" + strOutputPath);
		throw runtime_error(info.Get().c_str());
	}

	/*
	 FILE * fp;
	 string strFileName = strOutputPath + ".temp";
	 fp = fopen(strFileName.c_str(),"w");
	 if(fp == NULL)
	 {
	 CErrInfo info("CXLinkOpenFlow", "_WritePFDFile", "in the function fopen.");
	 throw runtime_error(info.Get().c_str());
	 }

	 for(size_t i=0; i<m_vSpectra.size(); ++i)
	 {
	 fprintf(fp,"name=%s	%d\n",m_vSpectra[i].m_strFilePath.c_str(),m_vResultsCandidateNum[i]);

	 fprintf(fp,"Low Mass Results : %d\n", m_vLowMassResults[i].size());
	 for(size_t j = 0 ; j < m_vLowMassResults[i].size(); ++ j)
	 {
	 fprintf(fp,"pep%d=%d	%s	%f	",j+1,m_vLowMassResults[i][j].m_nLinkSite,m_vLowMassResults[i][j].m_peptide.m_szSequence,m_vLowMassResults[i][j].m_lfScore);

	 double lfMass = m_vLowMassResults[i][j].m_peptide.m_lfMass;
	 lfMass += _GetWaterMass();
	 lfMass += _GetProtonMass();


	 fprintf(fp,"	%f	%f	",lfMass,m_vLowMassResults[i][j].m_lfOpenMass);
	 for(int k=0;k<m_vLowMassResults[i][j].m_peptide.m_tModCnt;++k)
	 {
	 fprintf(fp,"%d,%d	",m_vLowMassResults[i][j].m_peptide.m_tModSites[k][0],m_vLowMassResults[i][j].m_peptide.m_tModSites[k][1]+1);
	 }
	 fprintf(fp,"\n");
	 }

	 fprintf(fp,"\nHigh Mass Results : %d\n", m_vHighMassResults[i].size());
	 for(size_t j = 0 ; j < m_vHighMassResults[i].size(); ++ j)
	 {
	 fprintf(fp,"pep%d=%d	%s	%f	",j+1,m_vHighMassResults[i][j].m_nLinkSite,m_vHighMassResults[i][j].m_peptide.m_szSequence,m_vHighMassResults[i][j].m_lfScore);

	 double lfMass = m_vHighMassResults[i][j].m_peptide.m_lfMass;
	 lfMass += _GetWaterMass();
	 lfMass += _GetProtonMass();


	 fprintf(fp,"	%f	%f	",lfMass,m_vHighMassResults[i][j].m_lfOpenMass);
	 for(int k=0;k<m_vHighMassResults[i][j].m_peptide.m_tModCnt;++k)
	 {
	 fprintf(fp,"%d,%d	",m_vHighMassResults[i][j].m_peptide.m_tModSites[k][0],m_vHighMassResults[i][j].m_peptide.m_tModSites[k][1]+1);
	 }
	 fprintf(fp,"\n");
	 }
	 }
	 fclose(fp);
	 */
}

void CXLinkOpenFlow::Init(const CCondition &condition, CSearchState * pState,
		size_t tID)
{
	m_pTrace = CTrace::GetInstance();

	if (!condition.ValidateAll())
	{
		CErrInfo info("CXLinkOpenFlow", "Init",
				"try CCondition::ValidateAll failed.");
		throw runtime_error(info.Get().c_str());
	}

	Close();

	m_tID = tID;
	m_pState = pState;
	m_Condition = condition;

	try
	{
		CPreProcessFactory PreProcessFactory;
		m_pPreProc = PreProcessFactory.GetPreProcessor(
				m_Condition.m_eSimplePreProcMethod);
		m_pPreProc->Init(m_Condition);
		m_pRefinePreProc = PreProcessFactory.GetPreProcessor(
				m_Condition.m_ePreProcMethod);
		m_pRefinePreProc->Init(m_Condition);
	} catch (exception & e)
	{
		CErrInfo info("CXLinkOpenFlow", "Init",
				"in the function CPreProcess::GetPreProcessor&&Init");
		throw runtime_error(info.Get(e).c_str());
	} catch (...)
	{
		CErrInfo info("CXLinkOpenFlow", "Init",
				"cauth an unknown exception, in the function CPreProcess::GetPreProcessor&&Init");
		throw runtime_error(info.Get().c_str());
	}

	try
	{
		m_pScorer = new CXLinkOpenScorer();
		m_pScorer->Initialize(m_Condition);
	} catch (exception & e)
	{
		CErrInfo info("CXLinkOpenFlow", "Init",
				"in the function CXLinkOpenScorer::Init");
		throw runtime_error(info.Get(e).c_str());
	} catch (...)
	{
		CErrInfo info("CXLinkOpenFlow", "Init",
				"cauth an unknown exception, in the function CXLinkOpenScorer::Init");
		throw runtime_error(info.Get().c_str());
	}

	try
	{
		m_pRefineScorer = new CXLinkRefineScorer();
		m_pRefineScorer->Initialize(m_Condition);
	} catch (exception & e)
	{
		CErrInfo info("CXLinkOpenFlow", "Init",
				"in the function CXLinkKSDPScorer::Init");
		throw runtime_error(info.Get(e).c_str());
	} catch (...)
	{
		CErrInfo info("CXLinkOpenFlow", "Init",
				"cauth an unknown exception, in the function CXLinkKSDPScorer::Init");
		throw runtime_error(info.Get().c_str());
	}

	try
	{
		m_pEV = new CXLinkPeptideEvaluater();
		m_pEV->Init(m_Condition);
	} catch (exception & e)
	{
		CErrInfo info("CXLinkOpenFlow", "Init",
				"in the function CEvaluater::GetEV&&Init");
		throw runtime_error(info.Get(e).c_str());
	} catch (...)
	{
		CErrInfo info("CXLinkOpenFlow", "Init",
				"cauth an unknown exception, in the function CEvaluater::GetEV&&Init");
		throw runtime_error(info.Get().c_str());
	}

	CAAConf aa(m_Condition.m_strAAListPath);
	CPepGeneratorFactory factory;
	m_pGenerator = factory.GetGenerator(m_Condition.m_ePepGen);
	m_pGenerator->Init(aa, m_Condition.m_bPepMono);
	srand(RANDOM_TIME_SEED);

	_GenMapAA2Mod();
	_CreateLinkerMap();
}

void CXLinkOpenFlow::_CreateLinkerMap()
{
	size_t i, j, k;
	for (k = 0; k < MAX_LINKER_NUM; ++k)
	{
		for (i = 0; i < LINK_SITE_NUM; ++i)
		{
			m_bMonoLinked[k][i] = false;
			for (j = 0; j < LINK_SITE_NUM; ++j)
			{
				m_bLinked[k][i][j] = false;
			}
		}
	}

	char cSite1, cSite2;
	for (k = 0; k < m_Condition.m_vSelectedXLinker.size(); ++k)
	{
		for (i = 0; i < m_Condition.m_vSelectedXLinker[k].m_strAlphaAA.length();
				++i)
		{
			cSite1 = m_Condition.m_vSelectedXLinker[k].m_strAlphaAA[i];

			if (cSite1 == PROTEIN_N_CHAR)
				cSite1 = PROTEIN_N_SITE;
			else if (cSite1 == PROTEIN_C_CHAR)
				cSite1 = PROTEIN_C_SITE;
			else if (cSite1 == PEPTIDE_N_CHAR)
				cSite1 = PEPTIDE_N_SITE;
			else if (cSite1 == PEPTIDE_C_CHAR)
				cSite1 = PEPTIDE_C_SITE;
			else if (cSite1 > 'Z' || cSite1 < 'A')
				continue;

			m_bMonoLinked[k][cSite1 - 'A'] = true;

			for (j = 0;
					j < m_Condition.m_vSelectedXLinker[k].m_strBetaAA.length();
					++j)
			{
				cSite2 = m_Condition.m_vSelectedXLinker[k].m_strBetaAA[j];

				if (cSite2 == PROTEIN_N_CHAR)
					cSite2 = PROTEIN_N_SITE;
				else if (cSite2 == PROTEIN_C_CHAR)
					cSite2 = PROTEIN_N_SITE;
				else if (cSite2 == PEPTIDE_N_CHAR)
					cSite2 = PEPTIDE_N_SITE;
				else if (cSite2 == PEPTIDE_C_CHAR)
					cSite2 = PEPTIDE_C_SITE;
				else if (cSite2 > 'Z' || cSite2 < 'A')
					continue;

				m_bMonoLinked[k][cSite2 - 'A'] = true;
				m_bLinked[k][cSite1 - 'A'][cSite2 - 'A'] = true;
				m_bLinked[k][cSite2 - 'A'][cSite1 - 'A'] = true;
			}
		}
	}
}

void CXLinkOpenFlow::Run(const string & strTmpFile)
{
	_CreateSpecIndex();
	_InitMatchResult();

	m_pTrace->Info("PreProcessing....");
	m_pState->GetProgress(m_strInfo);
	m_pTrace->Info(m_strInfo);

	try
	{
		_Preproc();
	} catch (exception & e)
	{
		CErrInfo info("CXLinkOpenFlow", "Run",
				"in the function CXLinkOpenFlow::_PreProc");
		throw runtime_error(info.Get(e).c_str());
	} catch (...)
	{
		CErrInfo info("CXLinkOpenFlow", "Run", "caught an unknown exception");
		throw runtime_error(info.Get().c_str());
	}

	m_pTrace->Info("openning database...");
	m_pState->GetProgress(m_strInfo);
	m_pTrace->Info(m_strInfo);

	_OpenDatabase();

	m_pTrace->Info("Reading the peptide sequences...");
	m_pState->GetProgress(m_strInfo);
	m_pTrace->Info(m_strInfo);

	size_t nCurrent(0);
	m_tSalvoLine = 0;
	m_tCurrentSpec = 0;
	PEP_SQ stPepSQ;
	double lfPepSQMass = 0.0;

	size_t tTotalPepNum = 0;
	tTotalPepNum = m_pReader->GetSize();

	while (m_pReader->GetNext(stPepSQ))
	{

		lfPepSQMass = stPepSQ.dfMass;
		lfPepSQMass += _GetWaterMass();
		lfPepSQMass += _GetProtonMass();

		while (m_tSalvoLine < m_vSpectraIndex.size()
				&& m_vSpectraIndex[m_tSalvoLine].lfMassMax < lfPepSQMass)
		{

			_GotoRefineStage(m_tSalvoLine);
			m_tSalvoLine++;
		}

		if (m_tSalvoLine >= m_vSpectraIndex.size())
			break;


		if (nCurrent % 100000 == 0)
		{

			m_pState->SetThreadCurrent(m_tID, m_tCurrentSpec);
			m_pState->GetProgress(m_strInfo);
			m_pTrace->Info(m_strInfo);
			char szbuf[128];
			sprintf(szbuf, "%d/%d", nCurrent, tTotalPepNum);
			m_pTrace->Info(szbuf);
		}
		++nCurrent;


		m_AssignedPep.m_peptide.SetPeptideInfor(stPepSQ.strSQ.c_str(),
				stPepSQ.strSQ.length(), stPepSQ.dfMass, stPepSQ.cMiss,
				stPepSQ.cEnd, false);
		m_AssignedPep.m_peptide.m_tModCnt = 0;
		m_AssignedPep.m_peptide.m_tFixedModCnt = 0;
		m_AssignedPep.m_lfScore = 0.0;
		m_AssignedPep.m_lfOpenMass = 0.0;
		m_AssignedPep.m_nLinkSite = -1;
		m_AssignedPep.m_nLinkerId = -1;


		_GetAllPeptide();

	}

	while (m_tSalvoLine < m_vSpectraIndex.size())
	{

		_GotoRefineStage(m_tSalvoLine);
		m_tSalvoLine++;
	}

	m_pTrace->Info("Dispatching and matching complete ...");

	m_pTrace->Info("close database...");
	_CloseDatabase();

	m_pTrace->Info("write to temporary file...");
	_WritePFDFile(strTmpFile);





	vector<vector<CXLinkOpenPepResult> >().swap(m_vLowMassResults);
	vector<vector<CXLinkOpenPepResult> >().swap(m_vHighMassResults);
	vector<CXLinkMatchResult>().swap(m_vResults);
	vector<CSpectrum>().swap(m_vSpectra);
	vector<SPEC_WND_INFO>().swap(m_vSpectraIndex);

}

void CXLinkOpenFlow::Close(void)
{
	m_Condition.clear();

	if (m_pPreProc)
		delete m_pPreProc;
	m_pPreProc = NULL;

	if (m_pScorer)
		delete m_pScorer;
	m_pScorer = NULL;

	if (m_pRefinePreProc)
	{
		delete m_pRefinePreProc;
		m_pRefinePreProc = NULL;
	}

	if (m_pRefineScorer)
	{
		delete m_pRefineScorer;
		m_pRefineScorer = NULL;
	}

	if (m_pEV)
	{
		delete m_pEV;
		m_pEV = NULL;
	}

	if (m_pGenerator)
	{
		delete m_pGenerator;
		m_pGenerator = NULL;
	}
}

void CXLinkOpenFlow::_Preproc()
{
	for (size_t i = 0; i < m_vSpectra.size(); ++i)
	{
		if (0 == m_vSpectra[i].m_tPeaksNum)
			continue;


		CSpectrum Input(m_vSpectra[i]);
		try
		{

			m_pPreProc->Run(Input, m_vSpectra[i]);
		} catch (exception & e)
		{
			CErrInfo info("CXLinkOpenFlow", "_Preproc",
					"in the function CPreProcess::Run.");
			throw runtime_error(info.Get(e).c_str());
		} catch (...)
		{
			CErrInfo info("CXLinkOpenFlow", "_Preproc",
					"caught an unknown exception in the function CPreProcess::Run.");
			throw runtime_error(info.Get().c_str());
		}

	}
}

double CXLinkOpenFlow::_Score()
{
	double lfScore = 0.0;
	try
	{
#ifdef EMPTY_FLOW
		lfScore = 2 + rand()%20;
#else
		lfScore = m_pScorer->Score();
#endif
	} catch (exception & e)
	{
		CErrInfo info("CXLinkOpenFlow", "_Score");
		char temp[200] =
		{ 0 };
		sprintf(temp, "lfScore=%.6lf", lfScore);
		info.Append(temp);
		throw runtime_error(info.Get(e).c_str());
	} catch (...)
	{
		CErrInfo info("CXLinkOpenFlow", "_Score",
				"caught an unknown exception");
		char temp[200] =
		{ 0 };
		sprintf(temp, "lfScore=%.6lf", lfScore);
		info.Append(temp);
		throw runtime_error(info.Get().c_str());
	}
	return lfScore;
}

double CXLinkOpenFlow::_Calc_Theoretical_MH(const CSimplePeptide & pep_res)
{
	if (m_Condition.m_bPepMono)
		return pep_res.m_lfMass + (IonMass_Mono_O + IonMass_Mono_H * 2.0)
				+ IonMass_Proton;
	else
		return pep_res.m_lfMass + (IonMass_Aver_O + IonMass_Aver_H * 2.0)
				+ IonMass_Aver_H;
}

void CXLinkOpenFlow::_RefineTrack(size_t tIdx)
{
	CXLinkMatchResult & mr = m_vResults[tIdx];
	const CXLinkPepResult & pep_res = m_AssignedXLinkPep;
	++mr.m_tScore;

	mr.m_tRealCandidate++;
	if (mr.m_vlfScores.size() < m_Condition.m_tMaxScoreNum)
		mr.m_vlfScores.push_back(pep_res.m_lfScore);

	if ((m_Condition.m_nReportPep <= 0) && (pep_res.m_lfScore <= 0))
	{
		return;
	}

	if (mr.m_vPeptideResults.size() < m_Condition.m_nReportPep)
	{
		mr.m_vPeptideResults.push_back(pep_res);
		mr.m_vPeptideResults.back().m_vAlphaProteinID.clear();
		mr.m_vPeptideResults.back().m_vBetaProteinID.clear();


		sort(mr.m_vPeptideResults.begin(), mr.m_vPeptideResults.end(),
				XLINK_Score_Less);

		return;
	}
	else if (pep_res.m_lfScore >= mr.m_vPeptideResults[0].m_lfScore)
	{


		mr.m_vPeptideResults[0] = pep_res;
		mr.m_vPeptideResults[0].m_vAlphaProteinID.clear();
		mr.m_vPeptideResults[0].m_vBetaProteinID.clear();

		size_t tPos = 0;

		while (1)
		{
			if ((tPos << 1) + 1 < mr.m_vPeptideResults.size()
					&& mr.m_vPeptideResults[tPos].m_lfScore
							> mr.m_vPeptideResults[(tPos << 1) + 1].m_lfScore)
			{
				if ((tPos << 1) + 2 < mr.m_vPeptideResults.size()
						&& mr.m_vPeptideResults[tPos].m_lfScore
								> mr.m_vPeptideResults[(tPos << 1) + 2].m_lfScore)
				{
					if (mr.m_vPeptideResults[(tPos << 1) + 1].m_lfScore
							< mr.m_vPeptideResults[(tPos << 1) + 2].m_lfScore)
					{
						swap(mr.m_vPeptideResults[tPos],
								mr.m_vPeptideResults[(tPos << 1) + 1]);
						tPos = (tPos << 1) + 1;
					}
					else
					{
						swap(mr.m_vPeptideResults[tPos],
								mr.m_vPeptideResults[(tPos << 1) + 2]);
						tPos = (tPos << 1) + 2;
					}

				}
				else
				{
					swap(mr.m_vPeptideResults[tPos],
							mr.m_vPeptideResults[(tPos << 1) + 1]);
					tPos = (tPos << 1) + 1;
				}
			}
			else
			{
				if ((tPos << 1) + 2 < mr.m_vPeptideResults.size()
						&& mr.m_vPeptideResults[tPos].m_lfScore
								> mr.m_vPeptideResults[(tPos << 1) + 2].m_lfScore)
				{
					swap(mr.m_vPeptideResults[tPos],
							mr.m_vPeptideResults[(tPos << 1) + 2]);
					tPos = (tPos << 1) + 2;
				}
				else
				{
					break;
				}
			}
		}

		return;
	}
}

void CXLinkOpenFlow::_RefineDispatch(size_t idx)
{
	m_pRefineScorer->SetPeptide(m_AssignedXLinkPep);
	m_pRefineScorer->SetSpectrum(m_vSpectra[m_vSpectraIndex[idx].nIndex]);
	try
	{
#ifdef EMPTY_FLOW
		m_AssignedXLinkPep.m_lfScore = rand()%100;
#else
		m_AssignedXLinkPep.m_lfScore = m_pRefineScorer->Score();
#endif

	} catch (exception & e)
	{
		CErrInfo info("CXLinkOpenFlow", "_RefineDispatch");
		char temp[200] =
		{ 0 };
		sprintf(temp, "lfScore=%.6lf", m_AssignedXLinkPep.m_lfScore);
		info.Append(temp);
		throw runtime_error(info.Get(e).c_str());
	} catch (...)
	{
		CErrInfo info("CXLinkOpenFlow", "_RefineDispatch",
				"caught an unknown exception");
		char temp[200] =
		{ 0 };
		sprintf(temp, "lfScore=%.6lf", m_AssignedXLinkPep.m_lfScore);
		info.Append(temp);
		throw runtime_error(info.Get().c_str());
	}

	_RefineTrack(m_vSpectraIndex[idx].nIndex);

}

void CXLinkOpenFlow::_GenerateRedundantPepCandidates(
		const CXLinkPepResult & modelPepCandidate,
		vector<CXLinkPepResult> & vPepCandidates, int nMinNum)
{
	vPepCandidates.clear();

	if (nMinNum <= 0)
		return;

	double lfReductMass = 0.0;

	CAAConf aa(m_Condition.m_strAAListPath);
	CMapAAMass mapAAMass = aa.GetMapAAMass();

	char cSite1 = 0, cSite2 = 0;
	double lfMass1 = 0.0;
	double lfMass2 = 0.0;

	/*
	 if(modelPepCandidate.m_XLink.m_eXLinkType == MONO_LINK)
	 {
	 lfReductMass = _GetXLinkerMass(true);
	 }
	 else if(modelPepCandidate.m_XLink.m_eXLinkType == LOOP_LINK || modelPepCandidate.m_XLink.m_eXLinkType == X_LINK)
	 {
	 lfReductMass = _GetXLinkerMass(false);
	 }
	 */

	if (modelPepCandidate.m_XLink.m_eXLinkType != NONE_LINK)
	{
		cSite1 =
				modelPepCandidate.m_AlphaPeptide.m_szSequence[modelPepCandidate.m_XLink.m_tAlphaSite];
		lfReductMass += mapAAMass.m_mapAAMass[cSite1].m_lfAvrgMass;

		if (modelPepCandidate.m_XLink.m_eXLinkType == LOOP_LINK)
		{
			cSite2 =
					modelPepCandidate.m_AlphaPeptide.m_szSequence[modelPepCandidate.m_XLink.m_tBetaSite];
			lfReductMass += mapAAMass.m_mapAAMass[cSite2].m_lfAvrgMass;
		}
		else if (modelPepCandidate.m_XLink.m_eXLinkType == X_LINK)
		{
			cSite2 =
					modelPepCandidate.m_BetaPeptide.m_szSequence[modelPepCandidate.m_XLink.m_tBetaSite];
			lfReductMass += mapAAMass.m_mapAAMass[cSite2].m_lfAvrgMass;

			lfMass1 = CSimplePeptide::Calc_Theoretical_MH(
					modelPepCandidate.m_AlphaPeptide, false);
			lfMass1 -= mapAAMass.m_mapAAMass[cSite1].m_lfAvrgMass;
			lfMass1 -= IonMass_Aver_O + IonMass_Aver_H * 3.0;

			lfMass2 = CSimplePeptide::Calc_Theoretical_MH(
					modelPepCandidate.m_BetaPeptide, false);
			lfMass2 -= mapAAMass.m_mapAAMass[cSite2].m_lfAvrgMass;
			lfMass2 -= IonMass_Aver_O + IonMass_Aver_H * 3.0;
		}
	}

	double lfMass;
	lfMass = CXLinkMatchResult::Calc_Theoretical_MH(modelPepCandidate, false)
			- lfReductMass;



	if (modelPepCandidate.m_bPair)
	{
		lfMass -= 2 * (IonMass_Aver_O + IonMass_Aver_H * 2.0) + IonMass_Aver_H;
	}
	else
	{
		lfMass -= IonMass_Aver_O + IonMass_Aver_H * 3.0;
	}

	vector<string> vAASeq;
	vector<string> vAASeq1, vAASeq2;
	int nNum1 = 0, nNum2 = 0, nNum = 0;
	nNum1 = nMinNum / 3;
	nNum2 = nMinNum / 3;
	nNum = nMinNum - nNum1 - nNum2;
	if (nNum < 0)
		nNum = 0;

	if (modelPepCandidate.m_bPair)
	{
		m_pGenerator->Batch(lfMass2, nNum1, vAASeq1);
		m_pGenerator->Batch(lfMass1, nNum2, vAASeq2);
		m_pGenerator->Batch(lfMass, nNum, vAASeq);
	}
	else
	{
		m_pGenerator->Batch(lfMass, nMinNum, vAASeq);
	}

	CXLinkPepResult pepCandidate;
	pepCandidate.m_bPair = modelPepCandidate.m_bPair;
	pepCandidate.m_XLink = modelPepCandidate.m_XLink;

	int nLen1, nLen2;
	int nRandomPos;
	string str1, str2, str3;
	char szBuf[MAX_PEPTIDE_LENGTH + 1];

	if (modelPepCandidate.m_bPair)
	{
		pepCandidate.m_AlphaPeptide = modelPepCandidate.m_AlphaPeptide;
		pepCandidate.m_XLink.m_tAlphaSite =
				modelPepCandidate.m_XLink.m_tAlphaSite;
		for (size_t i = 0; i < vAASeq1.size(); ++i)
		{

			nLen2 = vAASeq1[i].length();
			nRandomPos = rand() % (nLen2);

			pepCandidate.m_XLink.m_tBetaSite = nRandomPos;


			str1 = vAASeq1[i].substr(0, nRandomPos);
			str2 = vAASeq1[i].substr(nRandomPos, nLen2 - nRandomPos);
			sprintf(szBuf, "%s%c%s", str1.c_str(), cSite2, str2.c_str());

			pepCandidate.m_BetaPeptide.SetPeptideInfor(szBuf, nLen2 + 1, 0, 0,
					0, false);
			pepCandidate.m_BetaPeptide.m_tModCnt = 0;
			pepCandidate.m_BetaPeptide.m_tFixedModCnt = 0;
			vPepCandidates.push_back(pepCandidate);
		}

		pepCandidate.m_BetaPeptide = modelPepCandidate.m_BetaPeptide;
		pepCandidate.m_XLink.m_tBetaSite =
				modelPepCandidate.m_XLink.m_tBetaSite;
		for (size_t i = 0; i < vAASeq2.size(); ++i)
		{
			nLen1 = vAASeq2[i].length();


			nRandomPos = rand() % nLen1;

			pepCandidate.m_XLink.m_tAlphaSite = nRandomPos;


			str1 = vAASeq2[i].substr(0, nRandomPos);
			str2 = vAASeq2[i].substr(nRandomPos, nLen1 - nRandomPos);
			sprintf(szBuf, "%s%c%s", str1.c_str(), cSite1, str2.c_str());

			pepCandidate.m_AlphaPeptide.SetPeptideInfor(szBuf, nLen1 + 1, 0, 0,
					0, false);
			pepCandidate.m_AlphaPeptide.m_tModCnt = 0;
			pepCandidate.m_AlphaPeptide.m_tFixedModCnt = 0;
			vPepCandidates.push_back(pepCandidate);
		}

	}

	pepCandidate.m_AlphaPeptide.m_tFixedModCnt = 0;
	pepCandidate.m_AlphaPeptide.m_tModCnt = 0;
	pepCandidate.m_BetaPeptide.m_tFixedModCnt = 0;
	pepCandidate.m_BetaPeptide.m_tModCnt = 0;

	for (size_t i = 0; i < vAASeq.size(); ++i)
	{
		if (modelPepCandidate.m_bPair)
		{

			nLen1 = int(vAASeq[i].length() / 2);
			nLen2 = vAASeq[i].length() - nLen1;

			if (nLen1 <= 0 || nLen2 <= 0)
				continue;


			nRandomPos = rand() % nLen1;

			pepCandidate.m_XLink.m_tAlphaSite = nRandomPos;


			str1 = vAASeq[i].substr(0, nRandomPos);
			str2 = vAASeq[i].substr(nRandomPos, nLen1 - nRandomPos);
			sprintf(szBuf, "%s%c%s", str1.c_str(), cSite1, str2.c_str());
			pepCandidate.m_AlphaPeptide.SetPeptideInfor(szBuf, nLen1 + 1, 0, 0,
					0, false);


			nRandomPos = rand() % nLen2;

			pepCandidate.m_XLink.m_tBetaSite = nRandomPos;


			str1 = vAASeq[i].substr(nLen1, nRandomPos);
			str2 = vAASeq[i].substr(nLen1 + nRandomPos, nLen2 - nRandomPos);
			sprintf(szBuf, "%s%c%s", str1.c_str(), cSite2, str2.c_str());
			pepCandidate.m_BetaPeptide.SetPeptideInfor(szBuf, nLen2 + 1, 0, 0,
					0, false);
		}
		else
		{
			nLen1 = vAASeq[i].length();
			if (modelPepCandidate.m_XLink.m_eXLinkType == LOOP_LINK)
			{

				nRandomPos = rand() % nLen1;

				pepCandidate.m_XLink.m_tAlphaSite = nRandomPos;
				nRandomPos = rand() % nLen1;

				pepCandidate.m_XLink.m_tBetaSite = nRandomPos;

				if (pepCandidate.m_XLink.m_tAlphaSite
						== pepCandidate.m_XLink.m_tBetaSite)
				{
					pepCandidate.m_XLink.m_tBetaSite =
							(pepCandidate.m_XLink.m_tAlphaSite + 1) % nLen1;
				}
				if (pepCandidate.m_XLink.m_tAlphaSite
						> pepCandidate.m_XLink.m_tBetaSite)
				{
					nRandomPos = pepCandidate.m_XLink.m_tAlphaSite;
					pepCandidate.m_XLink.m_tAlphaSite =
							pepCandidate.m_XLink.m_tBetaSite;
					pepCandidate.m_XLink.m_tBetaSite = nRandomPos;
				}


				str1 = vAASeq[i].substr(0, pepCandidate.m_XLink.m_tAlphaSite);
				str2 = vAASeq[i].substr(pepCandidate.m_XLink.m_tAlphaSite,
						pepCandidate.m_XLink.m_tBetaSite
								- pepCandidate.m_XLink.m_tAlphaSite);
				str3 = vAASeq[i].substr(pepCandidate.m_XLink.m_tBetaSite,
						nLen1 - pepCandidate.m_XLink.m_tBetaSite);
				sprintf(szBuf, "%s%c%s%c%s", str1.c_str(), cSite1, str2.c_str(),
						cSite2, str3.c_str());
				pepCandidate.m_AlphaPeptide.SetPeptideInfor(szBuf, nLen1 + 2, 0,
						0, 0, false);
				pepCandidate.m_BetaPeptide.SetPeptideInfor("", 0, 0, 0, 0,
						false);
			}
			else if (modelPepCandidate.m_XLink.m_eXLinkType == MONO_LINK)
			{

				nRandomPos = rand() % nLen1;

				pepCandidate.m_XLink.m_tAlphaSite = nRandomPos;
				pepCandidate.m_XLink.m_tBetaSite = -1;


				str1 = vAASeq[i].substr(0, pepCandidate.m_XLink.m_tAlphaSite);
				str2 = vAASeq[i].substr(pepCandidate.m_XLink.m_tAlphaSite,
						nLen1 - pepCandidate.m_XLink.m_tAlphaSite);
				sprintf(szBuf, "%s%c%s", str1.c_str(), cSite1, str2.c_str());
				pepCandidate.m_AlphaPeptide.SetPeptideInfor(szBuf, nLen1 + 1, 0,
						0, 0, false);
				pepCandidate.m_BetaPeptide.SetPeptideInfor("", 0, 0, 0, 0,
						false);
			}
			else
			{
				pepCandidate.m_XLink.m_tAlphaSite = -1;
				pepCandidate.m_XLink.m_tBetaSite = -1;
				pepCandidate.m_AlphaPeptide.SetPeptideInfor(vAASeq[i].c_str(),
						nLen1, 0, 0, 0, false);
				pepCandidate.m_BetaPeptide.SetPeptideInfor("", 0, 0, 0, 0,
						false);
			}
		}

		vPepCandidates.push_back(pepCandidate);
	}

}

void CXLinkOpenFlow::_Evaluate(size_t idx)
{

	CSpectrum & spec = m_vSpectra[idx];
	CXLinkMatchResult & result = m_vResults[idx];

	if (result.m_vPeptideResults.empty())
	{
		return;
	}

	sort(result.m_vPeptideResults.begin(), result.m_vPeptideResults.end(),
			XLINK_Score_Less);
	if (m_Condition.m_tMinScoreNum > result.m_vlfScores.size())
	{
		int nNum = m_Condition.m_tMinScoreNum - result.m_vlfScores.size();
		int nAvrgNum = nNum / result.m_vPeptideResults.size() + 1;

		vector<CXLinkPepResult> vPepCandidates;
		for (size_t i = 0; i < result.m_vPeptideResults.size(); ++i)
		{
			_GenerateRedundantPepCandidates(result.m_vPeptideResults[i],
					vPepCandidates, nAvrgNum);
			for (size_t j = 0; j < vPepCandidates.size(); ++j)
			{
				m_pRefineScorer->SetSpectrum(m_vSpectra[idx]);
				m_pRefineScorer->SetPeptide(vPepCandidates[j]);
				double lfScore = 0.0;
				try
				{
					lfScore = m_pRefineScorer->Score();
				} catch (exception & e)
				{
					CErrInfo info("CXLinkOpenFlow", "_Evaluate");
					char temp[200] =
					{ 0 };
					sprintf(temp, "lfScore=%.6lf",
							m_AssignedXLinkPep.m_lfScore);
					info.Append(temp);
					throw runtime_error(info.Get(e).c_str());
				} catch (...)
				{
					CErrInfo info("CXLinkOpenFlow", "_Evaluate",
							"caught an unknown exception");
					char temp[200] =
					{ 0 };
					sprintf(temp, "lfScore=%.6lf",
							m_AssignedXLinkPep.m_lfScore);
					info.Append(temp);
					throw runtime_error(info.Get().c_str());
				}

				++result.m_tScore;
				result.m_vlfScores.push_back(lfScore);
			}
		}
	}
	/*
	 fstream fout;
	 fout.open("
	 for (int i = 0; i < result.m_vlfScores.size(); i++)
	 {
	 fout<<result.m_vlfScores[i]<<" "<<endl;
	 }
	 fout.close();
	 */
	try
	{
		if (!result.m_vlfScores.empty())
		{



			m_pEV->Run(spec, &result);
		}
	}

	catch (exception & e)
	{
		CErrInfo info("CXLinkOpenFlow", "_Evaluate",
				"in the function CEvaluate::Run.");
		char temp[200];
		info.Append("Spectrum=" + spec.m_strFilePath);
		sprintf(temp, "MH=%lf", spec.m_lfMH);
		info.Append(temp);
		sprintf(temp, "Charge=%d", spec.m_nCharge);
		info.Append(temp);
		sprintf(temp, "Score Histogram=%d", result.m_vlfScores.size());
		info.Append(temp);
		throw runtime_error(info.Get(e).c_str());
	} catch (...)
	{
		CErrInfo info("CXLinkOpenFlow", "_Evaluate",
				"caught an unknown exception in the function CEvaluate::Run.");
		char temp[200];
		info.Append("Spectrum=" + spec.m_strFilePath);
		sprintf(temp, "MH=%lf", spec.m_lfMH);
		info.Append(temp);
		sprintf(temp, "Charge=%d", spec.m_nCharge);
		info.Append(temp);
		sprintf(temp, "Score Histogram=%d", result.m_vlfScores.size());
		info.Append(temp);
		throw runtime_error(info.Get().c_str());
	}



}

























































void CXLinkOpenFlow::_GotoRefineStage(size_t idx)
{
	for (size_t i = 0; i < m_vHighMassResults[idx].size(); ++i)
	{
		if (m_vHighMassResults[idx][i].m_lfOpenMass == 0)
		{

			m_AssignedXLinkPep.m_bPair = false;
			m_AssignedXLinkPep.m_AlphaPeptide =
					m_vHighMassResults[idx][i].m_peptide;
			m_AssignedXLinkPep.m_XLink.m_eXLinkType = NONE_LINK;
			m_AssignedXLinkPep.m_XLink.m_tAlphaSite = -1;
			m_AssignedXLinkPep.m_XLink.m_tBetaSite = -1;
			m_AssignedXLinkPep.m_XLink.m_nLinkerId = -1;
			m_AssignedXLinkPep.m_BetaPeptide.SetPeptideInfor("", 0, 0.0, 0, 0,
					false);
			m_AssignedXLinkPep.m_BetaPeptide.m_tModCnt = 0;
			m_AssignedXLinkPep.m_BetaPeptide.m_tFixedModCnt = 0;

			_RefineDispatch(idx);

		}
		else if (m_vHighMassResults[idx][i].m_lfOpenMass
				== _GetXLinkerMass(true,
						m_vHighMassResults[idx][i].m_nLinkerId))
		{

			m_AssignedXLinkPep.m_bPair = false;
			m_AssignedXLinkPep.m_AlphaPeptide =
					m_vHighMassResults[idx][i].m_peptide;
			m_AssignedXLinkPep.m_XLink.m_eXLinkType = MONO_LINK;
			m_AssignedXLinkPep.m_XLink.m_tAlphaSite =
					m_vHighMassResults[idx][i].m_nLinkSite;
			m_AssignedXLinkPep.m_XLink.m_tBetaSite = -1;
			m_AssignedXLinkPep.m_XLink.m_nLinkerId =
					m_vHighMassResults[idx][i].m_nLinkerId;
			m_AssignedXLinkPep.m_BetaPeptide.SetPeptideInfor("", 0, 0.0, 0, 0,
					false);
			m_AssignedXLinkPep.m_BetaPeptide.m_tModCnt = 0;
			m_AssignedXLinkPep.m_BetaPeptide.m_tFixedModCnt = 0;

			/*bool bIllegalLink = 0;
			 for (int m = 0; m <  m_AssignedXLinkPep.m_AlphaPeptide.m_tModCnt; m++)
			 {
			 if (m_AssignedXLinkPep.m_XLink.m_tAlphaSite == m_AssignedXLinkPep.m_AlphaPeptide.m_tModSites[m][0])
			 bIllegalLink = 1;
			 }

			 if (bIllegalLink)
			 continue;
			 */
			_RefineDispatch(idx);

		}
		else if (m_vHighMassResults[idx][i].m_lfOpenMass
				== _GetXLinkerMass(false,
						m_vHighMassResults[idx][i].m_nLinkerId))
		{


			m_AssignedXLinkPep.m_bPair = false;
			m_AssignedXLinkPep.m_AlphaPeptide =
					m_vHighMassResults[idx][i].m_peptide;
			m_AssignedXLinkPep.m_XLink.m_eXLinkType = LOOP_LINK;
			m_AssignedXLinkPep.m_XLink.m_tAlphaSite =
					m_vHighMassResults[idx][i].m_nLinkSite;
			m_AssignedXLinkPep.m_XLink.m_nLinkerId =
					m_vHighMassResults[idx][i].m_nLinkerId;
			m_AssignedXLinkPep.m_BetaPeptide.SetPeptideInfor("", 0, 0.0, 0, 0,
					false);
			m_AssignedXLinkPep.m_BetaPeptide.m_tModCnt = 0;
			m_AssignedXLinkPep.m_BetaPeptide.m_tFixedModCnt = 0;

			char cSite1 = _GetCandLinkSiteIndex(
					m_AssignedXLinkPep.m_AlphaPeptide,
					m_AssignedXLinkPep.m_XLink.m_tAlphaSite,
					m_AssignedXLinkPep.m_XLink.m_nLinkerId);
			for (size_t k = m_AssignedXLinkPep.m_XLink.m_tAlphaSite + 1;
					k < m_AssignedXLinkPep.m_AlphaPeptide.m_tLength; ++k)
			{
				char cSite2 = _GetCandLinkSiteIndex(
						m_AssignedXLinkPep.m_AlphaPeptide, k,
						m_AssignedXLinkPep.m_XLink.m_nLinkerId);

				if (m_bLinked[m_vHighMassResults[idx][i].m_nLinkerId][cSite1
						- 'A'][cSite2 - 'A'])
				{
					m_AssignedXLinkPep.m_XLink.m_tBetaSite = k;
					_RefineDispatch(idx);
				}
			}
		}
		else
		{
			/*	fstream fout;
			 fout.open("
			 fout<<"Here"<<endl;
			 fout.close();
			 */

			CXLinkOpenPepResult & pep_res1 = m_vHighMassResults[idx][i];
			for (size_t j = 0; j < m_vLowMassResults[idx].size(); ++j)
			{
				CXLinkOpenPepResult & pep_res2 = m_vLowMassResults[idx][j];

				if (pep_res1.m_nLinkerId != pep_res2.m_nLinkerId)
					continue;

				double lfMH1 = _Calc_Theoretical_MH(pep_res1.m_peptide);
				double lfMH2 = _Calc_Theoretical_MH(pep_res2.m_peptide);
				double lfCalMH = lfMH1 + lfMH2
						+ _GetXLinkerMass(false, pep_res1.m_nLinkerId)
						- _GetProtonMass();

				/*				fstream fout;
				 fout.open("
				 fout<<"pep_res2 "<<pep_res2.m_peptide.m_szSequence<<" "<<pep_res2.m_lfScore<<endl;
				 fout<<"pep_res1 "<<pep_res1.m_peptide.m_szSequence<<" "<<pep_res1.m_lfScore<<endl;
				 if (pep_res2.m_peptide.m_szSequence == "VSSPVSTMMACPDGK")
				 {
				 fout<<"Here we are"<<endl;
				 if (pep_res1.m_peptide.m_szSequence == "ELCYQ")
				 fout<<"Found the wrong sequence ELCYQ."<<endl;
				 else if (pep_res1.m_peptide.m_szSequence == "CMFVR")
				 fout<<"Found the right sequence CMFVR."<<endl;
				 }
				 fout.close();
				 */
				if (lfCalMH <= m_vSpectraIndex[idx].lfMassMax
						&& lfCalMH >= m_vSpectraIndex[idx].lfMassMin)
				{
					char cSite1 = _GetCandLinkSiteIndex(pep_res1.m_peptide,
							pep_res1.m_nLinkSite, pep_res1.m_nLinkerId);

					char cSite2 = _GetCandLinkSiteIndex(pep_res2.m_peptide,
							pep_res2.m_nLinkSite, pep_res2.m_nLinkerId);

					if (m_bLinked[m_vHighMassResults[idx][i].m_nLinkerId][cSite1
							- 'A'][cSite2 - 'A'])
					{
						m_AssignedXLinkPep.m_bPair = true;
						m_AssignedXLinkPep.m_AlphaPeptide = pep_res1.m_peptide;
						m_AssignedXLinkPep.m_XLink.m_eXLinkType = X_LINK;
						m_AssignedXLinkPep.m_XLink.m_tAlphaSite =
								pep_res1.m_nLinkSite;
						m_AssignedXLinkPep.m_BetaPeptide = pep_res2.m_peptide;
						m_AssignedXLinkPep.m_XLink.m_tBetaSite =
								pep_res2.m_nLinkSite;
						m_AssignedXLinkPep.m_XLink.m_nLinkerId =
								pep_res1.m_nLinkerId;

						/*					fstream fout;
						 fout.open("
						 fout<<"pep_res2 "<<pep_res2.m_peptide.m_szSequence<<endl;
						 fout<<"pep_res1 "<<pep_res1.m_peptide.m_szSequence<<endl;
						 if (pep_res2.m_peptide.m_szSequence == "VSSPVSTMMACPDGK")
						 {
						 fout<<"Here we are to the refine dispatch."<<endl;
						 if (pep_res1.m_peptide.m_szSequence == "ELCYQ")
						 fout<<"Found the wrong sequence ELCYQ."<<endl;
						 else if (pep_res1.m_peptide.m_szSequence == "CMFVR")
						 fout<<"Found the right sequence CMFVR."<<endl;
						 }
						 fout.close();
						 */
						/*bool bIllegalLink = 0;
						 for (int m = 0; m <  m_AssignedXLinkPep.m_AlphaPeptide.m_tModCnt; m++)
						 {
						 if (m_AssignedXLinkPep.m_XLink.m_tAlphaSite == m_AssignedXLinkPep.m_AlphaPeptide.m_tModSites[m][0] || m_AssignedXLinkPep.m_XLink.m_tBetaSite == m_AssignedXLinkPep.m_BetaPeptide.m_tModSites[m][0])
						 bIllegalLink = 1;
						 }

						 if (bIllegalLink)
						 continue;
						 */
						_RefineDispatch(idx);
					}
				}
			}
		}
	}

	m_vHighMassResults[idx].clear();
	m_vLowMassResults[idx].clear();
	vector<CXLinkOpenPepResult>().swap(m_vLowMassResults[idx]);
	vector<CXLinkOpenPepResult>().swap(m_vHighMassResults[idx]);

	m_vSpecWndCnt[m_vSpectraIndex[idx].nIndex]--;
	if (m_vSpecWndCnt[m_vSpectraIndex[idx].nIndex] == 0)
	{
		_Evaluate(m_vSpectraIndex[idx].nIndex);
		m_tCurrentSpec++;
		vector<double>().swap(
				m_vResults[m_vSpectraIndex[idx].nIndex].m_vlfScores);
	}

}

void CXLinkOpenFlow::_GenMapAA2Mod()
{
	m_mapFixMod.clear();
	m_mapVarMod.clear();

	bool bPepMono = m_Condition.m_bPepMono;
	size_t i = 0, j = 0;
	CSimpleMod smod;

	vector<CModification> &vVarMod = m_Condition.m_vSelectedVarMod;
	for (i = 0; i < vVarMod.size(); ++i)
	{
		smod.m_nIdx = i;
		smod.m_eModType = vVarMod[i].m_eModType;
		CModification & mod = vVarMod[i];
		if (bPepMono)
		{
			smod.m_lfMass = vVarMod[i].m_lfMonoMass_dif;
		}
		else
		{
			smod.m_lfMass = vVarMod[i].m_lfAvrgMass_dif;
		}

		if (MT_PEP_NTERM == mod.m_eModType || MT_PRO_NTERM == mod.m_eModType)
		{
			m_mapVarMod[0].push_back(smod);
		}
		else if (MT_PEP_CTERM == mod.m_eModType
				|| MT_PRO_CTERM == mod.m_eModType)
		{
			m_mapVarMod[1].push_back(smod);
		}
		else
		{
			for (j = 0; j < mod.m_strAA.length(); ++j)
			{
				m_mapVarMod[mod.m_strAA[j]].push_back(smod);
			}
		}
	}

	vector<CModification> &vFixMod = m_Condition.m_vSelectedFixMod;
	for (i = 0; i < vFixMod.size(); ++i)
	{
		smod.m_nIdx = i + vVarMod.size();
		smod.m_eModType = vFixMod[i].m_eModType;
		CModification & mod = vFixMod[i];
		if (bPepMono)
		{
			smod.m_lfMass = vFixMod[i].m_lfMonoMass_dif;
		}
		else
		{
			smod.m_lfMass = vFixMod[i].m_lfAvrgMass_dif;
		}

		if (MT_PEP_NTERM == mod.m_eModType || MT_PRO_NTERM == mod.m_eModType)
		{
			m_mapFixMod[0].push_back(smod);
		}
		else if (MT_PEP_CTERM == mod.m_eModType
				|| MT_PRO_CTERM == mod.m_eModType)
		{
			m_mapFixMod[1].push_back(smod);
		}
		else
		{
			for (j = 0; j < mod.m_strAA.length(); ++j)
				m_mapFixMod[mod.m_strAA[j]].push_back(smod);
		}
	}
}

void CXLinkOpenFlow::_OpenDatabase()
{
	_CloseDatabase();

	string strDBFilePath;
	CPeptideRandomReaderFactory factory;
	if (MEMORY_MAPPED == m_Condition.m_eIndex)
		m_pReader = factory.GetPeptideReader(Peptide_FileType_Index_Map);
	else
		m_pReader = factory.GetPeptideReader(Peptide_FileType_Index_Disk);

	string strMeta, strDBPath;

	try
	{
		CDBConf db(m_Condition.m_strDBConfPath.c_str());
		strMeta = db.GetMetaName(m_Condition.m_vSelectedDBName[0].c_str(),
				m_Condition.m_SelectedEnzyme.m_strName.c_str(),
				m_Condition.m_bPepMono);
		if (strMeta == "NULL")
		{
			CErrInfo info("CXLinkOpenFlow", "_OpenDatabase", "strMeta==NULL");
			info.Append("dbconfpath=" + m_Condition.m_strDBConfPath);
			info.Append("selecteddbname=" + m_Condition.m_vSelectedDBName[0]);
			info.Append(
					"selectedenzyme=" + m_Condition.m_SelectedEnzyme.m_strName);
			throw runtime_error(info.Get().c_str());
		}

		strDBPath = db.GetPath(m_Condition.m_vSelectedDBName[0].c_str(),
				m_Condition.m_SelectedEnzyme.m_strName.c_str());

		if (MEMORY_MAPPED == m_Condition.m_eIndex)
			m_pReader->Open(strDBPath, strMeta, Protein_Random_Index_Map,
					false);
		else
			m_pReader->Open(strDBPath, strMeta, Protein_Random_Index_Disk,
					false);
	} catch (exception & e)
	{
		CErrInfo info("CXLinkOpenFlow", "_OpenDatabase",
				"in the initialize of the database.");
		info.Append("strMeta=" + strMeta);
		info.Append("strDBPath=" + strDBPath);
		throw runtime_error(info.Get(e).c_str());
	} catch (...)
	{
		CErrInfo info("CXLinkOpenFlow", "_OpenDatabase",
				"caught an unknown exception in the initialize of the database.");
		info.Append("strMeta=" + strMeta);
		info.Append("strDBPath=" + strDBPath);
		throw runtime_error(info.Get().c_str());
	}
}

void CXLinkOpenFlow::_CloseDatabase()
{
	if (m_pReader == NULL)
		return;

	try
	{
		m_pReader->Close();
	} catch (exception & e)
	{
		CErrInfo info("CXLinkOpenFlow", "_CloseDatabase",
				"in the function CPeptideReader::Close.");
		throw runtime_error(info.Get(e).c_str());
	} catch (...)
	{
		CErrInfo info("CXLinkOpenFlow", "_CloseDatabase",
				"caught an unknown exception in the function CPeptideReader::Close.");
		throw runtime_error(info.Get().c_str());
	}

	delete m_pReader;
	m_pReader = NULL;
}

void CXLinkOpenFlow::_GetAllPeptide()
{

	_AttatchFixedMod();


	_AttatchVarMod(-1);
}

void CXLinkOpenFlow::_AttatchFixedMod()
{
	CSimplePeptide * pPeptide = NULL;
	pPeptide = &m_AssignedPep.m_peptide;

	if (pPeptide->m_tLength <= 0)
		return;

	pPeptide->m_tFixedModCnt = 0;
	pPeptide->m_tModCnt = 0;
	size_t i = 0;


	for (i = 0; i < m_mapFixMod[0].size(); ++i)
	{
		if (m_mapFixMod[0][i].m_eModType == MT_PRO_NTERM
				&& !pPeptide->IsProNTerm())
			continue;


		if (m_Condition.m_vSelectedFixMod[m_mapFixMod[0][i].m_nIdx
				- m_Condition.m_vSelectedVarMod.size()].m_strAA.find(
				pPeptide->m_szSequence[0]) == string::npos)
			continue;

		pPeptide->m_tModSites[pPeptide->m_tModCnt][0] = 0;
		pPeptide->m_tModSites[pPeptide->m_tModCnt][1] =
				m_mapFixMod[0][i].m_nIdx;
		pPeptide->m_lfMass += m_mapFixMod[0][i].m_lfMass;
		++pPeptide->m_tModCnt;
	}


	for (i = 0; i < m_mapFixMod[1].size(); ++i)
	{
		if (m_mapFixMod[1][i].m_eModType == MT_PRO_CTERM
				&& !pPeptide->IsProCTerm())
			continue;

		if (m_Condition.m_vSelectedFixMod[m_mapFixMod[1][i].m_nIdx
				- m_Condition.m_vSelectedVarMod.size()].m_strAA.find(
				pPeptide->m_szSequence[pPeptide->m_tLength - 1])
				== string::npos)
			continue;

		pPeptide->m_tModSites[pPeptide->m_tModCnt][0] = pPeptide->m_tLength - 1;
		pPeptide->m_tModSites[pPeptide->m_tModCnt][1] =
				m_mapFixMod[1][i].m_nIdx;
		pPeptide->m_lfMass += m_mapFixMod[1][i].m_lfMass;
		++pPeptide->m_tModCnt;
	}


	for (i = 0; i < pPeptide->m_tLength; ++i)
	{
		char cSite = pPeptide->m_szSequence[i];
		for (size_t j = 0; j < m_mapFixMod[cSite].size(); ++j)
		{
			pPeptide->m_tModSites[pPeptide->m_tModCnt][0] = i;
			pPeptide->m_tModSites[pPeptide->m_tModCnt][1] =
					m_mapFixMod[cSite][j].m_nIdx;
			pPeptide->m_lfMass += m_mapFixMod[cSite][j].m_lfMass;
			++pPeptide->m_tModCnt;
		}
	}
	pPeptide->m_tFixedModCnt = pPeptide->m_tModCnt;
}

void CXLinkOpenFlow::_AttatchVarMod(int nCurrLen)
{
	CSimplePeptide * pPeptide = NULL;
	pPeptide = &m_AssignedPep.m_peptide;

	if ((int) (pPeptide->m_tModCnt - pPeptide->m_tFixedModCnt)
			> m_Condition.m_nMaxModifyNumber)
		return;

	if ((int) pPeptide->m_tLength <= 0)
		return;

	if (nCurrLen == -1)
	{


		_AttatchVarMod(0);

		if ((int) (pPeptide->m_tModCnt - pPeptide->m_tFixedModCnt)
				< m_Condition.m_nMaxModifyNumber)
		{
			for (size_t w = 0; w < m_mapVarMod[0].size(); ++w)
			{
				if (m_mapVarMod[0][w].m_eModType == MT_PRO_NTERM
						&& !pPeptide->IsProNTerm())
					continue;


				if (m_Condition.m_vSelectedVarMod[m_mapVarMod[0][w].m_nIdx].m_strAA.find(
						pPeptide->m_szSequence[0]) == string::npos)
					continue;

				pPeptide->m_tModSites[pPeptide->m_tModCnt][0] = 0;
				pPeptide->m_tModSites[pPeptide->m_tModCnt][1] =
						m_mapVarMod[0][w].m_nIdx;

				pPeptide->m_lfMass += m_mapVarMod[0][w].m_lfMass;
				++pPeptide->m_tModCnt;

				_AttatchVarMod(0);

				pPeptide->m_lfMass -= m_mapVarMod[0][w].m_lfMass;
				--pPeptide->m_tModCnt;

			}
		}

	}
	else if (nCurrLen == (int) pPeptide->m_tLength)
	{


		_Dispatch();

		if ((int) (pPeptide->m_tModCnt - pPeptide->m_tFixedModCnt)
				< m_Condition.m_nMaxModifyNumber)
		{
			for (size_t i = 0; i < m_mapVarMod[1].size(); ++i)
			{
				if (m_mapVarMod[1][i].m_eModType == MT_PRO_CTERM
						&& !pPeptide->IsProCTerm())
					continue;


				if (m_Condition.m_vSelectedVarMod[m_mapVarMod[1][i].m_nIdx].m_strAA.find(
						pPeptide->m_szSequence[nCurrLen - 1]) == string::npos)
					continue;

				pPeptide->m_tModSites[pPeptide->m_tModCnt][0] = nCurrLen - 1;
				pPeptide->m_tModSites[pPeptide->m_tModCnt][1] =
						m_mapVarMod[1][i].m_nIdx;

				pPeptide->m_lfMass += m_mapVarMod[1][i].m_lfMass;
				++pPeptide->m_tModCnt;

				_Dispatch();

				pPeptide->m_lfMass -= m_mapVarMod[1][i].m_lfMass;
				--pPeptide->m_tModCnt;
			}
			return;
		}
	}
	else
	{
		while (nCurrLen < (int) pPeptide->m_tLength
				&& m_mapVarMod.find(pPeptide->m_szSequence[nCurrLen])
						== m_mapVarMod.end())
		{
			++nCurrLen;
		}

		if (nCurrLen == (int) pPeptide->m_tLength)
			_AttatchVarMod(nCurrLen);
		else
		{


			_AttatchVarMod(nCurrLen + 1);

			char cSite = pPeptide->m_szSequence[nCurrLen];

			if ((int) (pPeptide->m_tModCnt - pPeptide->m_tFixedModCnt)
					< m_Condition.m_nMaxModifyNumber)
			{
				for (size_t i = 0; i < m_mapVarMod[cSite].size(); ++i)
				{
					pPeptide->m_tModSites[pPeptide->m_tModCnt][0] = nCurrLen;
					pPeptide->m_tModSites[pPeptide->m_tModCnt][1] =
							m_mapVarMod[cSite][i].m_nIdx;

					pPeptide->m_lfMass += m_mapVarMod[cSite][i].m_lfMass;
					++pPeptide->m_tModCnt;

					_AttatchVarMod(nCurrLen + 1);

					pPeptide->m_lfMass -= m_mapVarMod[cSite][i].m_lfMass;
					--pPeptide->m_tModCnt;
				}

				return;
			}
		}
	}

}

double CXLinkOpenFlow::_Calc_Theoretical_MH()
{
	if (m_Condition.m_bPepMono)
		return m_AssignedPep.m_peptide.m_lfMass
				+ (IonMass_Mono_O + IonMass_Mono_H * 2.0) + IonMass_Proton;
	else
		return m_AssignedPep.m_peptide.m_lfMass
				+ (IonMass_Aver_O + IonMass_Aver_H * 2.0) + IonMass_Aver_H;
}

double CXLinkOpenFlow::_GetWaterMass()
{
	if (m_Condition.m_bPepMono)
	{
		return IonMass_Mono_O + IonMass_Mono_H * 2.0;
	}
	else
	{
		return IonMass_Aver_O + IonMass_Aver_H * 2.0;
	}
}
double CXLinkOpenFlow::_GetProtonMass()
{
	if (m_Condition.m_bPepMono)
	{
		return IonMass_Proton;
	}
	else
	{
		return IonMass_Aver_H;
	}
}

double CXLinkOpenFlow::_GetXLinkerMass(bool bMono, int nLinkerId)
{

	if (bMono)
	{
		if (m_Condition.m_bPepMono)
		{
			return m_Condition.m_vSelectedXLinker[nLinkerId].m_lfMLMonoMass_dif;
		}
		else
		{
			return m_Condition.m_vSelectedXLinker[nLinkerId].m_lfMLAvrgMass_dif;
		}
	}
	else
	{
		if (m_Condition.m_bPepMono)
		{
			return m_Condition.m_vSelectedXLinker[nLinkerId].m_lfMonoMass_dif;
		}
		else
		{
			return m_Condition.m_vSelectedXLinker[nLinkerId].m_lfAvrgMass_dif;
		}
	}
}

bool CXLinkOpenFlow::_IsXLinked(const CSimplePeptide &peptide)
{
	size_t i, j;
	for (j = 0; j < m_Condition.m_vSelectedXLinker.size(); ++j)
	{
		for (i = 0; i < peptide.m_tLength; ++i)
		{
			if (m_bMonoLinked[j][peptide.m_szSequence[i] - 'A'])
			{
				return true;
			}
		}

		if (m_bMonoLinked[j][PEPTIDE_N_SITE - 'A'])
		{
			return true;
		}
		if (m_bMonoLinked[j][PEPTIDE_C_SITE - 'A'])
		{
			return true;
		}


		if (peptide.IsProNTerm() && m_bMonoLinked[j][PROTEIN_N_SITE - 'A'])
		{
			return true;
		}

		if (peptide.IsProCTerm() && m_bMonoLinked[j][PROTEIN_C_SITE - 'A'])
		{
			return true;
		}
	}
	return false;
}

bool CXLinkOpenFlow::IsModLinkSiteConflict(const CSimplePeptide &peptide,
		const int nIdx)
{
#ifndef MODSITE_LINKABLE
	for (size_t i = 0; i < peptide.m_tModCnt; ++i)
	{

		int iIdx = peptide.m_tModSites[i][1] - m_Condition.m_vSelectedVarMod.size();
		if ( iIdx >=0 && iIdx < m_Condition.m_vSelectedFixMod.size() )
		{
			if (m_Condition.m_vSelectedFixMod[iIdx].m_strName == CANBELINKMODNAME)
				continue;
		}

		if (peptide.m_tModSites[i][0] == nIdx)
			return true;
	}
#endif
	return false;
}

bool CXLinkOpenFlow::IsLinkSiteCleave(const CSimplePeptide &peptide, const int nLinkSite)
{
#ifndef CLEAVESITE_LINKABLE
	if (nLinkSite >= peptide.m_tLength)
		cerr<<"Error found in CXLinkPepSalvoFlow::IsLinkSiteCleave."<<endl;

	string strCleaveString =m_Condition.m_SelectedEnzyme.GetCleaveString();
	if (m_Condition.m_SelectedEnzyme.m_bNTerm)
	{
		if (nLinkSite != 0)
			return false;
		if (peptide.IsProNTerm())
			return false;
		if (strCleaveString.find(peptide.m_szSequence[nLinkSite]) != -1)
			return true;
	}
	else
	{
		if (nLinkSite != peptide.m_tLength-1)
			return false;
		if (peptide.IsProCTerm())
			return false;
		if (strCleaveString.find(peptide.m_szSequence[nLinkSite]) != -1)
			return true;
	}
#endif

	return false;
}

char CXLinkOpenFlow::_GetCandLinkSiteIndex(const CSimplePeptide &peptide,
		int nIdx, int nLinkerId)
{
	if (nIdx == 0)
	{
		if (peptide.IsProNTerm()
				&& m_bMonoLinked[nLinkerId][PROTEIN_N_SITE - 'A'])
		{
			return PROTEIN_N_SITE;
		}
		else if (m_bMonoLinked[nLinkerId][PEPTIDE_N_SITE - 'A'])
		{
			return PEPTIDE_N_SITE;
		}
		else
		{
			return peptide.m_szSequence[nIdx];
		}
	}
	else if (nIdx == peptide.m_tLength - 1)
	{
		if (peptide.IsProCTerm()
				&& m_bMonoLinked[nLinkerId][PROTEIN_C_SITE - 'A'])
		{
			return PROTEIN_C_SITE;
		}
		else if (m_bMonoLinked[nLinkerId][PEPTIDE_C_SITE - 'A'])
		{
			return PEPTIDE_C_SITE;
		}
		else
		{
			return peptide.m_szSequence[nIdx];
		}
	}
	else
	{
		return peptide.m_szSequence[nIdx];
	}
}

void CXLinkOpenFlow::_Dispatch()
{


	double lfMass = _Calc_Theoretical_MH();
	double lfDeltaMass = 0.0;
	int nMassRegion = 1;

	size_t nId;
	for (nId = m_tSalvoLine; nId < m_vSpectraIndex.size(); ++nId)
	{
		if (lfMass > m_vSpectraIndex[nId].lfMassMax)
			continue;
		else
			break;
	}

	if (nId < m_vSpectraIndex.size())
	{

		m_pScorer->SetPeptide(m_AssignedPep);
	}

	for (; nId < m_vSpectraIndex.size(); ++nId)
	{
		m_pScorer->SetSpectrum(m_vSpectra[m_vSpectraIndex[nId].nIndex]);

		double lfWndMass = (m_vSpectraIndex[nId].lfMassMax
				+ m_vSpectraIndex[nId].lfMassMin) / 2;
		lfDeltaMass = lfWndMass - lfMass;

		if (lfMass <= m_vSpectraIndex[nId].lfMassMax
				&& lfMass >= m_vSpectraIndex[nId].lfMassMin)
		{

			lfDeltaMass = 0;

			m_AssignedPep.m_lfOpenMass = lfDeltaMass;
			m_AssignedPep.m_nLinkSite = -1;
			m_AssignedPep.m_nLinkerId = -1;
			m_AssignedPep.m_lfScore = _Score();

			if (m_AssignedPep.m_lfScore > m_Condition.m_lfSimpleScoreCutoff)
			{
				m_vResultsCandidateNum[m_vSpectraIndex[nId].nIndex]++;
				_track(nId, false);
			}
		}
		else
		{
			if (!_IsXLinked(m_AssignedPep.m_peptide))
				return;

			for (size_t nLinkerId = 0;
					nLinkerId < m_Condition.m_vSelectedXLinker.size();
					++nLinkerId)
			{
				m_AssignedPep.m_nLinkerId = nLinkerId;

				if (lfMass + _GetXLinkerMass(true, nLinkerId)
						<= m_vSpectraIndex[nId].lfMassMax
						&& lfMass + _GetXLinkerMass(true, nLinkerId)
								>= m_vSpectraIndex[nId].lfMassMin)
				{

					nMassRegion = 1;
					lfDeltaMass = _GetXLinkerMass(true, nLinkerId);
				}
				else if (lfMass + _GetXLinkerMass(false, nLinkerId)
						<= m_vSpectraIndex[nId].lfMassMax
						&& lfMass + _GetXLinkerMass(false, nLinkerId)
								>= m_vSpectraIndex[nId].lfMassMin)
				{

					nMassRegion = 1;
					lfDeltaMass = _GetXLinkerMass(false, nLinkerId);
				}
				else
				{

					if (lfMass + MIN_PEPTIDE_MASS
							+ _GetXLinkerMass(false, nLinkerId)
							> m_vSpectraIndex[nId].lfMassMax)

						continue;

					double lfMidMassLB = (m_vSpectraIndex[nId].lfMassMin
							- _GetXLinkerMass(false, nLinkerId)
							- _GetProtonMass()) / 2;
					double lfMidMassUB = (m_vSpectraIndex[nId].lfMassMax
							- _GetXLinkerMass(false, nLinkerId)
							- _GetProtonMass()) / 2;
					double lfNeutralMass = lfMass - _GetProtonMass();

					if (lfNeutralMass < lfMidMassLB) /* this is beta */
						nMassRegion = -1;
					else if (lfNeutralMass > lfMidMassUB) /* this is alpha */
						nMassRegion = 1;
					else
						/* two peptides almost have the same mass */
						nMassRegion = 0;
				}

				m_AssignedPep.m_lfOpenMass = lfDeltaMass;

				for (size_t i = 0; i < m_AssignedPep.m_peptide.m_tLength; ++i)
				{
					char cSite = _GetCandLinkSiteIndex(m_AssignedPep.m_peptide,
							i, nLinkerId);

					if (m_bMonoLinked[nLinkerId][cSite - 'A'])
					{
						if (IsModLinkSiteConflict(m_AssignedPep.m_peptide, i))
							continue;
						if (IsLinkSiteCleave(m_AssignedPep.m_peptide, i))
							continue;
						m_AssignedPep.m_nLinkSite = i;
						m_AssignedPep.m_lfScore = _Score();

						if (m_AssignedPep.m_lfScore
								> m_Condition.m_lfSimpleScoreCutoff)
						{
							m_vResultsCandidateNum[m_vSpectraIndex[nId].nIndex]++;
							if (nMassRegion <= 0)
								_track(nId, true);
							if (nMassRegion >= 0)
								_track(nId, false);
						}
					}
				}
			}
		}
	}
}

void CXLinkOpenFlow::_track(size_t tIdx, bool bLowOrHigh)
{
#ifdef _DEBUG2
	sprintf(m_debug.m_szBuf,"%s	%d",m_AssignedPep.m_peptide.m_szSequence,bLowOrHigh);
	m_debug.WriteLine();
#endif

	/*	fstream fout;
	 fout.open("
	 fout<<"m_AssignedPep.m_peptide.m_szSequence "<<m_AssignedPep.m_peptide.m_szSequence<<" "<<m_AssignedPep.m_lfScore<<endl;
	 fout.close();
	 */

	/*bool bIllegalLink = 0;
	for (int m = 0; m < m_AssignedPep.m_peptide.m_tModCnt; m++)
	{
		if (m_AssignedPep.m_nLinkSite
				== m_AssignedPep.m_peptide.m_tModSites[m][0])
			bIllegalLink = 1;
	}

	if (bIllegalLink)
		return;*/

	const CXLinkOpenPepResult & pep_res = m_AssignedPep;

	vector<vector<CXLinkOpenPepResult> > *pResults = NULL;

	if (bLowOrHigh)
		pResults = &m_vLowMassResults;
	else
		pResults = &m_vHighMassResults;

	if ((*pResults)[tIdx].size() < m_Condition.m_nSimpleReportPep)
	{
		(*pResults)[tIdx].push_back(pep_res);


		sort((*pResults)[tIdx].begin(), (*pResults)[tIdx].end(),
				XLINK_OPENFLOW_Score_Less);
		return;
	}
	else if (pep_res.m_lfScore >= (*pResults)[tIdx][0].m_lfScore)
	{



		(*pResults)[tIdx][0] = pep_res;

		size_t tPos = 0;

		while (1)
		{
			if ((tPos << 1) + 1 < (*pResults)[tIdx].size()
					&& (*pResults)[tIdx][tPos].m_lfScore
							> (*pResults)[tIdx][(tPos << 1) + 1].m_lfScore)
			{
				if ((tPos << 1) + 2 < (*pResults)[tIdx].size()
						&& (*pResults)[tIdx][tPos].m_lfScore
								> (*pResults)[tIdx][(tPos << 1) + 2].m_lfScore)
				{
					if ((*pResults)[tIdx][(tPos << 1) + 1].m_lfScore
							< (*pResults)[tIdx][(tPos << 1) + 2].m_lfScore)
					{
						swap((*pResults)[tIdx][tPos],
								(*pResults)[tIdx][(tPos << 1) + 1]);
						tPos = (tPos << 1) + 1;
					}
					else
					{
						swap((*pResults)[tIdx][tPos],
								(*pResults)[tIdx][(tPos << 1) + 2]);
						tPos = (tPos << 1) + 2;
					}

				}
				else
				{
					swap((*pResults)[tIdx][tPos],
							(*pResults)[tIdx][(tPos << 1) + 1]);
					tPos = (tPos << 1) + 1;
				}
			}
			else
			{
				if ((tPos << 1) + 2 < (*pResults)[tIdx].size()
						&& (*pResults)[tIdx][tPos].m_lfScore
								> (*pResults)[tIdx][(tPos << 1) + 2].m_lfScore)
				{
					swap((*pResults)[tIdx][tPos],
							(*pResults)[tIdx][(tPos << 1) + 2]);
					tPos = (tPos << 1) + 2;
				}
				else
				{
					break;
				}
			}
		}
		return;

	}

}

