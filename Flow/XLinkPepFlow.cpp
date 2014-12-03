#define CANBELINKMODNAME "C-1"
#define PEPMINLEN 4
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
#include "SpectraSearcher.h"
#include "../PreProcess/PreProcessFactory.h"
#include "../ProteinInfer/ProteinInferFactory.h"
#include "../ProteinIndex/ProteinIndex.h"
#include "../ProteinIndex/ProteinReader.h"
#include "../ProteinIndex/DiskRandomSQReader.h"
#include "../ProteinIndex/ProteinReaderFactory.h"
#include "../Mass2PepIndex/Mass2PepIndex.h"
#include "../Mass2PepIndex/MetaReader.h"
#include "../Mass2PepIndex/Mass2PepIndexReader.h"
#include "../Mass2PepIndex/DiskMass2PepIndexReader.h"
#include "../Mass2PepIndex/PeptideReaderFactory.h"
#include "../XLinkEvaluate/XLinkPeptideEvaluater.h"

#include "../XLinkScore/XLinkRefineScorer.h"
#include "../XLinkPFDIO/XLinkPFDFileIO.h"
#include "../PeptideGenerator/PepGeneratorFactory.h"
/*pep list*/
#include "../PepPairContainer/common.h"
#include "../PepPairContainer/PepPairContainer.h"
#include "PepPairList.h"
#include "XLinkPepFlow.h"

using namespace std;
using namespace proteomics_sdk;
using namespace proteomics_search;
using namespace ProteinIndex;
using namespace Mass2PepIndex;



bool CXLinkPepSalvoFlow::m_bPepPairListCreated = false;
/*pep list*/
#ifdef _USE_PEP_LIST
CPepPairList CXLinkPepSalvoFlow::m_CandidatePepPairs = CPepPairList();
#else
CPepPairContainer CXLinkPepSalvoFlow::m_CandidatePepPairs = CPepPairContainer();
#endif




CXLinkPepSalvoFlow::CXLinkPepSalvoFlow(vector<CSpectrum> & vSpec) :
		m_pPreProc(NULL), m_pScorer(NULL), m_pEV(NULL), m_pGenerator(NULL), m_tID(
				0), m_tSalvoLine(0), m_pReader(NULL), m_pState(NULL), m_vSpectra(
				vSpec), CFlow(vSpec)
{
}

CXLinkPepSalvoFlow::~CXLinkPepSalvoFlow(void)
{
	Close();
}

void CXLinkPepSalvoFlow::Init(const CCondition &condition,
		CSearchState * pState, size_t tID)
{
	m_pTrace = CTrace::GetInstance();
	m_pTrace->Debug("in CXLinkPepSalvoFlow::Init()");

	if (!condition.ValidateAll())
	{
		CErrInfo info("CXLinkPepSalvoFlow", "Init",
				"try CCondition::ValidateAll failed.");
		throw runtime_error(info.Get().c_str());
	}

	Close();

#ifdef _DEBUG2
	char szbuf[256];
	sprintf(szbuf,"xlinkdebug_%d.txt",tID);
	m_debug.Init(szbuf);
	m_debug.Open();
#endif

	m_tID = tID;
	m_pState = pState;
	m_Condition = condition;

	try
	{
		CPreProcessFactory PreProcessFactory;
		m_pPreProc = PreProcessFactory.GetPreProcessor(
				m_Condition.m_ePreProcMethod);
		m_pPreProc->Init(m_Condition);
	} catch (exception & e)
	{
		CErrInfo info("CXLinkPepSalvoFlow", "Init",
				"in the function CPreProcess::GetPreProcessor&&Init");
		throw runtime_error(info.Get(e).c_str());
	} catch (...)
	{
		CErrInfo info("CXLinkPepSalvoFlow", "Init",
				"cauth an unknown exception, in the function CPreProcess::GetPreProcessor&&Init");
		throw runtime_error(info.Get().c_str());
	}

	try
	{
		m_pScorer = new CXLinkRefineScorer();
		m_pScorer->Initialize(m_Condition);

	} catch (exception & e)
	{
		CErrInfo info("CXLinkPepSalvoFlow", "Init",
				"in the function CXLinkKSDPScorer::Init");
		throw runtime_error(info.Get(e).c_str());
	} catch (...)
	{
		CErrInfo info("CXLinkPepSalvoFlow", "Init",
				"cauth an unknown exception, in the function CXLinkRefineScorer::Init");
		throw runtime_error(info.Get().c_str());
	}

	try
	{
		m_pEV = new CXLinkPeptideEvaluater();
		m_pEV->Init(m_Condition);
	} catch (exception & e)
	{
		CErrInfo info("CXLinkPepSalvoFlow", "Init",
				"in the function CEvaluater::GetEV&&Init");
		throw runtime_error(info.Get(e).c_str());
	} catch (...)
	{
		CErrInfo info("CXLinkPepSalvoFlow", "Init",
				"cauth an unknown exception, in the function CEvaluater::GetEV&&Init");
		throw runtime_error(info.Get().c_str());
	}

	_GenMapAA2Mod();
	_CreateLinkerMap();

	CAAConf aa(m_Condition.m_strAAListPath);
	CPepGeneratorFactory factory;
	m_pGenerator = factory.GetGenerator(m_Condition.m_ePepGen);
	m_pGenerator->Init(aa, m_Condition.m_bPepMono);

	CMapAAMass mapAAMass = aa.GetMapAAMass();


	srand(RANDOM_TIME_SEED);

}

void CXLinkPepSalvoFlow::_GenMapAA2Mod()
{
	m_mapFixMod.clear();
	m_mapVarMod.clear();

	bool bPepMono = m_Condition.m_bPepMono;
	size_t i = 0, j = 0;
	vector<CModification> & vVarMod = m_Condition.m_vSelectedVarMod;

	CSimpleMod smod;
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
	vector<CModification> & vFixMod = m_Condition.m_vSelectedFixMod;
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

void CXLinkPepSalvoFlow::_AttatchVarMod(int nPeptideId, int nCurrLen)
{
	CSimplePeptide * pPeptide;
	if (nPeptideId == 0)
	{
		pPeptide = &m_AssignedPep.m_AlphaPeptide;
	}
	else
	{
		pPeptide = &m_AssignedPep.m_BetaPeptide;
	}

	bool bContinued = false;
	if (m_AssignedPep.m_bPair && !nPeptideId)
		bContinued = true;

	if ((int) (pPeptide->m_tModCnt - pPeptide->m_tFixedModCnt)
			> m_Condition.m_nMaxModifyNumber)
	{
		return;
	}

	if (pPeptide->m_tLength <= 0)
	{
		return;
	}

	if (nCurrLen == -1)
	{


		_AttatchVarMod(nPeptideId, 0);

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

				_AttatchVarMod(nPeptideId, 0);

				pPeptide->m_lfMass -= m_mapVarMod[0][w].m_lfMass;
				--pPeptide->m_tModCnt;

			}
		}

	}
	else if (nCurrLen == (int) pPeptide->m_tLength)
	{
		if (bContinued)
		{
			_AttatchVarMod(1, -1);
		}
		else
		{


			size_t k;
			for (k = 0; k < m_Condition.m_vSelectedXLinker.size(); k++)
				_AttatchXLink(k);

			if ((int) (pPeptide->m_tModCnt - pPeptide->m_tFixedModCnt)
					< m_Condition.m_nMaxModifyNumber)
			{
				for (size_t i = 0; i < m_mapVarMod[1].size(); ++i)
				{
					if (m_mapVarMod[1][i].m_eModType == MT_PRO_CTERM
							&& !pPeptide->IsProCTerm())
						continue;


					if (m_Condition.m_vSelectedVarMod[m_mapVarMod[1][i].m_nIdx].m_strAA.find(
							pPeptide->m_szSequence[nCurrLen - 1])
							== string::npos)
						continue;

					pPeptide->m_tModSites[pPeptide->m_tModCnt][0] = nCurrLen
							- 1;
					pPeptide->m_tModSites[pPeptide->m_tModCnt][1] =
							m_mapVarMod[1][i].m_nIdx;

					pPeptide->m_lfMass += m_mapVarMod[1][i].m_lfMass;
					++pPeptide->m_tModCnt;

					for (k = 0; k < m_Condition.m_vSelectedXLinker.size(); k++)
						_AttatchXLink(k);

					pPeptide->m_lfMass -= m_mapVarMod[1][i].m_lfMass;
					--pPeptide->m_tModCnt;
				}
				return;
			}
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
			_AttatchVarMod(nPeptideId, nCurrLen);
		else
		{


			_AttatchVarMod(nPeptideId, nCurrLen + 1);

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

					_AttatchVarMod(nPeptideId, nCurrLen + 1);

					pPeptide->m_lfMass -= m_mapVarMod[cSite][i].m_lfMass;
					--pPeptide->m_tModCnt;
				}
				return;
			}
		}

	}
}

void CXLinkPepSalvoFlow::_printAssignedPep()
{
	char szBuf[1024];
	string strTXT = "";
	sprintf(szBuf, "%s ", m_AssignedPep.m_AlphaPeptide.m_szSequence);
	strTXT += szBuf;

	if (m_AssignedPep.m_bPair)
	{
		sprintf(szBuf, "-%s ", m_AssignedPep.m_BetaPeptide.m_szSequence);
		strTXT += szBuf;
	}

	if (m_AssignedPep.m_XLink.m_eXLinkType == 0)
	{
		sprintf(szBuf, "%d ", m_AssignedPep.m_XLink.m_eXLinkType);
		strTXT += szBuf;
	}
	else if (m_AssignedPep.m_XLink.m_eXLinkType == 1)
	{
		sprintf(szBuf, "%d %d ", m_AssignedPep.m_XLink.m_eXLinkType,
				m_AssignedPep.m_XLink.m_tAlphaSite);
		strTXT += szBuf;
	}
	else
	{
		sprintf(szBuf, "%d %d %d ", m_AssignedPep.m_XLink.m_eXLinkType,
				m_AssignedPep.m_XLink.m_tAlphaSite,
				m_AssignedPep.m_XLink.m_tBetaSite);
		strTXT += szBuf;
	}

	for (size_t j = 0; j < m_AssignedPep.m_AlphaPeptide.m_tModCnt; ++j)
	{
		sprintf(szBuf, "%d,%d ", m_AssignedPep.m_AlphaPeptide.m_tModSites[j][0],
				m_AssignedPep.m_AlphaPeptide.m_tModSites[j][1] + 1);
		strTXT += szBuf;
	}

	if (m_AssignedPep.m_bPair)
	{
		sprintf(szBuf, "|	");
		strTXT += szBuf;
		for (size_t j = 0; j < m_AssignedPep.m_BetaPeptide.m_tModCnt; ++j)
		{
			sprintf(szBuf, "%d,%d ",
					m_AssignedPep.m_BetaPeptide.m_tModSites[j][0],
					m_AssignedPep.m_BetaPeptide.m_tModSites[j][1] + 1);
			strTXT += szBuf;
		}
	}

	/*
	 sprintf(m_debug.m_szBuf,"PEP(PAIR):");
	 m_debug.WriteLine();

	 sprintf(m_debug.m_szBuf,"link type = %d",m_AssignedPep.m_XLink.m_eXLinkType);
	 m_debug.WriteLine();
	 sprintf(m_debug.m_szBuf,"site1 = %d",m_AssignedPep.m_XLink.m_tAlphaSite);
	 m_debug.WriteLine();
	 sprintf(m_debug.m_szBuf,"site2 = %d",m_AssignedPep.m_XLink.m_tBetaSite);
	 m_debug.WriteLine();
	 sprintf(m_debug.m_szBuf,"seq1 = %s",m_AssignedPep.m_AlphaPeptide.m_szSequence);
	 m_debug.WriteLine();
	 if(m_AssignedPep.m_bPair)
	 {
	 sprintf(m_debug.m_szBuf,"seq2 = %s",m_AssignedPep.m_BetaPeptide.m_szSequence);
	 m_debug.WriteLine();
	 }
	 sprintf(m_debug.m_szBuf,"mod on pep1 :");
	 m_debug.WriteLine();
	 for(size_t i = 0;i < m_AssignedPep.m_AlphaPeptide.m_tModCnt;++i)
	 {
	 sprintf(m_debug.m_szBuf,"%d:%d	",m_AssignedPep.m_AlphaPeptide.m_tModSites[i][0],m_AssignedPep.m_AlphaPeptide.m_tModSites[i][1]);
	 m_debug.WriteLine();
	 }

	 if(m_AssignedPep.m_bPair)
	 {
	 sprintf(m_debug.m_szBuf,"mod on pep2 :");
	 m_debug.WriteLine();
	 for(size_t i = 0;i < m_AssignedPep.m_BetaPeptide.m_tModCnt;++i)
	 {
	 sprintf(m_debug.m_szBuf,"%d:%d	",m_AssignedPep.m_BetaPeptide.m_tModSites[i][0],m_AssignedPep.m_BetaPeptide.m_tModSites[i][1]);
	 m_debug.WriteLine();
	 }
	 }
	 */

}

void CXLinkPepSalvoFlow::_track(size_t tIdx, size_t nLinkerId)
{



	CXLinkMatchResult & mr = m_vResults[tIdx];
	const CXLinkPepResult & pep_res = m_AssignedPep;
	++mr.m_tScore;

	const CSpectrum & spectrum = m_vSpectra[tIdx];
	double lfCalc_MH = 0.0;
	lfCalc_MH = mr.Calc_Theoretical_MH(pep_res, m_Condition.m_bPepMono);

	lfCalc_MH += _GetXLinkerMass(nLinkerId);

	double lfMZ = m_Condition.GetMZ(spectrum.m_lfMH, spectrum.m_nCharge);
	double lfPepTol = m_Condition.GetPepTol(lfMZ, spectrum.m_nCharge);
	double lfPepTolBase = m_Condition.GetPepTolBase(lfMZ, spectrum.m_nCharge);

	if (mr.m_vlfScores.size() < m_Condition.m_tMaxScoreNum)
		mr.m_vlfScores.push_back(pep_res.m_lfScore);

	if (fabs(spectrum.m_lfMH + lfPepTolBase - lfCalc_MH) <= lfPepTol)
		++mr.m_tRealCandidate;
	else
		return;

	if ((m_Condition.m_nReportPep <= 0) && (pep_res.m_lfScore <= 0))
		return;

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

void CXLinkPepSalvoFlow::_CreateLinkerMap()
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

bool CXLinkPepSalvoFlow::IsModLinkSiteConflict(const CSimplePeptide &peptide,
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
char CXLinkPepSalvoFlow::_GetCandLinkSiteIndex(const CSimplePeptide &peptide,
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

bool CXLinkPepSalvoFlow::IsLinkSiteCleave(const CSimplePeptide &peptide,
		const int nLinkSite)
{
#ifndef CLEAVESITE_LINKABLE
	if (nLinkSite >= peptide.m_tLength)
		cerr << "Error found in CXLinkPepSalvoFlow::IsLinkSiteCleave." << endl;

	string strCleaveString = m_Condition.m_SelectedEnzyme.GetCleaveString();
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
		if (nLinkSite != peptide.m_tLength - 1)
			return false;
		if (peptide.IsProCTerm())
			return false;
		if (strCleaveString.find(peptide.m_szSequence[nLinkSite]) != -1)
			return true;
	}
#endif

	return false;
}

void CXLinkPepSalvoFlow::_AttatchXLink(size_t nLinkerId)
{
	size_t i, j;

	char cSite1, cSite2;
	if (m_AssignedPep.m_bPair)
	{


		for (i = 0; i < m_AssignedPep.m_AlphaPeptide.m_tLength; ++i)
		{
			for (j = 0; j < m_AssignedPep.m_BetaPeptide.m_tLength; ++j)
			{


				cSite1 = _GetCandLinkSiteIndex(m_AssignedPep.m_AlphaPeptide, i,
						nLinkerId);
				cSite2 = _GetCandLinkSiteIndex(m_AssignedPep.m_BetaPeptide, j,
						nLinkerId);

				if (m_bLinked[nLinkerId][cSite1 - 'A'][cSite2 - 'A'])
				{
					if (IsModLinkSiteConflict(m_AssignedPep.m_AlphaPeptide, i)
							|| IsModLinkSiteConflict(
									m_AssignedPep.m_BetaPeptide, j))
						continue;
					if (IsLinkSiteCleave(m_AssignedPep.m_AlphaPeptide, i)
							|| IsLinkSiteCleave(m_AssignedPep.m_BetaPeptide, j))
						continue;
					m_AssignedPep.m_XLink.m_eXLinkType = X_LINK;
					m_AssignedPep.m_XLink.m_tAlphaSite = i;
					m_AssignedPep.m_XLink.m_tBetaSite = j;
					_Dispatch(nLinkerId);
				}
			}
		}
	}
	else
	{
#ifndef SUMO

		for (i = 0; i < m_AssignedPep.m_AlphaPeptide.m_tLength; ++i)
		{
			for (j = i + 1; j < m_AssignedPep.m_AlphaPeptide.m_tLength; ++j)
			{
				cSite1 = _GetCandLinkSiteIndex(m_AssignedPep.m_AlphaPeptide, i,
						nLinkerId);
				cSite2 = _GetCandLinkSiteIndex(m_AssignedPep.m_AlphaPeptide, j,
						nLinkerId);

				if (m_bLinked[nLinkerId][cSite1 - 'A'][cSite2 - 'A'])
				{
					if (IsModLinkSiteConflict(m_AssignedPep.m_AlphaPeptide, i)
							|| IsModLinkSiteConflict(
									m_AssignedPep.m_AlphaPeptide, j))
						continue;
					if (IsLinkSiteCleave(m_AssignedPep.m_AlphaPeptide, i)
							|| IsLinkSiteCleave(m_AssignedPep.m_AlphaPeptide,
									j))
						continue;
					m_AssignedPep.m_XLink.m_eXLinkType = LOOP_LINK;
					m_AssignedPep.m_XLink.m_tAlphaSite = i;
					m_AssignedPep.m_XLink.m_tBetaSite = j;
					_Dispatch(nLinkerId);
				}
			}
		}


		for (i = 0; i < m_AssignedPep.m_AlphaPeptide.m_tLength; ++i)
		{
			cSite1 = m_AssignedPep.m_AlphaPeptide.m_szSequence[i];

			cSite1 = _GetCandLinkSiteIndex(m_AssignedPep.m_AlphaPeptide, i,
					nLinkerId);


			if (m_bMonoLinked[nLinkerId][cSite1 - 'A'])
			{
				if (IsModLinkSiteConflict(m_AssignedPep.m_AlphaPeptide, i))
					continue;
				if (IsLinkSiteCleave(m_AssignedPep.m_AlphaPeptide, i))
					continue;
				m_AssignedPep.m_XLink.m_eXLinkType = MONO_LINK;
				m_AssignedPep.m_XLink.m_tAlphaSite = i;
				m_AssignedPep.m_XLink.m_tBetaSite = -1;
				_Dispatch(nLinkerId);
			}
		}
#endif

		m_AssignedPep.m_XLink.m_eXLinkType = NONE_LINK;
		m_AssignedPep.m_XLink.m_tAlphaSite = -1;
		m_AssignedPep.m_XLink.m_tBetaSite = -1;

		_Dispatch(nLinkerId);
	}

}

double CXLinkPepSalvoFlow::_GetXLinkerMass(size_t nLinkerId)
{
	switch (m_AssignedPep.m_XLink.m_eXLinkType)
	{
	case NONE_LINK:
		return 0.0;
	case MONO_LINK:
		if (m_Condition.m_bPepMono)
		{
			return m_Condition.m_vSelectedXLinker[nLinkerId].m_lfMLMonoMass_dif;
		}
		else
		{
			return m_Condition.m_vSelectedXLinker[nLinkerId].m_lfMLAvrgMass_dif;
		}

	case LOOP_LINK:

	case X_LINK:
		if (m_Condition.m_bPepMono)
		{
			return m_Condition.m_vSelectedXLinker[nLinkerId].m_lfMonoMass_dif;
		}
		else
		{
			return m_Condition.m_vSelectedXLinker[nLinkerId].m_lfAvrgMass_dif;
		}
	default:
		return 0.0;
	}
}

void CXLinkPepSalvoFlow::_Dispatch(size_t nLinkerId)
{
	double lfMass = 0.0;


	lfMass = CXLinkMatchResult::Calc_Theoretical_MH(m_AssignedPep,
			m_Condition.m_bPepMono);



	lfMass += _GetXLinkerMass(nLinkerId);

	m_AssignedPep.m_XLink.m_nLinkerId = nLinkerId;

	if (lfMass < m_vInfo[m_tSalvoLine].lfMin
			|| lfMass > m_vInfo[m_vInfo.size() - 1].lfMax)
	{
		return;
	}

	size_t tPos = 0;
	if ((int) lfMass >= m_SpecSearcher.m_nStart)
		tPos = m_SpecSearcher.m_pData[(int) lfMass - m_SpecSearcher.m_nStart];

#ifdef _DEBUG2

#endif

	m_pScorer->SetPeptide(m_AssignedPep);

	while (tPos < m_vInfo.size() && m_vInfo[tPos].lfMax < lfMass)
	{
		++tPos;
	}

	size_t tTemp = tPos;

	while (tTemp < m_vInfo.size())
	{
		int nIdx = m_vInfo[tTemp].nIndex;
		if (m_vInfo[tTemp].lfMin > lfMass)
		{
			break;
		}

		m_pScorer->SetSpectrum(m_vSpectra[nIdx]);

		m_AssignedPep.m_lfScore = _Score();









#ifdef	_DEBUG2

#endif






		_track((size_t) nIdx, nLinkerId);

		++tTemp;
	}
}

double CXLinkPepSalvoFlow::_Score()
{
	double lfScore = -1.0;
	try
	{
#ifdef _TEST_FLOW_RANDOM_SCORE
		lfScore = rand();
#else

#ifdef I2SCORE
		lfScore = m_pScorer->i2Score();
#else
		lfScore = m_pScorer->Score();
#endif

#endif
	} catch (exception & e)
	{
		CErrInfo info("CXLinkPepSalvoFlow", "_Score");
		char temp[200] =
		{ 0 };
		sprintf(temp, "lfScore=%.6lf", lfScore);
		info.Append(temp);
		throw runtime_error(info.Get(e).c_str());
	} catch (...)
	{
		CErrInfo info("CXLinkPepSalvoFlow", "_Score",
				"caught an unknown exception");
		char temp[200] =
		{ 0 };
		sprintf(temp, "lfScore=%.6lf", lfScore);
		info.Append(temp);
		throw runtime_error(info.Get().c_str());
	}
	return lfScore;

}

void CXLinkPepSalvoFlow::_CreateSpecIndex()
{
	m_vInfo.clear();
	SPEC_SORTER_INDEX_INFO stTemp;
	m_vSpecWndCnt.clear();
	for (size_t i = 0; i < m_vSpectra.size(); ++i)
	{
		double lfMHMin, lfMHMax;
		for (size_t j = 0; j < m_Condition.m_vPepTolWnds.size(); ++j)
		{
			stTemp.nIndex = i;
			stTemp.lfMass = m_vSpectra[i].m_lfMH;
			stTemp.nCharge = m_vSpectra[i].m_nCharge;
			m_vInfo.push_back(stTemp);
			m_Condition.GetMHMassBorder(m_vInfo[i].lfMass, m_vInfo[i].nCharge,
					m_vInfo[i].lfMin, m_vInfo[i].lfMax);
		}
		m_vSpecWndCnt.push_back(m_Condition.m_vPepTolWnds.size());
	}

	/*
	 * sort all spectra index items by their lower_bound mass
	 */

	sort(m_vInfo.begin(), m_vInfo.end(), MIN_MASS_LESS);

	m_SpecSearcher.Construct(m_vInfo);

}

void CXLinkPepSalvoFlow::_GetAllPeptide()
{


	_AttatchFixedMod(&m_AssignedPep.m_AlphaPeptide);
	if (m_AssignedPep.m_bPair)
		_AttatchFixedMod(&m_AssignedPep.m_BetaPeptide);


	_AttatchVarMod(0, -1);

}
void CXLinkPepSalvoFlow::_AttatchFixedMod(CSimplePeptide * const pPeptide)
{
	if (pPeptide->m_tLength <= 0)
		return;

	pPeptide->m_tFixedModCnt = 0;
	pPeptide->m_tModCnt = 0;
	size_t i = 0;


	if (m_mapFixMod.count(0))
	{
		for (i = 0; i < m_mapFixMod[0].size(); ++i)
		{
			if (m_mapFixMod[0][i].m_eModType == MT_PRO_NTERM
					&& !pPeptide->IsProNTerm())
				continue;


			size_t nIndex = m_mapFixMod[0][i].m_nIdx
					- m_Condition.m_vSelectedVarMod.size();

			if (nIndex < 0 || nIndex >= m_Condition.m_vSelectedFixMod.size())
				continue;

			if (m_Condition.m_vSelectedFixMod[nIndex].m_strAA.find(
					pPeptide->m_szSequence[0]) == string::npos)
				continue;

			if (pPeptide->m_tModCnt >= MAX_MODIFY_NUM)
				continue;

			pPeptide->m_tModSites[pPeptide->m_tModCnt][0] = 0;
			pPeptide->m_tModSites[pPeptide->m_tModCnt][1] =
					m_mapFixMod[0][i].m_nIdx;
			pPeptide->m_lfMass += m_mapFixMod[0][i].m_lfMass;
			++pPeptide->m_tModCnt;
		}
	}


	if (m_mapFixMod.count(1))
	{
		for (i = 0; i < m_mapFixMod[1].size(); ++i)
		{
			if (m_mapFixMod[1][i].m_eModType == MT_PRO_CTERM
					&& !pPeptide->IsProCTerm())
				continue;

			size_t nIndex = m_mapFixMod[1][i].m_nIdx
					- m_Condition.m_vSelectedVarMod.size();
			if (nIndex < 0 || nIndex >= m_Condition.m_vSelectedFixMod.size())
				continue;

			if (m_Condition.m_vSelectedFixMod[nIndex].m_strAA.find(
					pPeptide->m_szSequence[pPeptide->m_tLength - 1])
					== string::npos)
				continue;

			if (pPeptide->m_tModCnt >= MAX_MODIFY_NUM)
				continue;

			pPeptide->m_tModSites[pPeptide->m_tModCnt][0] = pPeptide->m_tLength
					- 1;
			pPeptide->m_tModSites[pPeptide->m_tModCnt][1] =
					m_mapFixMod[1][i].m_nIdx;
			pPeptide->m_lfMass += m_mapFixMod[1][i].m_lfMass;
			++pPeptide->m_tModCnt;
		}
	}


	for (i = 0; i < pPeptide->m_tLength; ++i)
	{
		char cSite = pPeptide->m_szSequence[i];
		if (m_mapFixMod.count(cSite))
		{
			for (size_t j = 0; j < m_mapFixMod[cSite].size(); ++j)
			{
				if (pPeptide->m_tModCnt >= MAX_MODIFY_NUM)
					continue;
				pPeptide->m_tModSites[pPeptide->m_tModCnt][0] = i;
				pPeptide->m_tModSites[pPeptide->m_tModCnt][1] =
						m_mapFixMod[cSite][j].m_nIdx;
				pPeptide->m_lfMass += m_mapFixMod[cSite][j].m_lfMass;
				++pPeptide->m_tModCnt;
			}
		}
	}
	pPeptide->m_tFixedModCnt = pPeptide->m_tModCnt;
}

void CXLinkPepSalvoFlow::_InitMatchResult()
{
	m_vResults.clear();
	m_vResults.reserve(m_vSpectra.size());
	CXLinkMatchResult mr;

	for (size_t i = 0; i < m_vSpectra.size(); ++i)
	{
		m_vResults.push_back(mr);
		m_vResults.back().m_vPeptideResults.reserve(m_Condition.m_nReportPep);
	}
}

void CXLinkPepSalvoFlow::_CloseDatabase()
{
	if (m_pReader == NULL)
		return;

	try
	{
		m_pReader->Close();
	} catch (exception & e)
	{
		CErrInfo info("CXLinkPepSalvoFlow", "_CloseDatabase",
				"in the function CPeptideReader::Close.");
		throw runtime_error(info.Get(e).c_str());
	} catch (...)
	{
		CErrInfo info("CXLinkPepSalvoFlow", "_CloseDatabase",
				"caught an unknown exception in the function CPeptideReader::Close.");
		throw runtime_error(info.Get().c_str());
	}

	delete m_pReader;
	m_pReader = NULL;

}

void CXLinkPepSalvoFlow::_OpenDatabase()
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
			CErrInfo info("CXLinkPepSalvoFlow", "_OpenDatabase",
					"strMeta==NULL");
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
		CErrInfo info("CXLinkPepSalvoFlow", "_OpenDatabase",
				"in the initialize of the database.");
		info.Append("strMeta=" + strMeta);
		info.Append("strDBPath=" + strDBPath);
		throw runtime_error(info.Get(e).c_str());
	} catch (...)
	{
		CErrInfo info("CXLinkPepSalvoFlow", "_OpenDatabase",
				"caught an unknown exception in the initialize of the database.");
		info.Append("strMeta=" + strMeta);
		info.Append("strDBPath=" + strDBPath);
		throw runtime_error(info.Get().c_str());
	}

}

void CXLinkPepSalvoFlow::_PrintPepPairList()
{
#ifdef _DEBUG2

	PEP_SQ stPepSQ,stPepSQ1;
#ifdef _USE_PEP_LIST
	CPepPairListReader listreader;
#else
	CPepPairContainerReader listreader;
#endif
	listreader.Init(&CXLinkPepSalvoFlow::m_CandidatePepPairs);
	listreader.Begin();
#ifdef _USE_PEP_LIST
	struct PEP_PAIR stPepPair;
#else
	struct PEP_PAIR_ITEM stPepPair;
#endif

	while(listreader.GetNext(stPepPair))
	{

		m_pReader->GetByID(stPepSQ,stPepPair.tPepId1);

		if(stPepPair.bPair)
		{
			m_pReader->GetByID(stPepSQ1,stPepPair.tPepId2);


			sprintf(m_debug.m_szBuf,"%f	%s %s",stPepPair.lfMass,stPepSQ.strSQ.c_str(),stPepSQ1.strSQ.c_str());
		}
		else
		{
			sprintf(m_debug.m_szBuf,"%f	%s",stPepPair.lfMass,stPepSQ.strSQ.c_str());
		}

		m_debug.WriteLine();
	}
#endif

}

void CXLinkPepSalvoFlow::CreatePepPairList()
{
	_OpenDatabase();
	_CreatePepPairList();
	_CloseDatabase();
}

void CXLinkPepSalvoFlow::_CreatePepPairList()
{
#ifdef _USE_PEP_LIST	

	if(CXLinkPepSalvoFlow::m_bPepPairListCreated)
	return;
	CXLinkPepSalvoFlow::m_bPepPairListCreated = true;
	char tmpch;

	size_t tCurId;
	tCurId = 0;

	PEP_SQ stPepSQ,stPepSQ1;

	CPepPairList curPepList;
	CPepPairList oldPepList;
	tCurId = 0;

	oldPepList.Init();
	while(m_pReader->GetNext(stPepSQ))
	{
		oldPepList.Attach(false,tCurId,0,stPepSQ.dfMass);
		++tCurId;
	}

	CPepPairListReader listreader1,listreader2;

	listreader1.Init(&oldPepList);
	listreader2.Init(&oldPepList);
	listreader1.Begin();
	listreader2.Begin();

	struct PEP_PAIR stPepPair1 , stPepPair2;
	tCurId = 0;

	while(listreader1.GetNext(stPepPair1))
	{
		curPepList.Init();
		curPepList.Attach(false,stPepPair1.tPepId1,0,stPepPair1.lfMass);
		listreader2.Begin();
		for(size_t i=0;i<=tCurId;++i)
		{
			listreader2.GetNext(stPepPair2);
			curPepList.Attach(true,stPepPair1.tPepId1,stPepPair2.tPepId1,stPepPair1.lfMass + stPepPair2.lfMass);
		}
		++tCurId;



		CXLinkPepSalvoFlow::m_CandidatePepPairs.Merge(&curPepList);


	}




#else



	if (CXLinkPepSalvoFlow::m_bPepPairListCreated)
		return;
	CXLinkPepSalvoFlow::m_bPepPairListCreated = true;

	CXLinkPepSalvoFlow::m_CandidatePepPairs.Init();

	size_t tCurId;
	tCurId = 0;

	PEP_SQ stPepSQ, stPepSQ1;

	CPepPairContainer curPepList;
	CPepPairContainer oldPepList;

#ifdef SUMO	CPepPairContainer oldPepList2;

	oldPepList2.Init();
#endif

	oldPepList.Init();

#ifdef PEPMASSNOCOUNT
	int iTotalPep = 0;
	int SingleMassDistribution[20001] =
	{	0};
	vector<double> singleMass;
#endif

	while (m_pReader->GetNext(stPepSQ))
	{

#ifdef SUMO
		if (stPepSQ.strSQ.length() >= PEPMINLEN)
		{
			oldPepList.Attach(false, tCurId, 0, stPepSQ.dfMass);
		}

#else
		oldPepList.Attach(false, tCurId, 0, stPepSQ.dfMass);
#endif

#ifdef PEPMASSNOCOUNT
		singleMass.push_back(stPepSQ.dfMass);
#endif
#ifdef SUMO





		if (stPepSQ.strSQ == SUMOSEQ)
		{






			oldPepList2.Attach(false,tCurId,0,stPepSQ.dfMass);
		}
#endif

		++tCurId;

#ifdef PEPMASSNOCOUNT
		++iTotalPep;
		if ( (int)(stPepSQ.dfMass) > 20000 )
		SingleMassDistribution[20000] += 1;
		else
		SingleMassDistribution[(int)(stPepSQ.dfMass)] += 1;
#endif
	}

#ifdef PEPMASSNOCOUNT
	fstream fout;
	fout.open("E:\\SinglePep.txt", ios::out | ios::app);
	fout<<"Single pep total: "<<iTotalPep<<endl;
	for (int i = 0; i < 20001; i++)
	fout<<i<<"\t"<<SingleMassDistribution[i]<<endl;

	fout.close();
#endif

	oldPepList.EndAttach();
#ifdef SUMO
	oldPepList2.EndAttach();
#endif
	CPepPairContainerReader listreader1, listreader2;








	listreader1.Init(&oldPepList);
#ifdef SUMO
	listreader2.Init(&oldPepList2);
#else
	listreader2.Init(&oldPepList);
#endif

	listreader1.Begin();
	listreader2.Begin();

	struct PEP_PAIR_ITEM stPepPair1, stPepPair2;
	tCurId = 0;

#ifdef PEPMASSNOCOUNT
	int iTotalCrossPep = 0;
	int CrossMassDistribution[20001] =
	{	0};
#endif

#ifdef SUMO
	bool bHaveCommonAdded = 0;
#endif
	while (listreader1.GetNext(stPepPair1))
	{
		CXLinkPepSalvoFlow::m_CandidatePepPairs.Attach(false,
				stPepPair1.tPepId1, 0, stPepPair1.lfMass);

		listreader2.Begin();
#ifdef SUMO
		while (listreader2.GetNext(stPepPair2))
		{
			if (!bHaveCommonAdded)
			CXLinkPepSalvoFlow::m_CandidatePepPairs.Attach(false,
					stPepPair2.tPepId1, 0, stPepPair2.lfMass);

			CXLinkPepSalvoFlow::m_CandidatePepPairs.Attach(true,
					stPepPair1.tPepId1, stPepPair2.tPepId1,
					stPepPair1.lfMass + stPepPair2.lfMass);
		}
		bHaveCommonAdded = 1;
#else
		for (size_t i = 0; i <= tCurId; ++i)
		{
			listreader2.GetNext(stPepPair2);

			CXLinkPepSalvoFlow::m_CandidatePepPairs.Attach(true,
					stPepPair1.tPepId1, stPepPair2.tPepId1,
					stPepPair1.lfMass + stPepPair2.lfMass);

#ifdef PEPMASSNOCOUNT
			++iTotalCrossPep;
			if ( (int)(stPepPair1.lfMass + stPepPair2.lfMass) > 20000 )
			CrossMassDistribution[20000] += 1;
			else
			CrossMassDistribution[(int)(stPepPair1.lfMass + stPepPair2.lfMass)] += 1;
#endif
		}
		++tCurId;
#endif
	}

#ifdef PEPMASSNOCOUNT
	fstream fout2;
	fout2.open("CrossPep.txt", std::fstream::in | std::fstream::out | std::fstream::app);
	fout2<<"Cross pep total: "<<iTotalCrossPep<<endl;
	for (int i = 0; i < 20001; i++)
	fout2<<i<<"\t"<<CrossMassDistribution[i]<<endl;

	fout2.close();
#endif

	CXLinkPepSalvoFlow::m_CandidatePepPairs.EndAttach();
	CXLinkPepSalvoFlow::m_CandidatePepPairs.Sort();



#endif
}

void CXLinkPepSalvoFlow::Run(const string & strTmpFile)
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
		CErrInfo info("CXLinkPepSalvoFlow", "Run",
				"in the function CXLinkPepSalvoFlow::_PreProc");
		throw runtime_error(info.Get(e).c_str());
	} catch (...)
	{
		CErrInfo info("CXLinkPepSalvoFlow", "Run",
				"caught an unknown exception");
		throw runtime_error(info.Get().c_str());
	}

	m_pTrace->Info("openning database...");
	m_pState->GetProgress(m_strInfo);
	m_pTrace->Info(m_strInfo);
	_OpenDatabase();

	m_pTrace->Info("creating the peptide pair list...");
	m_pState->GetProgress(m_strInfo);
	m_pTrace->Info(m_strInfo);
	_CreatePepPairList();

	m_pTrace->Info("Reading the peptide sequences...");
	m_pState->GetProgress(m_strInfo);
	m_pTrace->Info(m_strInfo);

	/*pep list*/
#ifdef _USE_PEP_LIST
	CPepPairListReader pepPairListReader;
#else
	CPepPairContainerReader pepPairListReader;
#endif

	pepPairListReader.Init(&CXLinkPepSalvoFlow::m_CandidatePepPairs);
	pepPairListReader.Begin();

	size_t nCurrent(0);
	m_tSalvoLine = 0;
	bool bPepEnd(false);
	/*pep list*/
#ifdef _USE_PEP_LIST
	struct PEP_PAIR stPepPair;
#else
	struct PEP_PAIR_ITEM stPepPair;
#endif
	double lfPepPairMass(0);
	double lfLeastNegativeDif = 0;
	lfLeastNegativeDif = m_Condition.GetLeastNegativeMod()
			* m_Condition.m_nMaxModifyNumber;

	while (1)
	{
		try
		{
			bPepEnd = !pepPairListReader.GetNext(stPepPair);
			if (bPepEnd && (m_tSalvoLine >= m_vInfo.size()))
				break;

			lfPepPairMass = stPepPair.lfMass;
		} catch (exception & e)
		{
			CErrInfo info("CXLinkPepSalvoFlow", "Run",
					"in the function pepPairListReader::GetNext.");
			throw runtime_error(info.Get(e).c_str());
		} catch (...)
		{
			CErrInfo info("CXLinkPepSalvoFlow", "Run",
					"caught an unknown exception in the function pepPairListReader::GetNext.");
			throw runtime_error(info.Get().c_str());
		}

		while (m_tSalvoLine < m_vInfo.size())
		{
			double lfMinMass;
			lfMinMass = lfPepPairMass + lfLeastNegativeDif;
			if (stPepPair.bPair)
				lfMinMass += lfLeastNegativeDif;










			if ((lfMinMass > m_vInfo[m_tSalvoLine].lfMax) || bPepEnd)
			{
				int nIdx = 0;
				nIdx = m_vInfo[m_tSalvoLine].nIndex;
				try
				{
					_Evaluate(nIdx);

					vector<double>().swap(m_vResults[nIdx].m_vlfScores);
				} catch (exception & e)
				{
					CErrInfo info("CXLinkPepSalvoFlow", "Run",
							"in the function CXLinkPepSalvoFlow::Evaluate.");
					char temp[200];
					sprintf(temp, "current=%d", nIdx);
					info.Append(temp);
					throw runtime_error(info.Get(e).c_str());
				} catch (...)
				{
					CErrInfo info("CXLinkPepSalvoFlow", "Run",
							"caught an unknown exception in the function CXLinkPepSalvoFlow::Evaluate.");
					char temp[200];
					sprintf(temp, "current=%d", nIdx);
					info.Append(temp);
					throw runtime_error(info.Get().c_str());
				}
				++m_tSalvoLine;
			}
			else
				break;
		}

		if (m_tSalvoLine >= m_vInfo.size())
			break;

		if (nCurrent % 10000 == 0)
		{
			m_pState->SetThreadCurrent(m_tID, m_tSalvoLine);
			m_pState->GetProgress(m_strInfo);
			m_pTrace->Info(m_strInfo);
		}

		++nCurrent;

		try
		{
			PEP_SQ stPepSQ;

			m_pReader->GetByID(stPepSQ, stPepPair.tPepId1);
			m_AssignedPep.m_AlphaPeptide.SetPeptideInfor(stPepSQ.strSQ.c_str(),
					stPepSQ.strSQ.length(), stPepSQ.dfMass, stPepSQ.cMiss,
					stPepSQ.cEnd, false);
			m_AssignedPep.m_AlphaPeptide.m_tFixedModCnt = 0;
			m_AssignedPep.m_AlphaPeptide.m_tModCnt = 0;

			if (stPepPair.bPair)
			{
				m_AssignedPep.m_bPair = true;
				m_pReader->GetByID(stPepSQ, stPepPair.tPepId2);
				m_AssignedPep.m_BetaPeptide.SetPeptideInfor(
						stPepSQ.strSQ.c_str(), stPepSQ.strSQ.length(),
						stPepSQ.dfMass, stPepSQ.cMiss, stPepSQ.cEnd, false);
				m_AssignedPep.m_BetaPeptide.m_tFixedModCnt = 0;
				m_AssignedPep.m_BetaPeptide.m_tModCnt = 0;
			}
			else
			{
				m_AssignedPep.m_bPair = false;
				m_AssignedPep.m_BetaPeptide.SetPeptideInfor("", 0, 0.0, 0, 0,
						false);
				m_AssignedPep.m_BetaPeptide.m_tFixedModCnt = 0;
				m_AssignedPep.m_BetaPeptide.m_tModCnt = 0;
			}
#ifdef	_DEBUG2

			if (stPepPair.bPair)
			sprintf(m_debug.m_szBuf, "%f %s %s", stPepPair.lfMass,
					m_AssignedPep.m_AlphaPeptide.m_szSequence,
					m_AssignedPep.m_BetaPeptide.m_szSequence);

			else
			sprintf(m_debug.m_szBuf, "%f %s", stPepPair.lfMass,
					m_AssignedPep.m_AlphaPeptide.m_szSequence);

			m_debug.WriteLine();

#endif
		} catch (exception & e)
		{
			CErrInfo info("CXLinkPepSalvoFlow", "Run",
					"in the function CMass2PepIndexReader::GetByID.");
			throw runtime_error(info.Get(e).c_str());
		} catch (...)
		{
			CErrInfo info("CXLinkPepSalvoFlow", "Run",
					"caught an unknown exception in the function CMass2PepIndexReader::GetByID.");
			throw runtime_error(info.Get().c_str());
		}




		try
		{
#ifdef _TEST_FLOW_WITHOUT_DISPATCH

#else

			if (m_AssignedPep.m_bPair)
			{
				if (_IsXLinked(m_AssignedPep.m_AlphaPeptide)
						&& _IsXLinked(m_AssignedPep.m_BetaPeptide))
				{
					_GetAllPeptide();
				}
			}
			else
			{
				_GetAllPeptide();
			}

#endif        	
		} catch (exception & e)
		{
			CErrInfo info("CXLinkPepSalvoFlow", "Run",
					"in the function _GetAllPeptide().");
		} catch (...)
		{
			CErrInfo info("CXLinkPepSalvoFlow", "Run",
					"caught an unknown exception in the function _GetAllPeptide().");
		}

	}

	m_pTrace->Info("Dispatching and matching complete ...");

	try
	{

		_WritePFDFile(strTmpFile);
	} catch (exception & e)
	{
		CErrInfo info("CXLinkPepSalvoFlow", "Run",
				"in the function CXLinkPepSalvoFlow::_WritePFDFile.");
		info.Append("strTmpFile=" + strTmpFile);
		throw runtime_error(info.Get(e).c_str());
	} catch (...)
	{
		CErrInfo info("CXLinkPepSalvoFlow", "Run",
				"caught an unknown exception in the function CXLinkPepSalvoFlow::_WritePFDFile.");
		info.Append("strTmpFile=" + strTmpFile);
		throw runtime_error(info.Get().c_str());
	}

	_CloseDatabase();

	vector<CXLinkMatchResult>().swap(m_vResults);
	vector<CSpectrum>().swap(m_vSpectra);
	vector<SPEC_SORTER_INDEX_INFO>().swap(m_vInfo);
}

bool CXLinkPepSalvoFlow::_IsXLinked(const CSimplePeptide &peptide)
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

void CXLinkPepSalvoFlow::_ResultCheck()
{
	for (size_t tIdx = 0; tIdx < m_vResults.size(); ++tIdx)
	{
		for (size_t tOrder = 0;
				tOrder < m_vResults[tIdx].m_vPeptideResults.size(); ++tOrder)
		{
			if (m_vResults[tIdx].m_vPeptideResults[tOrder].m_XLink.m_eXLinkType
					!= X_LINK)
				continue;
			CXLinkPepResult &pep = m_vResults[tIdx].m_vPeptideResults[tOrder];
			if (pep.m_AlphaPeptide.m_lfMass < pep.m_BetaPeptide.m_lfMass)
			{
				pep.SwapAlphaBeta();
			}
		}
	}
}

void CXLinkPepSalvoFlow::_WritePFDFile(const string & strOutputPath)
{
	try
	{
		CXLinkPFDFileIO pfdIO;
		_ResultCheck();
		pfdIO.WriteFile(strOutputPath, m_vResults, m_vSpectra);
	} catch (exception & e)
	{
		CErrInfo info("CXLinkPepSalvoFlow", "_WritePFDFile",
				"in the function CPFDIO::WriteFile.");
		info.Append("strOutputPath=" + strOutputPath);
		throw runtime_error(info.Get(e).c_str());
	} catch (...)
	{
		CErrInfo info("CXLinkPepSalvoFlow", "_WritePFDFile",
				"caught an unknown exception in the function CPFDIO::WriteFile.");
		info.Append("strOutputPath=" + strOutputPath);
		throw runtime_error(info.Get().c_str());
	}

}

void CXLinkPepSalvoFlow::_Preproc()
{
	for (size_t i = 0; i < m_vSpectra.size(); ++i)
	{
		if (0 == m_vSpectra[i].m_tPeaksNum)
			continue;

		CSpectrum Output = m_vSpectra[i];
		try
		{
			m_pPreProc->Run(m_vSpectra[i], Output);
		} catch (exception & e)
		{
			CErrInfo info("CXLinkPepSalvoFlow", "_Preproc",
					"in the function CPreProcess::Run.");
			throw runtime_error(info.Get(e).c_str());
		} catch (...)
		{
			CErrInfo info("CXLinkPepSalvoFlow", "_Preproc",
					"caught an unknown exception in the function CPreProcess::Run.");
			throw runtime_error(info.Get().c_str());
		}

		m_vSpectra[i] = Output;
	}
}

void CXLinkPepSalvoFlow::_GenerateRedundantPepCandidates(
		const CXLinkPepResult & modelPepCandidate,
		vector<CXLinkPepResult> & vPepCandidates, int nMinNum)
{
	vPepCandidates.clear();

	if (nMinNum <= 0)
		return;


	double lfReductMass = 0.0;

	char cSite1 = 0, cSite2 = 0;
	double lfMass1 = 0.0;
	double lfMass2 = 0.0;


	/*
	 if(modelPepCandidate.m_XLink.m_eXLinkType == MONO_LINK)
	 {
	 if(m_Condition.m_bPepMono)
	 {
	 lfReductMass = m_Condition.m_SelectedXLinker.m_lfMLMonoMass_dif;
	 }
	 else
	 {
	 lfReductMass = m_Condition.m_SelectedXLinker.m_lfMLAvrgMass_dif;
	 }
	 }
	 else if(modelPepCandidate.m_XLink.m_eXLinkType == LOOP_LINK || modelPepCandidate.m_XLink.m_eXLinkType == X_LINK)
	 {
	 if(m_Condition.m_bPepMono)
	 {
	 lfReductMass = m_Condition.m_SelectedXLinker.m_lfMonoMass_dif;
	 }
	 else
	 {
	 lfReductMass = m_Condition.m_SelectedXLinker.m_lfAvrgMass_dif;
	 }
	 }
	 */


#ifndef SUMO
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
#endif

	double lfMass;
	lfMass = CXLinkMatchResult::Calc_Theoretical_MH(modelPepCandidate, false)
			- lfReductMass;

#ifndef SUMO
	if (modelPepCandidate.m_bPair)
	{
		lfMass -= 2 * (IonMass_Aver_O + IonMass_Aver_H * 2.0) + IonMass_Aver_H;
	}
	else
	{
		lfMass -= IonMass_Aver_O + IonMass_Aver_H * 3.0;
	}
#else
	lfMass -= IonMass_Aver_O + IonMass_Aver_H * 3.0;
#endif

	vector<string> vAASeq;
	vector<string> vAASeq1, vAASeq2;
	int nNum1 = 0, nNum2 = 0, nNum = 0;
	nNum1 = nMinNum / 3;
	nNum2 = nMinNum / 3;
	nNum = nMinNum - nNum1 - nNum2;
	if (nNum < 0)
		nNum = 0;

#ifndef SUMO
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
#else
	m_pGenerator->Batch(lfMass, nMinNum, vAASeq);
#endif



	CXLinkPepResult pepCandidate;
#ifndef SUMO
	pepCandidate.m_bPair = modelPepCandidate.m_bPair;
	pepCandidate.m_XLink = modelPepCandidate.m_XLink;
#else
	pepCandidate.m_bPair = 0;
	pepCandidate.m_XLink.m_eXLinkType = NONE_LINK;
#endif

	int nLen1, nLen2;
	int nRandomPos;
	string str1, str2, str3;
	char szBuf[MAX_PEPTIDE_LENGTH + 1];

	if (pepCandidate.m_bPair)
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
		if (pepCandidate.m_bPair)
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
			if (pepCandidate.m_XLink.m_eXLinkType == LOOP_LINK)
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
			else if (pepCandidate.m_XLink.m_eXLinkType == MONO_LINK)
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

void CXLinkPepSalvoFlow::_Evaluate(size_t idx)
{
	CSpectrum & spec = m_vSpectra[idx];
	CXLinkMatchResult & result = m_vResults[idx];

	if (result.m_vPeptideResults.empty())
		return;


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
				m_pScorer->SetSpectrum(m_vSpectra[idx]);
				m_pScorer->SetPeptide(vPepCandidates[j]);
#ifdef _DEBUG2	
				char szBuf[1024];
				string strTXT = "";
				sprintf(szBuf,"%s ",vPepCandidates[j].m_AlphaPeptide.m_szSequence);
				strTXT += szBuf;

				if(vPepCandidates[j].m_bPair)
				{
					sprintf(szBuf,"-%s ",vPepCandidates[j].m_BetaPeptide.m_szSequence);
					strTXT += szBuf;
				}

				if(vPepCandidates[j].m_XLink.m_eXLinkType == 0)
				{
					sprintf(szBuf,"%d ",vPepCandidates[j].m_XLink.m_eXLinkType);
					strTXT += szBuf;
				}
				else if(vPepCandidates[j].m_XLink.m_eXLinkType == 1)
				{
					sprintf(szBuf,"%d %d ",vPepCandidates[j].m_XLink.m_eXLinkType,vPepCandidates[j].m_XLink.m_tAlphaSite);
					strTXT += szBuf;
				}
				else
				{
					sprintf(szBuf,"%d %d %d ",vPepCandidates[j].m_XLink.m_eXLinkType,vPepCandidates[j].m_XLink.m_tAlphaSite,vPepCandidates[j].m_XLink.m_tBetaSite);
					strTXT += szBuf;
				}

				for(size_t j=0;j<vPepCandidates[j].m_AlphaPeptide.m_tModCnt;++j)
				{
					sprintf(szBuf,"%d,%d ",vPepCandidates[j].m_AlphaPeptide.m_tModSites[j][0],vPepCandidates[j].m_AlphaPeptide.m_tModSites[j][1]+1);
					strTXT += szBuf;
				}

				if(vPepCandidates[j].m_bPair)
				{
					sprintf(szBuf,"|	");
					strTXT += szBuf;
					for(size_t j=0;j<vPepCandidates[j].m_BetaPeptide.m_tModCnt;++j)
					{
						sprintf(szBuf,"%d,%d ",vPepCandidates[j].m_BetaPeptide.m_tModSites[j][0],vPepCandidates[j].m_BetaPeptide.m_tModSites[j][1]+1);
						strTXT += szBuf;
					}
				}

				sprintf(m_debug.m_szBuf,"%s",strTXT.c_str());
				m_debug.WriteLine();
#endif

				double lfScore = _Score();
				++result.m_tScore;
				result.m_vlfScores.push_back(lfScore);
			}
		}
	}









#ifdef I2SCORE
	for (int i = 0; i <result.m_vPeptideResults.size(); i++)
	result.m_vPeptideResults[i].m_lfEvalue = 1 - result.m_vPeptideResults[i].m_lfScore;
	return;
#endif
	try
	{
		if (!result.m_vlfScores.empty())
		{









			m_pEV->Run(spec, &result);


		}
#ifdef _DEBUG2
		/*
		 cout << m_vResults[idx].m_vPeptideResults.size() << endl;
		 for(size_t t = 0;t < m_vResults[idx].m_vPeptideResults.size();++t)
		 {
		 cout << m_vResults[idx].m_vPeptideResults[t].m_peptide.m_szSequence << ' ' << m_vResults[idx].m_vPeptideResults[t].m_lfEvalue << ' '
		 << m_vResults[idx].m_vPeptideResults[t].m_peptide.m_lfMass << ' '
		 << m_vResults[idx].m_vPeptideResults[t].m_bEV << endl;
		 }
		 */
#endif
	} catch (exception & e)
	{
		CErrInfo info("CXLinkPepSalvoFlow", "_Evaluate",
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
		CErrInfo info("CXLinkPepSalvoFlow", "_Evaluate",
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

void CXLinkPepSalvoFlow::Close(void)
{
	m_Condition.clear();

	if (m_pPreProc)
		delete m_pPreProc;
	m_pPreProc = NULL;

	if (m_pScorer)
		delete m_pScorer;
	m_pScorer = NULL;

	if (m_pEV)
		delete m_pEV;
	m_pEV = NULL;

	if (m_pGenerator)
		delete m_pGenerator;
	m_pGenerator = NULL;
#ifdef _DEBUG2
	m_debug.Close();
#endif

}

