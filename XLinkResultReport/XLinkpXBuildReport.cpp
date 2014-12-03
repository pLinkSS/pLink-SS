#include <iostream>
#include <string>
 
#include "../include/predefine.h"
#include "../include/sdk.h"
#include "../include/interface.h"
#include "../include/option.h"

#include "XLinkResultReport.h"
#include "XLinkpXBuildReport.h"

using namespace std;
using namespace proteomics_sdk;

CXLinkpXBuildReport::CXLinkpXBuildReport()
{
	char szBuf[1024] = {0};
	getcwd(szBuf, 1024);
    m_strWorkDir = szBuf;
    if(m_strWorkDir[m_strWorkDir.length() - 1] != SLASH)
    {
    	m_strWorkDir += SLASH;
    }
}

CXLinkpXBuildReport::~CXLinkpXBuildReport()
{
}

void CXLinkpXBuildReport::Init(string strpFindFile,time_t tmStartTime)
{
	m_strpFindFile = strpFindFile;
	m_tmStartTime = tmStartTime;
	
	CConditionReader reader(m_strpFindFile, m_strWorkDir);
	try
	{
		reader.Read();
	}
	catch(exception & e)
	{
		CErrInfo info("CXLinkpXBuildReport", "Init", "in the function CConditionReader::Read.");
		info.Append("strOption=" + strpFindFile);
		info.Append("m_strWorkDir=" + m_strWorkDir);
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CXLinkpXBuildReport", "Init", "caught an unknown exception in the function CConditionReader::Read.");
		info.Append("strOption=" + strpFindFile);
		info.Append("m_strWorkDir=" + m_strWorkDir);
		throw runtime_error(info.Get().c_str());
	}
	
	m_Condition = reader.m_Condition;
		
    for(size_t i = 0;i < m_Condition.m_vSelectedFixMod.size();++i)
    	m_Condition.m_vSelectedVarMod.push_back(m_Condition.m_vSelectedFixMod[i]);
    
}

double CXLinkpXBuildReport::_GetDeltaMass(double lfCalcMH, double lfExpMH)
{
	double lfDeltaMass = lfExpMH - lfCalcMH;

	/*
	if(int(lfDeltaMass + 0.5) == 1)
	{
		lfDeltaMass -= 1.007825035;
	}
	else if(int(lfDeltaMass + 0.5) == 2)
	{
		lfDeltaMass -= 2.01565007;
	}
	else if(int(lfDeltaMass + 0.5) == 3)
	{
		lfDeltaMass -= 3.023475105;
	}
	else if(int(lfDeltaMass + 0.5) == 4)
	{
		lfDeltaMass -= 4.03130014;
	}
	*/
	
	/*
	if ( 0 == (m_Condition.m_vPepTolWnds[0].m_strPepTolType.compare("%"))) 
		return lfDeltaMass / lfExpMH * 100;
	else if ( 0 == (m_Condition.m_vPepTolWnds[0].m_strPepTolType.compare("mmu")))
		return lfDeltaMass / lfExpMH * 1000;
	else if ( 0 == (m_Condition.m_vPepTolWnds[0].m_strPepTolType.compare("ppm")))
		return lfDeltaMass / lfExpMH * 1000000;
	else
		return lfDeltaMass;
	*/
	
	return lfDeltaMass ;
}

double CXLinkpXBuildReport::_GetXLinkerMass(const CXLinkPepResult & pep_res)
{
	switch(pep_res.m_XLink.m_eXLinkType)
	{
	case NONE_LINK:
		return 0.0;
	case MONO_LINK:
		if(m_Condition.m_bPepMono)
		{
			return m_Condition.m_vSelectedXLinker[pep_res.m_XLink.m_nLinkerId].m_lfMLMonoMass_dif;
		}
		else
		{
			return m_Condition.m_vSelectedXLinker[pep_res.m_XLink.m_nLinkerId].m_lfMLAvrgMass_dif;
		}
		
	case LOOP_LINK:
		
	case X_LINK:
		if(m_Condition.m_bPepMono)
		{
			return m_Condition.m_vSelectedXLinker[pep_res.m_XLink.m_nLinkerId].m_lfMonoMass_dif;
		}
		else
		{
			return m_Condition.m_vSelectedXLinker[pep_res.m_XLink.m_nLinkerId].m_lfAvrgMass_dif;
		}
	default:
		return 0.0;
	}
}

void CXLinkpXBuildReport::GetLines(const vector<CXLinkMatchResult> & vResults , const vector<CSpectrum> & vSpectra ,string & strTXT)
{
	char szBuf[1024] = {'\0'};
	
	strTXT.clear();
	size_t tCnt = 0;
	for(size_t i=0; i<vResults.size(); ++i)
	{
		if(vResults[i].m_vPeptideResults.empty())
			continue;	

		strTXT += vSpectra[i].m_strFilePath;
		strTXT += "	";
		
		const CXLinkPepResult & pep_res = vResults[i].m_vPeptideResults[0];
		
		sprintf(szBuf,"%E	%f	",pep_res.m_lfEvalue,pep_res.m_lfScore);
		strTXT += szBuf;
		
		
		// calculate continous amino num (the minimum of pep1 and pep2)
		int nMaxLen = 0 , nCurLen = 0;
		for(size_t k=0;k<pep_res.m_AlphaPeptide.m_tLength;++k)
		{
			if(pep_res.m_stMatchInfo.aPepConf[k] > 0)
			{
				nCurLen ++;
			}
			else
			{
				if(nCurLen > nMaxLen)
					nMaxLen = nCurLen;
				nCurLen = 0;
			}	
		}
		if(nCurLen > nMaxLen)
			nMaxLen = nCurLen;
		
		if(pep_res.m_bPair)
		{
			int nMaxLen1 = 0;
			nCurLen = 0;
			for(size_t k=0;k<pep_res.m_BetaPeptide.m_tLength;++k)
			{
				if(pep_res.m_stMatchInfo.aPepConf[pep_res.m_AlphaPeptide.m_tLength + k] > 0)
				{
					nCurLen ++;
				}
				else
				{
					if(nCurLen > nMaxLen1)
						nMaxLen1 = nCurLen;
					nCurLen = 0;
				}	
			}
			if(nCurLen > nMaxLen1)
				nMaxLen1 = nCurLen;

			if(nMaxLen1 < nMaxLen)
				nMaxLen = nMaxLen1;
		}
			
		sprintf(szBuf,"%d	%f	",nMaxLen,pep_res.m_stMatchInfo.lfMatchedSpecInt/(pep_res.m_stMatchInfo.lfUnMatchedSpecInt+pep_res.m_stMatchInfo.lfMatchedSpecInt));
		strTXT += szBuf;
		
		
		// calcuate delta mass
		double lfMass1 = 0,lfMass2 = 0;
		if(pep_res.m_lfCalc_MH < 0.0001)
		{
			lfMass1 = CXLinkMatchResult::Calc_Theoretical_MH(pep_res, m_Condition.m_bPepMono);
			lfMass1 += _GetXLinkerMass(pep_res);
		}
		else
			lfMass1 = pep_res.m_lfCalc_MH;
		
		lfMass2 = _GetDeltaMass(lfMass1,vSpectra[i].m_lfMH);
		sprintf(szBuf,"%f	",lfMass2);
		strTXT += szBuf;
		
		// get charge 
		sprintf(szBuf,"%d	",vSpectra[i].m_nCharge);
		strTXT += szBuf;
		// see the tag : T U F
		char cTag = 'T';
		
		bool bForward1 = false,bForward2 = false;
		for(size_t k=0 ; k < pep_res.m_vAlphaProteinAC.size(); k++)
		{
			
			if(m_Condition.m_strFPSign != pep_res.m_vAlphaProteinAC[k].substr(0, m_Condition.m_strFPSign.size()))
			{
				bForward1 = true;
				break;
			}
		}
		
		if(pep_res.m_bPair)
		{
			bForward2 = false;
			for(size_t k=0 ; k < pep_res.m_vBetaProteinAC.size(); k++)
			{
				
				if(m_Condition.m_strFPSign != pep_res.m_vBetaProteinAC[k].substr(0, m_Condition.m_strFPSign.size()))
				{
					bForward2 = true;
					break;
				}
			}
			
			if(bForward1 && bForward2)
			{
				// true
				cTag = 'T';
			}
			else if(bForward1 || bForward2)
			{
				// undefined
				cTag = 'U';
			}
			else
			{
				// false
				cTag = 'F';
			}
			
		}
		else
		{
			if(bForward1)
			{
				// true
				cTag = 'T';
			}
			else
			{
				// false
				cTag = 'F';
			}
			
		}
		
		sprintf(szBuf,"%c	\n",cTag);
		strTXT += szBuf;
	}

}

void CXLinkpXBuildReport::WriteFile(const vector<CXLinkMatchResult> & vResults , const vector<CSpectrum> & vSpectra , string strOutputPath , string strTitle)
{
	string strTXTContent;

	size_t t = m_strpFindFile.find(".pfind");
	if(string::npos == t)
	{
		t = m_strpFindFile.find(".PFIND");
		if(string::npos == t)
		{
			t = m_strpFindFile.find(".");
		}
	}
	int nStart = m_strpFindFile.length();
	while(nStart >= 0 && m_strpFindFile[nStart] != '\\' && m_strpFindFile[nStart] != '/')
	{
		--nStart;
	}
	++nStart;
	if(string::npos != t)
		t -= nStart;
	string strFileHead = m_strpFindFile.substr((size_t)nStart, t);
	if('\"' == strOutputPath[0] && '\"' == strOutputPath[strOutputPath.length() - 1])
	{
		strOutputPath = strOutputPath.substr(1, strOutputPath.length() - 2);
	}
	
	if(SLASH != strOutputPath[strOutputPath.length()-1])
	{
		strOutputPath += SLASH;
	}
	string strFile = strOutputPath + strFileHead + strTitle + ".pXBuild";
	
	GetLines(vResults,vSpectra,strTXTContent);
	ofstream of(strFile.c_str());
	of<<strTXTContent;
	of.close();
}


void CXLinkpXBuildReport::Close(void)
{
	m_Condition.clear();
}

bool CXLinkpXBuildReport::LoadFile(string strReportFile, vector<CXLinkMatchResult> & vResults , vector<CSpectrum> & vSpectra )
{
	vResults.clear();
	vSpectra.clear();
	return true;
}
