
#include <iostream>
#include <string>
 
#include "../include/predefine.h"
#include "../include/sdk.h"
#include "../include/interface.h"
#include "../include/option.h"

#include "XLinkResultReport.h"
#include "XLinkPepReport.h"


using namespace std;
using namespace proteomics_sdk;

CXLinkPepReport::CXLinkPepReport()
{
	char szBuf[1024] = {0};
	getcwd(szBuf, 1024);
    m_strWorkDir = szBuf;
    if(m_strWorkDir[m_strWorkDir.length() - 1] != SLASH)
    {
    	m_strWorkDir += SLASH;
    }
}

CXLinkPepReport::~CXLinkPepReport()
{
}

void CXLinkPepReport::Init(string strpFindFile,time_t tmStartTime)
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
		CErrInfo info("CXLinkPepReport", "Init", "in the function CConditionReader::Read.");
		info.Append("strOption=" + strpFindFile);
		info.Append("m_strWorkDir=" + m_strWorkDir);
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CXLinkPepReport", "Init", "caught an unknown exception in the function CConditionReader::Read.");
		info.Append("strOption=" + strpFindFile);
		info.Append("m_strWorkDir=" + m_strWorkDir);
		throw runtime_error(info.Get().c_str());
	}
	
	m_Condition = reader.m_Condition;
		
}

void CXLinkPepReport::GetLines(const vector<CXLinkMatchResult> & vResults , const vector<CSpectrum> & vSpectra ,string & strTXT)
{
	strTXT = "" ;
	map<string ,double> mapPepsScore;
	for(size_t i = 0 ;i < vResults.size(); ++i)
	{
		if(vResults[i].m_vPeptideResults.empty())
				continue;
		
		const CXLinkPepResult & pep_res = vResults[i].m_vPeptideResults[0];
		
		int nReverse = 0; 
		for(size_t i = 0;i < pep_res.m_vAlphaProteinAC.size() ; ++i)
		{
			if(strncmp ( pep_res.m_vAlphaProteinAC[i].c_str(), m_Condition.m_strFPSign.c_str(),m_Condition.m_strFPSign.length())== 0)
			{
				nReverse ++;
				break;
			}
		}
		if(pep_res.m_bPair)
		{
			for(size_t i = 0;i < pep_res.m_vBetaProteinAC.size() ; ++i)
			{
				if(strncmp ( pep_res.m_vBetaProteinAC[i].c_str(), m_Condition.m_strFPSign.c_str(),m_Condition.m_strFPSign.length())== 0)
				{
					nReverse ++;
					break;
				}
			}
		}
		
		if(nReverse > 0)
			continue;
		
		
		string strPep = _getPepString(pep_res);
		char szbuf[10];
		sprintf(szbuf,"	%d+",vSpectra[i].m_nCharge);
		strPep += szbuf;
				
		if(mapPepsScore.find(strPep) == mapPepsScore.end() )
		{
			mapPepsScore[strPep] = pep_res.m_lfScore;
		}
		else
		{
			if(mapPepsScore[strPep] < pep_res.m_lfScore)
			{
				mapPepsScore[strPep] = pep_res.m_lfScore;
			}
		}
	}
	
	map<string, double >::iterator iter;
	char szbuf[1024];
	
	for(iter = mapPepsScore.begin();iter!= mapPepsScore.end(); ++iter)
	{
		sprintf(szbuf,"%s	%f\n",iter->first.c_str(),iter->second);
		strTXT += szbuf;
	}
	
}

bool CXLinkPepReport::LoadFile(string strReportFile, vector<CXLinkMatchResult> & vResults , vector<CSpectrum> & vSpectra )
{
	return true;
}

string CXLinkPepReport::_getPepString(const CXLinkPepResult & pep_res)
{
	string strPep;
	char szbuf[10];
	if(pep_res.m_XLink.m_eXLinkType < 0)
		return strPep;
	strPep += pep_res.m_AlphaPeptide.m_szSequence;
	if(pep_res.m_XLink.m_eXLinkType < 1)
		return strPep;
	sprintf(szbuf,"(%d)",pep_res.m_XLink.m_tAlphaSite);
	strPep += szbuf;
	if(pep_res.m_XLink.m_eXLinkType < 2)
		return strPep;
	if(pep_res.m_XLink.m_eXLinkType == 2)
	{
		sprintf(szbuf,"(%d)",pep_res.m_XLink.m_tBetaSite);
		strPep += szbuf;
	}
	else if(pep_res.m_XLink.m_eXLinkType == 3)
	{
		strPep += "-";
		strPep += pep_res.m_BetaPeptide.m_szSequence;
		sprintf(szbuf,"(%d)",pep_res.m_XLink.m_tBetaSite);
		strPep += szbuf;
		size_t tSite ;
		if(strcmp(pep_res.m_BetaPeptide.m_szSequence,pep_res.m_AlphaPeptide.m_szSequence)==0 && pep_res.m_XLink.m_tAlphaSite > pep_res.m_XLink.m_tBetaSite)
		{
			strPep = pep_res.m_AlphaPeptide.m_szSequence;
			sprintf(szbuf,"(%d)",pep_res.m_XLink.m_tBetaSite);
			strPep += szbuf;
			strPep += "-";
			strPep += pep_res.m_BetaPeptide.m_szSequence;
			sprintf(szbuf,"(%d)",pep_res.m_XLink.m_tAlphaSite);
			strPep += szbuf;
		}
	}
	
	return strPep;
}

void CXLinkPepReport::WriteFile(const vector<CXLinkMatchResult> & vResults , const vector<CSpectrum> & vSpectra , string strOutputPath , string strTitle)
{

	string strTXTContent;
	GetLines(vResults,vSpectra,strTXTContent);

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
	string strFile = strOutputPath + strFileHead + strTitle + ".peptide";
	ofstream of(strFile.c_str());
	of<<strTXTContent;
	of.close();
}

void CXLinkPepReport::Close(void)
{
	
}


