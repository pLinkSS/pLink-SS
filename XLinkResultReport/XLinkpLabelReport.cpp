#include <iostream>
#include <string>
 
#include "../include/predefine.h"
#include "../include/sdk.h"
#include "../include/interface.h"
#include "../include/option.h"

#include "XLinkResultReport.h"
#include "XLinkpLabelReport.h"
using namespace std;
using namespace proteomics_sdk;

CXLinkpLabelReport::CXLinkpLabelReport()
{
	char szBuf[1024] = {0};
	getcwd(szBuf, 1024);
    m_strWorkDir = szBuf;
    if(m_strWorkDir[m_strWorkDir.length() - 1] != SLASH)
    {
    	m_strWorkDir += SLASH;
    }
}

CXLinkpLabelReport::~CXLinkpLabelReport()
{
}

void CXLinkpLabelReport::Init(string strpFindFile,time_t  tmStartTime)
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
		CErrInfo info("CXLinkpLabelReport", "Init", "in the function CConditionReader::Read.");
		info.Append("strOption=" + strpFindFile);
		info.Append("m_strWorkDir=" + m_strWorkDir);
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CXLinkpLabelReport", "Init", "caught an unknown exception in the function CConditionReader::Read.");
		info.Append("strOption=" + strpFindFile);
		info.Append("m_strWorkDir=" + m_strWorkDir);
		throw runtime_error(info.Get().c_str());
	}
	
	m_Condition = reader.m_Condition;
    for(size_t i = 0;i < m_Condition.m_vSelectedFixMod.size();++i)
    	m_Condition.m_vSelectedVarMod.push_back(m_Condition.m_vSelectedFixMod[i]);
 
    COptionTool * pOption = new COptionTool("spectrum", m_strpFindFile.c_str());
    m_strSpectraPath = pOption->GetString("spec_path", "null");
    delete pOption;
}

string CXLinkpLabelReport::_2Uppercase(string str)
{
	string strRet = "";
	for(size_t i = 0;i < str.length() ; ++ i)
	{
		if(str[i] >= 'a' && str[i] <= 'z')
			strRet += str[i] - ('a'-'A');
		else
			strRet += str[i];
	}
	return strRet;
}

void CXLinkpLabelReport::GetLines(const vector<CXLinkMatchResult> & vResults , const vector<CSpectrum> & vSpectra ,string & strTXT)
{
	char szBuf[1024] = {'\0'};
	string strHead;
	strHead += "[FilePath]\n";
		
	// only DTA is supported 
	strHead += "File_Path=" + m_strSpectraPath + "\n";

	// todo : wrong for mod info
	strHead += "[Modification]\n";
	for(size_t i = 0;i < m_Condition.m_vSelectedVarMod.size();++i)
	{
		sprintf(szBuf,"%d=%s\n",i+1,m_Condition.m_vSelectedVarMod[i].m_strName.c_str());
		strHead += szBuf;
	}
		
	strHead += "[xlink]\n";
	strHead += "xlink=" + m_Condition.m_vSelectedXLinker[0].m_strName + "\n";
	strHead += "[Total]\n";

	
	strTXT.clear();
		
	size_t tCnt = 0;
	for(size_t i=0; i<vResults.size(); ++i)
	{
		if(vResults[i].m_vPeptideResults.empty())
			continue;	

		// choose linker
		/*
		if(vResults[i].m_vPeptideResults[0].m_XLink.m_nLinkerId !=1)
			continue;
		*/
		
		++tCnt;
		sprintf(szBuf, "[Spectrum%d]\n", tCnt);
		strTXT += szBuf;
		
		string strFilePath = _2Uppercase(vSpectra[i].m_strFilePath);
		sprintf(szBuf,"name=%s\n",strFilePath.c_str());
		strTXT += szBuf;
		
		if(vResults[i].m_vPeptideResults[0].m_XLink.m_eXLinkType == 0)
		{
			sprintf(szBuf,"pep1=%d ",vResults[i].m_vPeptideResults[0].m_XLink.m_eXLinkType);
		}
		else if(vResults[i].m_vPeptideResults[0].m_XLink.m_eXLinkType == 1)
		{
			sprintf(szBuf,"pep1=%d %d ",vResults[i].m_vPeptideResults[0].m_XLink.m_eXLinkType,vResults[i].m_vPeptideResults[0].m_XLink.m_tAlphaSite);
		}
		else
		{
			sprintf(szBuf,"pep1=%d %d %d ",vResults[i].m_vPeptideResults[0].m_XLink.m_eXLinkType,vResults[i].m_vPeptideResults[0].m_XLink.m_tAlphaSite,vResults[i].m_vPeptideResults[0].m_XLink.m_tBetaSite);
		}
		strTXT += szBuf;
		
		sprintf(szBuf,"%s %f ",vResults[i].m_vPeptideResults[0].m_AlphaPeptide.m_szSequence,vResults[i].m_vPeptideResults[0].m_lfScore);
		strTXT += szBuf;
		
		size_t tSite = 0;
		size_t tModId = 0;
		for(int j=0;j<vResults[i].m_vPeptideResults[0].m_AlphaPeptide.m_tModCnt;++j)
		{
			tModId = vResults[i].m_vPeptideResults[0].m_AlphaPeptide.m_tModSites[j][1];
			if(MT_PEP_NTERM == m_Condition.m_vSelectedVarMod[tModId].m_eModType || MT_PRO_NTERM == m_Condition.m_vSelectedVarMod[tModId].m_eModType)
				tSite = 0;
			else if(MT_PEP_CTERM == m_Condition.m_vSelectedVarMod[tModId].m_eModType || MT_PRO_CTERM == m_Condition.m_vSelectedVarMod[tModId].m_eModType)
				tSite = vResults[i].m_vPeptideResults[0].m_AlphaPeptide.m_tLength + 1;
			else
				tSite = vResults[i].m_vPeptideResults[0].m_AlphaPeptide.m_tModSites[j][0]; // 2014.7.16 Modification site no loner plus 1, since plused earlier
			
			//tSite = vResults[i].m_vPeptideResults[0].m_AlphaPeptide.m_tModSites[j][0]; 
			sprintf(szBuf,"%d,%d ",tSite,vResults[i].m_vPeptideResults[0].m_AlphaPeptide.m_tModSites[j][1] + 1);
			strTXT += szBuf;
		}
		
		if(vResults[i].m_vPeptideResults[0].m_bPair)
		{
			sprintf(szBuf,"%s %E ",vResults[i].m_vPeptideResults[0].m_BetaPeptide.m_szSequence,vResults[i].m_vPeptideResults[0].m_lfEvalue);
			strTXT += szBuf;
			
			size_t tBaseSite = vResults[i].m_vPeptideResults[0].m_AlphaPeptide.m_tLength + 3;
			
			for(int j=0;j<vResults[i].m_vPeptideResults[0].m_BetaPeptide.m_tModCnt;++j)
			{
				tModId = vResults[i].m_vPeptideResults[0].m_BetaPeptide.m_tModSites[j][1];
				if(MT_PEP_NTERM == m_Condition.m_vSelectedVarMod[tModId].m_eModType || MT_PRO_NTERM == m_Condition.m_vSelectedVarMod[tModId].m_eModType)
					tSite = tBaseSite + 0;
				else if(MT_PEP_CTERM == m_Condition.m_vSelectedVarMod[tModId].m_eModType || MT_PRO_CTERM == m_Condition.m_vSelectedVarMod[tModId].m_eModType)
					tSite = tBaseSite + 1 + vResults[i].m_vPeptideResults[0].m_BetaPeptide.m_tLength;
				else
					tSite = tBaseSite + vResults[i].m_vPeptideResults[0].m_BetaPeptide.m_tModSites[j][0]; // 2014.7.16 Modification site no loner plus 1, since plused earlier
				
				//tSite = vResults[i].m_vPeptideResults[0].m_BetaPeptide.m_tModSites[j][0];
				sprintf(szBuf,"%d,%d ",tSite,vResults[i].m_vPeptideResults[0].m_BetaPeptide.m_tModSites[j][1] + 1);
				strTXT += szBuf;
			}
		}
		strTXT += "\n";
	}

	sprintf(szBuf,"total=%d\n",tCnt);
	strHead += szBuf;
	
	strTXT = strHead + strTXT; 
}

bool CXLinkpLabelReport::LoadFile(string strReportFile, vector<CXLinkMatchResult> & vResults , vector<CSpectrum> & vSpectra )
{
	vResults.clear();
	vSpectra.clear();
	return true;
}

void CXLinkpLabelReport::WriteFile(const vector<CXLinkMatchResult> & vResults , const vector<CSpectrum> & vSpectra , string strOutputPath , string strTitle)
{
	string strTXTContent;
	GetLines(vResults,vSpectra, strTXTContent);
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
	
	string strFile = strOutputPath + strFileHead + strTitle + ".plabel";
	ofstream of(strFile.c_str());
	of<<strTXTContent;
	of.close();
}

void CXLinkpLabelReport::Close(void)
{
	m_Condition.clear();
}

