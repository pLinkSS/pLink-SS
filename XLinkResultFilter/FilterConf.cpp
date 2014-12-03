#include <string>
#include <iostream>

#include "../include/predefine.h"
#include "../include/sdk.h"
#include "../include/interface.h"
#include "../include/option.h"
//#include "../ConditionReader/ConditionReader.h"
#include "../XLinkPepResultFilter/XLinkPepResultFilter.h"
#include "FilterConf.h"

using namespace std;
CFilterConf::CFilterConf()
{
	// get working dir
	char szBuf[1024] = {0};
	getcwd(szBuf, 1024);
    m_strWorkDir = szBuf;
    if(m_strWorkDir[m_strWorkDir.length() - 1] != SLASH)
    {
    	m_strWorkDir += SLASH;
    }

    m_lfTol = 0;
	m_vTolBase.clear();
	
	m_strInclusionList = "";
	m_bInclusionListAvail = false;
	
	m_nLinkerId= -1;

}

CFilterConf::~CFilterConf()
{
	m_Condition.clear();
}

void CFilterConf::Load(string strConfigFile)
{
	m_strOption = strConfigFile;
	
	_LoadIO();
	_LoadFilter();
}

void CFilterConf::Output()
{
	/*
	cout << "=======IN XLINK RESULT FILTER =======" << endl
	<< "pfind : " << m_strpFindFile << endl
	<< "input file : " << m_strInputFile << endl
	<< "output file : " << m_strOutputPath << endl
	<< "xlink type : " << m_nXLinkType << endl
	<< "linker id : " << m_nLinkerId << endl
	<< "precursor tolerance : " << m_lfTol << endl
	<< "precursor tolerance type : " << m_strTolType << endl;
	for(size_t i = 0 ; i < m_vTolBase.size() ; ++ i)
	{
		cout << "precursor tolerance base " << i << " : " << m_vTolBase[i] << endl;
	}
	cout << "precursor tolerance base type : " << m_strTolBaseType << endl
	<< "fdr : " << m_lfFDR << endl
	<< "reverse tag " << m_strReverseTag << endl
	<< "Inclusion List " << m_strInclusionList << endl;
	m_pepFilter.output();
	*/
}

void CFilterConf::_LoadIO()
{
	COptionTool * pOption(NULL);
    pOption = new COptionTool("IO", m_strOption.c_str());
    if(!pOption)
    {
    	cout << "can't open pOption file : " << m_strOption << " on [IO]" << endl;
    }

    m_strpFindFile = pOption->GetString("pFindFile", "");
    m_strInputFile = pOption->GetString("InputReportFile", "");
	m_strOutputPath = pOption->GetString("OutputPath", "");
	m_strTitle= pOption->GetString("title","filter");

	CConditionReader reader(m_strpFindFile, m_strWorkDir);

    try
    {
    	reader.Read();
    }
    catch(exception & e)
    {
    	CErrInfo info("CFilterConf", "_LoadIO", 
    			"in the function CConditionReader::Read");
    	info.Append("m_strParamFile=" + m_strpFindFile);
    	info.Append("m_strWorkDir=" + m_strWorkDir);
    	throw runtime_error(info.Get(e).c_str());
    }
    catch(...)
    {
    	CErrInfo info("CFilterConf", "_LoadIO", 
    			"caught an unknown exception in the function CConditionReader::Read");
    	info.Append("m_strParamFile=" + m_strpFindFile);
    	info.Append("m_strWorkDir=" + m_strWorkDir);
    	throw runtime_error(info.Get().c_str());
    }
    
	m_Condition = reader.m_Condition;
	for(size_t i = 0;i < m_Condition.m_vSelectedFixMod.size();++i)
	   m_Condition.m_vSelectedVarMod.push_back(m_Condition.m_vSelectedFixMod[i]);
	
    delete pOption;
    
}


void CFilterConf::_LoadFilter()
{
	COptionTool * pOption(NULL);
	pOption = new COptionTool("filter", m_strOption.c_str());
	if(!pOption)
    {
    	cout << "can't open pOption file : " << m_strOption << " on [filter]" << endl;
    }
	
	string strTmp;
	m_strTolType = pOption->GetString("PepTolType", "ppm");
	strTmp = pOption->GetString("PepTol", "50");
	m_lfTol = atof(strTmp.c_str());
	m_strTolBaseType = pOption->GetString("PepTolBaseType", "Da");
	int nTmp = pOption->GetInteger("PepTolBaseTotal", 1);
	char szbuf[128];
	m_vTolBase.clear();
	for(int i = 0 ;i < nTmp ; ++i)
	{
		sprintf(szbuf,"PepTolBase%d",i);
		strTmp = pOption->GetString(szbuf, "0");
		m_vTolBase.push_back(atof(strTmp.c_str()));
	}
	
	m_strReverseTag = pOption->GetString("ReverseTag", "REVERSE");
	strTmp = pOption->GetString("FDR", "0.01");
	m_lfFDR = atof(strTmp.c_str());
	
	m_nXLinkType = pOption->GetInteger("XLinkType", 3);
	
	// retain results with either linker as default
	m_nLinkerId = pOption->GetInteger("LinkerId", -1);
	
	strTmp = pOption->GetString("MaxEvalue", "1");
	m_pepFilter.lfEvalue = atof(strTmp.c_str());
	strTmp = pOption->GetString("MinScore", "0");
	m_pepFilter.lfScore = atof(strTmp.c_str());
	m_pepFilter.nPepLength = pOption->GetInteger("MinPepLength", 0);
	m_pepFilter.nContiAAnum = pOption->GetInteger("MinContiAAnum", 0);
	m_pepFilter.nAAConfLevel = pOption->GetInteger("AAConfLevel", 1);
	strTmp = pOption->GetString("MinMatchOverUnMatch", "0.0");
	m_pepFilter.lfMatchOverUnMatch = atof(strTmp.c_str());
	m_pepFilter.bCterm = pOption->GetInteger("SaveCTerm", 1);
	
	strTmp = pOption->GetString("InclusionList", "");
	m_strInclusionList = "";
	FILE * fp;
	fp = fopen(strTmp.c_str(),"r");
	if(fp == NULL)
	{
		m_bInclusionListAvail = false;
	}
	else
	{
		m_bInclusionListAvail = true;
		m_strInclusionList = strTmp;
		fclose(fp);
	}
	
	
	delete pOption;
}