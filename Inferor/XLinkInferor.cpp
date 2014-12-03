#ifdef WIN32
#include <direct.h>
#endif
#include <iostream>
#include "../include/predefine.h"
#include "../include/sdk.h"
#include "../include/option.h"
#include "../include/interface.h"
#include "../XLinkPFDIO/XLinkPFDFileIO.h"
#include "../XLinkProteinInfer/XLinkACInfer.h"
#include "../XLinkResultReport/XLinkResultReport.h"
#include "../XLinkResultReport/XLinkResultReportFactory.h"
#include "../ProteinIndex/ProteinIndex.h"
#include "../ProteinIndex/ProteinReaderFactory.h"
#include "../ProteinIndex/ProteinWriterFactory.h"
#include "../ProteinIndex/ProteinReader.h"
#include "../ProteinIndex/ProteinWriter.h"
#include "Inferor.h"
#include "XLinkInferor.h"
using namespace ProteinIndex;
using namespace proteomics_sdk;
using namespace std;
//#define _DEBUG
#define PROTEININFER

CXLinkInferor::CXLinkInferor()
{
	char szBuf[BUFFER_SIZE] = {0};
	getcwd(szBuf, BUFFER_SIZE);
    m_strWorkDir = szBuf;
    if(m_strWorkDir[m_strWorkDir.length() - 1] != SLASH)
    {
    	m_strWorkDir += SLASH;
    }

    m_pTrace = CTrace::GetInstance();
    m_pTrace->Debug("CXLinkInferor::CXLinkInferor()");
}

CXLinkInferor::~CXLinkInferor()
{
}

bool CXLinkInferor::Init(string strOption,time_t tmStartTime)
{
	m_strOption = strOption;
	m_tmStartTime = tmStartTime;
	
	m_pTrace->Info("Inferor initializing...");
	CConditionReader reader(strOption, m_strWorkDir);
	try
	{
		reader.Read();
	}
	catch(exception & e)
	{
		CErrInfo info("CXLinkInferor", "Init", "in the function CConditionReader::Read.");
		info.Append("strOption=" + strOption);
		info.Append("m_strWorkDir=" + m_strWorkDir);
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CXLinkInferor", "Init", "caught an unknown exception in the function CConditionReader::Read.");
		info.Append("strOption=" + strOption);
		info.Append("m_strWorkDir=" + m_strWorkDir);
		throw runtime_error(info.Get().c_str());
	}
	
	m_Condition = reader.m_Condition;
    for(size_t i = 0;i < m_Condition.m_vSelectedFixMod.size();++i)
    	m_Condition.m_vSelectedVarMod.push_back(m_Condition.m_vSelectedFixMod[i]);
    
#ifdef _DEBUG
	CConditionReader::Test(m_Condition);
#endif
	
	
    COptionTool * pOption = new COptionTool("spectrum", m_strOption.c_str());

    m_strSpectraPath = pOption->GetString("spec_path", "null");
    if("null" == m_strSpectraPath)
    {
        delete pOption;
        return false;
    }
    
    delete pOption;
    
    return true;
	
}

void CXLinkInferor::_WritePfd2PBuildConfigFile(string strParams, string strOutputPath, int nFileTotal, int nSpectraTotal, string strIdentifier)
{
	FILE * fp;
	
	size_t t = strParams.find(".pfind");
	if(string::npos == t)
	{
		t = strParams.find(".PFIND");
		if(string::npos == t)
		{
			t = strParams.find(".");
		}
	}
	int nStart = strParams.length();
	while(nStart >= 0 && strParams[nStart] != '\\' && strParams[nStart] != '/')
	{
		--nStart;
	}
	++nStart;
	if(string::npos != t)
		t -= nStart;
	string strFile = strParams.substr((size_t)nStart, t);
	
	strFile = strFile + ".pfd2pbuild";
	
	if( (fp = fopen(strFile.c_str(),"w")) == NULL )
	{
		return ;
	}
	
	fprintf(fp,"[pfd2pbuild]\n");
	fprintf(fp,"pFindFile=%s\n",strParams.c_str());
	fprintf(fp,"OutputPath=%s\n",strOutputPath.c_str());
	fprintf(fp,"TmpFileTotal=%d\n",nFileTotal);
	fprintf(fp,"TmpFileTag=%s\n",strIdentifier.c_str());
	fprintf(fp,"SpectraTotal=%d\n",nSpectraTotal);
	
	fclose(fp);

	
}
void CXLinkInferor::_WritePfd2PXBuildConfigFile(string strParams, string strOutputPath, int nFileTotal, int nSpectraTotal, string strIdentifier)
{
	FILE * fp;
	
	size_t t = strParams.find(".pfind");
	if(string::npos == t)
	{
		t = strParams.find(".PFIND");
		if(string::npos == t)
		{
			t = strParams.find(".");
		}
	}
	int nStart = strParams.length();
	while(nStart >= 0 && strParams[nStart] != '\\' && strParams[nStart] != '/')
	{
		--nStart;
	}
	++nStart;
	if(string::npos != t)
		t -= nStart;
	string strFile = strParams.substr((size_t)nStart, t);
	
	strFile = strFile + ".pfd2pxbuild";
	
	if( (fp = fopen(strFile.c_str(),"w")) == NULL )
	{
		return ;
	}
	
	fprintf(fp,"[pfd2pxbuild]\n");
	fprintf(fp,"pFindFile=%s\n",strParams.c_str());
	fprintf(fp,"OutputPath=%s\n",strOutputPath.c_str());
	fprintf(fp,"TmpFileTotal=%d\n",nFileTotal);
	fprintf(fp,"TmpFileTag=%s\n",strIdentifier.c_str());
	fprintf(fp,"SpectraTotal=%d\n",nSpectraTotal);
	
	fclose(fp);

	
}

void CXLinkInferor::_WritePfd2PlabelConfigFile(string strParams, string strOutputPath, int nFileTotal, int nSpectraTotal, string strIdentifier)
{
	FILE * fp;
	
	size_t t = strParams.find(".pfind");
	if(string::npos == t)
	{
		t = strParams.find(".PFIND");
		if(string::npos == t)
		{
			t = strParams.find(".");
		}
	}
	int nStart = strParams.length();
	while(nStart >= 0 && strParams[nStart] != '\\' && strParams[nStart] != '/')
	{
		--nStart;
	}
	++nStart;
	if(string::npos != t)
		t -= nStart;
	string strFile = strParams.substr((size_t)nStart, t);
	
	strFile = strFile + ".pfd2plabel";
	
	if( (fp = fopen(strFile.c_str(),"w")) == NULL )
	{
		return ;
	}
	
	fprintf(fp,"[pfd2plabel]\n");
	fprintf(fp,"pFindFile=%s\n",strParams.c_str());
	fprintf(fp,"OutputPath=%s\n",strOutputPath.c_str());
	fprintf(fp,"TmpFileTotal=%d\n",nFileTotal);
	fprintf(fp,"TmpFileTag=%s\n",strIdentifier.c_str());
	fprintf(fp,"SpectraTotal=%d\n",nSpectraTotal);
	fprintf(fp,"XLinkType=3\n");
	fprintf(fp,"MaxEvalue=1.0\n");
	fprintf(fp,"MinScore=0.0\n");
	fprintf(fp,"MinPepLength=4\n");
	fprintf(fp,"MinContiAAnum=3\n");
	fprintf(fp,"AAConfLevel=1\n");
	fprintf(fp,"MinMatchOverUnMatch=0.0\n");
	
	fclose(fp);
	
}

void CXLinkInferor::Run(string strParams, string strOutputPath, int nFileTotal, int nSpectraTotal, string strIdentifier)
{
	
#ifdef _DEBUG
	cout << strOutputPath << endl << nSpectraTotal << endl
		<< nFileTotal << endl << strIdentifier << endl;
#endif
	
	vector<CXLinkMatchResult> vResults;
	vector<CSpectrum> vSpectra;
	CXLinkPFDFileIO pfd;
	size_t tSize = 0;

	pfd.LoadAll(strOutputPath, nSpectraTotal, nFileTotal, vResults, vSpectra, tSize, strIdentifier);

	vector<CAssignedProtein> vProteins;

#ifdef PROTEININFER

	CXLinkACInfer proteinInfer;
	
	try{
		FILTER_CRITERIA_INFO cr = {m_Condition.m_lfFPRThreshold, 800, 5, 1, 1};
		m_pTrace->Debug("proteinInfer.Initialize(&m_Condition)");
		proteinInfer.Initialize(&m_Condition);
		m_pTrace->Debug("proteinInfer.SetSign(m_Condition.m_strFPSign)");
		proteinInfer.SetSign(m_Condition.m_strFPSign);
		m_pTrace->Debug("proteinInfer.Infer(vResults, vSpectra, cr)");

		vProteins = proteinInfer.Infer(vResults, vSpectra, cr);
		m_pTrace->Debug("end proteinInfer.Infer(vResults, vSpectra, cr)");
	}
	catch(exception & e)
	{
		CErrInfo info("CXLinkInferor", "Run", "in the function CXLinkACInfer::Infer");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CXLinkInferor", "Run", "in the function CXLinkACInfer::Infer");
		throw runtime_error(info.Get().c_str());
	}

#endif

//	cout << "Inferring completed." << endl;
//	cout << "Protein Num: " << vProteins.size() << endl;
	
	// todo : protein report
	CXLinkResultReportFactory reportFactory;
	CXLinkResultReport * report = reportFactory.GetReport(Xlink_STANDARD);
	
	try{
		report->Init(m_strOption,m_tmStartTime);
		report->WriteFile(vResults,vSpectra,strOutputPath,"");
		report->Close();
		delete report;
	}
	catch(exception & e)
	{
		
		CErrInfo info("CXLinkInferor", "Run", "in the function CXLinkProteinReport::Init");
		info.Append("strOutputPath=" + strOutputPath);
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		
		CErrInfo info("CXLinkInferor", "Run", "in the function CXLinkProteinReport::Init");
		info.Append("strOutputPath=" + strOutputPath);
		throw runtime_error(info.Get().c_str());
	}
}

