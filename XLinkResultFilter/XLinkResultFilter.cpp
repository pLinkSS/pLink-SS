#include <iostream>
#include <string>
 
#include "../include/predefine.h"
#include "../include/sdk.h"
#include "../include/interface.h"
#include "../include/option.h"

#include "../XLinkResultReport/XLinkResultReport.h"
#include "../XLinkResultReport/XLinkResultReportFactory.h"
#include "../XLinkPepResultFilter/XLinkPepResultFilter.h"
#include "FilterConf.h"
#include "XLinkResultFilterInterface.h"
#include "XLinkResultFilter.h"

using namespace std;
using namespace proteomics_sdk;


CXLinkResultFilter::CXLinkResultFilter()
{
}

CXLinkResultFilter::~CXLinkResultFilter()
{
}

void CXLinkResultFilter::_LoadInclusionList()
{
	m_mapInclusionList.clear();
	
	if(!m_filerconf.m_bInclusionListAvail)
		return ;
	
	FILE * fp;
	fp = fopen(m_filerconf.m_strInclusionList.c_str(),"r");
	if(fp == NULL)
	{
		m_filerconf.m_bInclusionListAvail = false;
		return;
	}
	
	char szbuf[1024];
	while(fgets(szbuf,1024,fp))
	{
		if(szbuf[strlen(szbuf)-1]==0x0a)
			szbuf[strlen(szbuf)-1]=0;
		if(szbuf[strlen(szbuf)-1]==0x0d)
			szbuf[strlen(szbuf)-1]=0;
		
		m_mapInclusionList[szbuf] = true;
	}
	
	fclose(fp);
	return;
}

void CXLinkResultFilter::Init(string strConf)
{
	m_filerconf.Load(strConf);
	m_pepResultFilter.SetFilterStandand(m_filerconf.m_pepFilter);
	_LoadInclusionList();
	
	m_filerconf.Output();
	
}

char CXLinkResultFilter::_GetReverseTag(const CXLinkPepResult & pep_res)
{
	int nForward = 0;
	for(size_t i = 0;i < pep_res.m_vAlphaProteinAC.size() ; ++i)
	{
		if(strncmp ( pep_res.m_vAlphaProteinAC[i].c_str(), m_filerconf.m_strReverseTag.c_str(),m_filerconf.m_strReverseTag.length()))
		{
			nForward ++;
			break;
		}
	}
	if(pep_res.m_bPair)
	{
		for(size_t i = 0;i < pep_res.m_vBetaProteinAC.size() ; ++i)
		{
			if(strncmp ( pep_res.m_vBetaProteinAC[i].c_str(), m_filerconf.m_strReverseTag.c_str(),m_filerconf.m_strReverseTag.length()))
			{
				nForward ++;
				break;
			}
		}
	}
	if(nForward == 0)
		// forward = 0
		return 'F';
	else if(nForward == 1)
	{
		// forward = 1
		if(pep_res.m_bPair)
			return 'U';
		else
			return 'T';
	}
	else
		// forward = 2
		return 'T';

}

void CXLinkResultFilter::_FilterByFDR(vector<CXLinkMatchResult> & vResults, vector<CSpectrum> & vSpectra,int nTolBase)
{ 
	if(m_vBeOutput.size()!=vResults.size())
		return;
	
	const int CHARGE_SIZE = 20;
	vector<vector<pair<double,char> > > vAllChgEvaHist;
	
	vAllChgEvaHist.clear();
	for(int i=0;i<CHARGE_SIZE;++i)
	{
		vAllChgEvaHist.push_back(vector<pair<double,char> >());
	}
	
	int nCharge;
	double lfEvalue;
	bool bTrue;
	string strProAC;
	
	for(size_t i=0;i<vResults.size();++i)
	{
		if(vResults[i].m_vPeptideResults.size() && m_vBeOutput[i] == nTolBase)
		{
			nCharge = vSpectra[i].m_nCharge;
			if(nCharge < 0)
			{
				continue;
			}
			else if(nCharge > (CHARGE_SIZE-2) )
			{
				nCharge = CHARGE_SIZE-1;
			}
			lfEvalue = vResults[i].m_vPeptideResults[0].m_lfEvalue;
			
			char cTag = _GetReverseTag(vResults[i].m_vPeptideResults[0]);
			vAllChgEvaHist[nCharge].push_back(pair<double,char>(lfEvalue,cTag));
		}
	}
	
	vector<double> vEvalueThreashold;
	vEvalueThreashold.clear();
	
	for(int i=0;i<CHARGE_SIZE;++i)
	{
		if(vAllChgEvaHist[i].size()==0)
		{
			vEvalueThreashold.push_back(1);
		}
		else
		{
			vEvalueThreashold.push_back(_GetEvaThreashold(vAllChgEvaHist[i]));
		}
	}
	
	/*
	cout << "======Filter by FDR Begin=====" << endl
	<< "Tol Base :" << nTolBase << endl
	<< "Evalue threshold : " << endl;
	*/
	for(int i=0;i<CHARGE_SIZE;++i)
	{
		cout << "charge " << i << " " << vEvalueThreashold[i] << endl
		<< "size " << vAllChgEvaHist[i].size() << endl;
	}
	//cout << "======Filter by FDR End=====" << endl;
	
	for(size_t i=0;i<vResults.size();++i)
	{
		if(vResults[i].m_vPeptideResults.size() && m_vBeOutput[i] == nTolBase)
		{
			if(_GetReverseTag(vResults[i].m_vPeptideResults[0])!='T')
			{
				m_vBeOutput[i] = -1;
				continue;
			}
			
			nCharge = vSpectra[i].m_nCharge;
			if(nCharge < 0)
			{
				continue;
			}
			else if(nCharge > (CHARGE_SIZE-2) )
			{
				nCharge = CHARGE_SIZE-1;
			}
			
			if(vResults[i].m_vPeptideResults[0].m_lfEvalue > vEvalueThreashold[nCharge])
				m_vBeOutput[i] = -1;
		}
	}
	
}


double CXLinkResultFilter::_GetEvaThreashold(vector<pair<double,char> > & vEvaHist)
{
	if(vEvaHist.size()<=0)
		return 1;
	
	sort(vEvaHist.begin(), vEvaHist.end(), EVLesser);

	double lfThreshold = vEvaHist[0].first;
	double lfCurFDR = 0;
	double lfFDR = m_filerconf.m_lfFDR;
	
	int nTCnt = 0;
	int nUCnt = 0;
	int nFCnt = 0;
	int i ;
	for(i=0 ; i<vEvaHist.size(); ++i)
	{
		if(vEvaHist[i].second == 'T')
			nTCnt ++;
		else if(vEvaHist[i].second == 'U')
			nUCnt ++;
		else if(vEvaHist[i].second == 'F')
			nFCnt ++ ;

		int nFalseCount = 0;
		if(nUCnt - 2*nFCnt > 0)
			nFalseCount = nUCnt - nFCnt;
		else
			nFalseCount = nFCnt ;
		
		if(nTCnt)
			lfCurFDR = (double)(nFalseCount) / (double) nTCnt;
		else
			lfCurFDR = 1;

		/*
		if(lfFDR > lfCurFDR && vEvaHist[i].second && lfThreshold < vEvaHist[i].first)
			lfThreshold = vEvaHist[i].first;
		*/
		
		if(lfCurFDR > lfFDR)
			break;
	}
	i--;
	if(i<0)
		i=0;
	lfThreshold = vEvaHist[i].first;
	/*
	cout << i << endl
		<< "TrueCnt = " << nTrueCnt << endl
		<< "threshold = " << lfThreshold << endl;
	*/
	
	if(0 >= lfThreshold)
		lfThreshold = 1;
	return lfThreshold;
}

void CXLinkResultFilter::Close()
{
	m_vBeOutput.clear();
	vector<char>().swap(m_vBeOutput);
}

int CXLinkResultFilter::_IsInTolWindow(double lfExpMH,double lfCalMH)
{
	double lfTolBase;
	double lfTol;
	if(m_filerconf.m_strTolType == "ppm")
	{
		lfTol = lfCalMH * 0.000001 * m_filerconf.m_lfTol;
	}
	else
		lfTol = m_filerconf.m_lfTol;
	
	for(size_t i = 0 ;i < m_filerconf.m_vTolBase.size(); ++ i)
	{
		if(m_filerconf.m_strTolBaseType == "ppm")
		{
			lfTolBase = lfCalMH * 0.000001 * m_filerconf.m_vTolBase[i];
		}
		else
			lfTolBase = m_filerconf.m_vTolBase[i];
		
		double lfTolLB ,lfTolUB;	
		lfTolLB = lfCalMH + lfTolBase - lfTol;
		lfTolUB = lfCalMH + lfTolBase + lfTol;
		
		if(lfExpMH <= lfTolUB && lfExpMH >= lfTolLB)
			return i;
	}
	return -1;
}

void CXLinkResultFilter::Run()
{
	CXLinkResultReportFactory reportFactory;
	CXLinkResultReport * report = reportFactory.GetReport(Xlink_STANDARD);
	if(report==NULL)
	{
		cout << "can't get the report" << endl;
		return ;
	}
	report->Init(m_filerconf.m_strpFindFile);

	vector<CXLinkMatchResult> vResults;
	vector<CSpectrum> vSpectra;

	cout << "loading from input protein report.." << endl;
	bool bOK = report->LoadFile(m_filerconf.m_strInputFile,vResults,vSpectra);

	if(bOK)
	{
		/*
		cout << "loading complete .." << endl
			<< "result num = " << vResults.size() << endl;
		*/
	}
	else
	{
		cout << "errors occur while loading .. " << endl;
		cout << m_filerconf.m_strInputFile<<endl;
		report->Close();
		return;
	}
	size_t tResultSize = 0;
	tResultSize = vResults.size();

	m_vBeOutput.clear();
	vector<char>().swap(m_vBeOutput);
	for(size_t i = 0 ;i < tResultSize; ++ i)
	{
		if(vResults[i].m_vPeptideResults.size()<=0)
		{
			m_vBeOutput.push_back(-1);
			continue;
		}
		
		if(m_filerconf.m_bInclusionListAvail)
		{
			if(m_mapInclusionList.find(vSpectra[i].m_strFilePath)== m_mapInclusionList.end())
			{
				m_vBeOutput.push_back(-1);
				continue;
			}
		}
		
		if(vResults[i].m_vPeptideResults[0].m_XLink.m_eXLinkType != m_filerconf.m_nXLinkType)
		{
			m_vBeOutput.push_back(-1);
			continue;
		}
		
		if(m_filerconf.m_nLinkerId >=0 && vResults[i].m_vPeptideResults[0].m_XLink.m_nLinkerId != m_filerconf.m_nLinkerId)
		{
			m_vBeOutput.push_back(-1);
			continue;
		}
		
		const CXLinkPepResult & pep_res = vResults[i].m_vPeptideResults[0];
		// filter by standard
		
		if(false == m_pepResultFilter.Filter(pep_res))
		{
			m_vBeOutput.push_back(-1);
			continue;
		}
				
		// filter by mass range
		int nRet = _IsInTolWindow(vSpectra[i].m_lfMH,vResults[i].m_vPeptideResults[0].m_lfCalc_MH);
		m_vBeOutput.push_back(nRet);
		
		
	}
	
	//cout << "filtering in the context of peptides by FDR and evalue (relative filtering)" << endl;
	if(m_filerconf.m_lfFDR >= 0)
	{
		for(size_t i = 0 ;i < m_filerconf.m_vTolBase.size(); ++ i)
			_FilterByFDR(vResults,vSpectra,i);
	}

	for(size_t i = 0 ;i < tResultSize; ++ i)
	{
		if(vResults[i].m_vPeptideResults.size()<=0)
		{
			continue;
		}
		
		if(m_vBeOutput[i] == -1)
		{
			vResults[i].m_vPeptideResults.clear();
			vector<CXLinkPepResult>().swap(vResults[i].m_vPeptideResults);
		}
	}
	
	//cout << "saving to output protein report ..." << endl;
	string strTitle = "_" + m_filerconf.m_strTitle;
	report->WriteFile(vResults,vSpectra,m_filerconf.m_strOutputPath,strTitle);
	report->Close();
	//cout << "saving complete .. " << endl;
	
	vResults.clear();
	vector<CXLinkMatchResult>().swap(vResults);
	vSpectra.clear();
	vector<CSpectrum>().swap(vSpectra);
	
	delete report;
	
}


	
