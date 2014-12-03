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
#include "XLinkResultFDRFilter.h"

using namespace std;
using namespace proteomics_sdk;


CXLinkResultFDRFilter::CXLinkResultFDRFilter()
{
}

CXLinkResultFDRFilter::~CXLinkResultFDRFilter()
{
}

void CXLinkResultFDRFilter::Init(string strConf)
{
	m_filerconf.Load(strConf);
	m_pepResultFilter.SetFilterStandand(m_filerconf.m_pepFilter);
	m_filerconf.Output();
}


char CXLinkResultFDRFilter::_GetReverseTag(const CXLinkPepResult & pep_res)
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
			return 'F';
	}
	else
		// forward = 2
		return 'T';
	
}

void CXLinkResultFDRFilter::Close()
{
	m_vBeOutput.clear();
	vector<int>().swap(m_vBeOutput);
	m_vTUFTag.clear();
	vector<char>().swap(m_vTUFTag);
	m_vFeatureList.clear();
	vector<pair<double,size_t> >().swap(m_vFeatureList);
}

void CXLinkResultFDRFilter::_SetTUFTag()
{
	m_vTUFTag.clear();
	for(size_t i = 0 ;i < m_vResults.size(); ++i)
	{
		m_vTUFTag.push_back(0);
	}
	
	for(size_t i = 0 ;i < m_vResults.size(); ++ i)
	{
		if(m_vResults[i].m_vPeptideResults.size() <= 0 || m_vBeOutput[i] == -1)
			continue;
		const CXLinkPepResult & pep_res = m_vResults[i].m_vPeptideResults[0];
		m_vTUFTag[i] = _GetReverseTag(pep_res);
	}
}

void CXLinkResultFDRFilter::_EstimateTargetNum()
{
	int nCountT = 0, nCountU = 0, nCountF = 0;
	 
	for(size_t i = 0 ;i < m_vResults.size(); ++ i)
	{
		if(m_vTUFTag[i] == 'T')
		{
			nCountT ++ ;
		}
		else if(m_vTUFTag[i] == 'U')
		{
			nCountU ++ ;
		}
		else if(m_vTUFTag[i] == 'F')
		{
			nCountF ++ ;
		}
	}
	
	cout << "Total num = " << nCountT + nCountU + nCountF << endl
		<< "T = " << nCountT << endl
		<< "U = " << nCountU << endl
		<< "F = " << nCountF << endl;
	int nFalseCount = nCountU - nCountF;
	if(nCountU < 2*nCountF)
		nFalseCount = nCountF;
	else
		nFalseCount = nCountU - nCountF;
	m_nTargetNum = nCountT - nFalseCount;
	cout << "estimated fdr = " << (nCountT!=0 ? ((nFalseCount + 0.0) / nCountT):1) << endl;
	if(m_nTargetNum <= 0)
		m_nTargetNum = 0;
	
	char tmpch;
	//cin >> tmpch;

}

int CXLinkResultFDRFilter::_GetTagLength(const CXLinkPepResult & pep_res)
{
	int nMaxLen = 0 , nCurLen = 0;
	for(size_t i=0;i<pep_res.m_AlphaPeptide.m_tLength;++i)
	{
		
		if(pep_res.m_stMatchInfo.aPepConf[i] >= m_filerconf.m_pepFilter.nAAConfLevel)
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
	
	int nFinalMaxLen = nMaxLen;
	
	if(pep_res.m_bPair)
	{
		nMaxLen = 0;
		nCurLen = 0;
		for(size_t i=0;i<pep_res.m_BetaPeptide.m_tLength;++i)
		{
			
			if(pep_res.m_stMatchInfo.aPepConf[pep_res.m_AlphaPeptide.m_tLength + i] >=  m_filerconf.m_pepFilter.nAAConfLevel)
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
		
		if(nFinalMaxLen > nMaxLen)
			nFinalMaxLen = nMaxLen;
	}
	return nFinalMaxLen;
}

int CXLinkResultFDRFilter::_CalCharge(bool bSet )
{
	m_vFeatureList.clear();
	m_vFeatureList.reserve(m_vResults.size());
	for(size_t i = 0;i < m_vResults.size(); ++i)
	{
		if(m_vResults[i].m_vPeptideResults.size() <= 0)
			continue;
		const CXLinkPepResult & pep_res = m_vResults[i].m_vPeptideResults[0];
		if(!m_vBeOutput[i])
		{
			m_vFeatureList.push_back(pair<double,size_t>(m_vSpectra[i].m_nCharge,i));
		}
	}
	sort(m_vFeatureList.begin(),m_vFeatureList.end(),_FeatureLesser);
	
	return _CalFeature(bSet);
}


int CXLinkResultFDRFilter::_CalDeltaMass(bool bSet )
{
	m_vFeatureList.clear();
	m_vFeatureList.reserve(m_vResults.size());
	for(size_t i = 0;i < m_vResults.size(); ++i)
	{
		if(m_vResults[i].m_vPeptideResults.size() <= 0)
			continue;
		const CXLinkPepResult & pep_res = m_vResults[i].m_vPeptideResults[0];
		if(!m_vBeOutput[i])
		{
			m_vFeatureList.push_back(pair<double,size_t>(m_vSpectra[i].m_lfMH-m_vResults[i].m_vPeptideResults[0].m_lfCalc_MH,i));
		}
	}
	sort(m_vFeatureList.begin(),m_vFeatureList.end(),_FeatureLesser);
	
	return _CalFeature(bSet);
}

int CXLinkResultFDRFilter::_CalIntenRatio(bool bSet )
{
	m_vFeatureList.clear();
	m_vFeatureList.reserve(m_vResults.size());
	for(size_t i = 0;i < m_vResults.size(); ++i)
	{
		if(m_vResults[i].m_vPeptideResults.size() <= 0)
			continue;
		const CXLinkPepResult & pep_res = m_vResults[i].m_vPeptideResults[0];
		if(!m_vBeOutput[i])
		{
			m_vFeatureList.push_back(pair<double,size_t>(pep_res.m_stMatchInfo.lfMatchedSpecInt/(pep_res.m_stMatchInfo.lfUnMatchedSpecInt+1.0E-5),i));
		}
	}
	sort(m_vFeatureList.begin(),m_vFeatureList.end(),_FeatureLesser);
	
	return _CalFeature(bSet);
}

int CXLinkResultFDRFilter::_CalTagLength(bool bSet )
{
	m_vFeatureList.clear();
	m_vFeatureList.reserve(m_vResults.size());
	for(size_t i = 0;i < m_vResults.size(); ++i)
	{
		if(m_vResults[i].m_vPeptideResults.size() <= 0)
			continue;
		const CXLinkPepResult & pep_res = m_vResults[i].m_vPeptideResults[0];
		if(!m_vBeOutput[i])
		{
			m_vFeatureList.push_back(pair<double,size_t>(_GetTagLength(pep_res),i));
		}
	}
	
	sort(m_vFeatureList.begin(),m_vFeatureList.end(),_FeatureLesser);
	
	return _CalFeature(bSet);
}

int CXLinkResultFDRFilter::_CalScore(bool bSet )
{
	m_vFeatureList.clear();
	m_vFeatureList.reserve(m_vResults.size());
	for(size_t i = 0;i < m_vResults.size(); ++i)
	{
		if(m_vResults[i].m_vPeptideResults.size() <= 0)
			continue;
		const CXLinkPepResult & pep_res = m_vResults[i].m_vPeptideResults[0];
		if(!m_vBeOutput[i])
			m_vFeatureList.push_back(pair<double,size_t>(pep_res.m_lfScore,i));
	}
	
	sort(m_vFeatureList.begin(),m_vFeatureList.end(),_FeatureGreater);
	
	return _CalFeature(bSet);
}

int CXLinkResultFDRFilter::_CalFeature(bool bSet)
{
	int nOutCount = 0;
	double lfCurFDR = 0.0;
	int nTCount = 0,nFCount = 0,nUCount =0;
	int nWindowCount = int(m_vFeatureList.size()/MIN_WINDOW);
	
	if(nWindowCount == 0)
	{
		nTCount = 0;
		nFCount = 0;
		nUCount =0;
		for(size_t j = 0 ; j < m_vFeatureList.size() ; ++ j)
		{
			if(m_vTUFTag[m_vFeatureList[j].second] == 'T')
			{
				nTCount ++;
			}
			else if(m_vTUFTag[m_vFeatureList[j].second] == 'U')
			{
				nUCount ++;
			}
			else if(m_vTUFTag[m_vFeatureList[j].second] == 'F')
				nFCount ++;
		}
		
		int nFPCount = 0;
		if(nUCount - 2* nFCount < 0)
			nFPCount = nFCount ;
		else
			nFPCount = nUCount - nFCount ;
		
		lfCurFDR = nTCount != 0 ? ((nFPCount + 0.0 )/nTCount): 1.0;
		
		if(bSet)
		{
			cout << "[" << m_vFeatureList[0].first << "," << m_vFeatureList[m_vFeatureList.size()-1].first << "]:" << endl;
			cout << "lfCurFDR = " << lfCurFDR << endl
			<< "T =" << nTCount << endl
			<< "U =" << nUCount << endl
			<< "F =" << nFCount << endl;
		}
		if(lfCurFDR <= m_filerconf.m_lfFDR)
		{
			if(bSet)
			{
				cout << "yes.." << endl;
				for(size_t j = 0; j < m_vFeatureList.size(); ++ j )
				{
					if(m_vTUFTag[m_vFeatureList[j].second] == 'T')
						m_vBeOutput[m_vFeatureList[j].second] = 1;
					else
						m_vBeOutput[m_vFeatureList[j].second] = -1;
				}
			}
			nOutCount += (nTCount);
			
			char tmpch;
			//cin >> tmpch;
		}
		else if(lfCurFDR >= 1 - SPECIFICITY)
		{
			if(bSet)
			{
				cout << "no.." << endl;
				for(size_t j = 0 ; j < m_vFeatureList.size(); ++ j)
				{
					m_vBeOutput[m_vFeatureList[j].second] = -1;
				}
			}
		}
		return nOutCount;
	}
	
	
	for(int i = 0 ; i < nWindowCount ; ++ i)
	{
		int nBegin = i * MIN_WINDOW;
		int nEnd = nBegin + MIN_WINDOW;
		if(m_vFeatureList.size() - nEnd < MIN_WINDOW)
			nEnd = m_vFeatureList.size();
		nTCount = 0;
		nFCount = 0;
		nUCount =0;
		for(int j = nBegin ; j < nEnd ; ++ j)
		{
			if(m_vTUFTag[m_vFeatureList[j].second] == 'T')
			{
				nTCount ++;
			}
			else if(m_vTUFTag[m_vFeatureList[j].second] == 'U')
			{
				nUCount ++;
			}
			else if(m_vTUFTag[m_vFeatureList[j].second] == 'F')
				nFCount ++;
		}
		
		int nFPCount = 0;
		if(nUCount - 2* nFCount < 0)
			nFPCount = nFCount ;
		else
			nFPCount = nUCount - nFCount ;
		
		lfCurFDR = nTCount != 0 ? ((nFPCount + 0.0 )/nTCount): 1.0;
		
		if(bSet)
		{
			cout << "[" << m_vFeatureList[nBegin].first << "," << m_vFeatureList[nEnd].first << "]:" << endl;
			cout << "lfCurFDR = " << lfCurFDR << endl
			<< "T =" << nTCount << endl
			<< "U =" << nUCount << endl
			<< "F =" << nFCount << endl;
			char tmpch;
			//cin >> tmpch;
		}
		if(lfCurFDR <= m_filerconf.m_lfFDR)
		{
			if(bSet)
			{
				cout << "yes.." << endl;
				for(int k = nBegin ; k < nEnd ; ++ k)
				{
					if(m_vTUFTag[m_vFeatureList[k].second] == 'T')
						m_vBeOutput[m_vFeatureList[k].second] = 1;
					else
						m_vBeOutput[m_vFeatureList[k].second] = -1;
				}
			}
			nOutCount += (nTCount);
			
			char tmpch;
			//cin >> tmpch;
		}
		else if(lfCurFDR >= 1 - SPECIFICITY)
		{
			if(bSet)
			{
				cout << "no.." << endl;
				for(int k = nBegin ; k < nEnd ; ++ k)
				{
					m_vBeOutput[m_vFeatureList[k].second] = -1;
				}
			}
		}
	}
	
	return nOutCount;
}


int CXLinkResultFDRFilter::_CalEvalue(bool bSet )
{
	m_vFeatureList.clear();
	m_vFeatureList.reserve(m_vResults.size());
	for(size_t i = 0;i < m_vResults.size(); ++i)
	{
		if(m_vResults[i].m_vPeptideResults.size() <= 0)
			continue;
		const CXLinkPepResult & pep_res = m_vResults[i].m_vPeptideResults[0];
		if(!m_vBeOutput[i])
		{
			m_vFeatureList.push_back(pair<double,size_t>(pep_res.m_lfEvalue,i));
		}
	}
	
	sort(m_vFeatureList.begin(),m_vFeatureList.end(),_FeatureLesser);
	
	return _CalFeature(bSet);

}

void CXLinkResultFDRFilter::_FilterByFDR()
{
	int nFound = 0;
	int nWho = 0;
	while(nFound < m_nTargetNum)
	{
		// for evalue
		int nMaxFound = 0;
		cout << "evalue :" << endl;
		int nCurFound = _CalEvalue();
		if(nCurFound > nMaxFound)
		{
			nMaxFound = nCurFound;
			nWho = 0;
		}
		// for score
		cout << "score :" << endl;
		nCurFound = _CalScore();
		if(nCurFound > nMaxFound)
		{
			nMaxFound = nCurFound;
			nWho = 1;
		}
		// for tag length
		cout << "tag length" << endl;
		nCurFound = _CalTagLength();
		if(nCurFound > nMaxFound)
		{
			nMaxFound = nCurFound;
			nWho = 2;
		}
		// for intensity
		cout << "intensity " << endl;
		nCurFound = _CalIntenRatio();
		if(nCurFound > nMaxFound)
		{
			nMaxFound = nCurFound;
			nWho = 3;
		}
		// for delta mass
		cout << "delta mass" << endl;
		nCurFound = _CalDeltaMass();
		if(nCurFound > nMaxFound)
		{
			nMaxFound = nCurFound;
			nWho = 4;
		}
		// for charge
		cout << "charge" << endl;
		nCurFound = _CalCharge();
		if(nCurFound > nMaxFound)
		{
			nMaxFound = nCurFound;
			nWho = 5;
		}
		
		if(nMaxFound <= 0)
		{
			// can't find any more
			break;
		}
		cout << "Choose " ;
		
		switch(nWho)
		{
		case 0:
			cout << "evalue .." << endl
			<< "Found " << nMaxFound;
			_CalEvalue(true);
			break;
		case 1:
			cout << "score .." << endl
			<< "Found " << nMaxFound;
			_CalScore(true);
			break;
		case 2:
			cout << "tag length .." << endl
			<< "Found " << nMaxFound;
			_CalTagLength(true);
			break;
		case 3:
			cout << "intensity .." << endl
			<< "Found " << nMaxFound;
			_CalIntenRatio(true);
			break;
		case 4:
			cout << "delta mass .." << endl
			<< "Found " << nMaxFound;
			_CalDeltaMass(true);
			break;
		case 5:
			cout << "charge .." << endl
			<< "Found " << nMaxFound;
			_CalCharge(true);
			break;
		default:
			break;
		}
		
		nFound += nMaxFound;
	}
	
	cout << "TargetNum = " << m_nTargetNum << endl
		<< "Found = " << nFound << endl;
}

bool CXLinkResultFDRFilter::_IsInTolWindow(double lfExpMH,double lfCalMH)
{
	double lfTolBase;
	double lfTol;
	if(m_filerconf.m_strTolBaseType == "ppm")
	{
		lfTolBase = lfCalMH * 0.000001 * m_filerconf.m_vTolBase[0];
	}
	else
		lfTolBase = m_filerconf.m_vTolBase[0];
	
	if(m_filerconf.m_strTolType == "ppm")
	{
		lfTol = lfCalMH * 0.000001 * m_filerconf.m_lfTol;
	}
	else
		lfTol = m_filerconf.m_lfTol;
	
	double lfTolLB ,lfTolUB;	
	lfTolLB = lfCalMH + lfTolBase - lfTol;
	lfTolUB = lfCalMH + lfTolBase + lfTol;
	
	if(lfExpMH <= lfTolUB && lfExpMH >= lfTolLB)
		return true;
	else
		return false;
	
}

void CXLinkResultFDRFilter::Run()
{
	CXLinkResultReportFactory reportFactory;
	CXLinkResultReport * report = reportFactory.GetReport(Xlink_STANDARD);
	if(report==NULL)
	{
		cout << "can't get the report" << endl;
		return ;
	}
	report->Init(m_filerconf.m_strpFindFile);

	cout << "loading from input protein report.." << endl;
	bool bOK = report->LoadFile(m_filerconf.m_strInputFile,m_vResults,m_vSpectra);
	if(bOK)
	{
		cout << "loading complete .." << endl
			<< "result num = " << m_vResults.size() << endl;
	}
	else
	{
		cout << "errors occur while loading .. " << endl;
		cout << m_filerconf.m_strInputFile<<endl;
		report->Close();
		return;
	}
	size_t tResultSize = 0;
	tResultSize = m_vResults.size();
	
	m_vBeOutput.clear();
	vector<int>().swap(m_vBeOutput);
	for(size_t i = 0 ;i < tResultSize; ++ i)
	{
		if(m_vResults[i].m_vPeptideResults.size()<=0 || m_vResults[i].m_vPeptideResults[0].m_XLink.m_eXLinkType != m_filerconf.m_nXLinkType)
		{
			
			m_vBeOutput.push_back(-1);
		}
		else
		{
			// filter by mass range
			if(!_IsInTolWindow(m_vSpectra[i].m_lfMH,m_vResults[i].m_vPeptideResults[0].m_lfCalc_MH))
			{
				m_vBeOutput.push_back(-1);
				continue;
			}

			if(false == m_pepResultFilter.Filter(m_vResults[i].m_vPeptideResults[0]))
			{
				m_vBeOutput.push_back(-1);
			}
			else
			{
				m_vBeOutput.push_back(0);
			}
		}
	}

	_SetTUFTag();
	_EstimateTargetNum();
	_FilterByFDR();
	
	for(size_t i = 0 ;i < tResultSize; ++ i)
	{
		if(m_vResults[i].m_vPeptideResults.size()<=0)
		{
			continue;
		}
		
		if(m_vBeOutput[i] != 1)
		{
			m_vResults[i].m_vPeptideResults.clear();
			vector<CXLinkPepResult>().swap(m_vResults[i].m_vPeptideResults);
		}
	}
	
	
	cout << "saving to output protein report ..." << endl;
	string strTitle = "_" + m_filerconf.m_strTitle;
	report->WriteFile(m_vResults,m_vSpectra,m_filerconf.m_strOutputPath,strTitle);
	report->Close();
	cout << "saving complete .. " << endl;
	
	m_vResults.clear();
	vector<CXLinkMatchResult>().swap(m_vResults);
	m_vSpectra.clear();
	vector<CSpectrum>().swap(m_vSpectra);
	
	delete report;
	
	return ;
}





