#include "../include/sdk.h"
#include "XLinkPepResultFilter.h"

using namespace proteomics_sdk;

CXLinkPepResultFilter::CXLinkPepResultFilter()
{
}

CXLinkPepResultFilter::~CXLinkPepResultFilter()
{
}


void CXLinkPepResultFilter::SetFilterStandand(XLinkFilterStandand stFilterStandand)
{
	m_stFilterStandand = stFilterStandand;
}

bool CXLinkPepResultFilter::Filter(const CXLinkPepResult & pep_res)
{
	
	if(!m_stFilterStandand.bCterm)
	{
		switch(pep_res.m_XLink.m_eXLinkType)
		case 1:
		{
			if(pep_res.m_XLink.m_tAlphaSite == pep_res.m_AlphaPeptide.m_tLength - 1)
				return false;
			else
				break;
		case 2:	
			if((pep_res.m_XLink.m_tAlphaSite == pep_res.m_AlphaPeptide.m_tLength - 1)||(pep_res.m_XLink.m_tBetaSite == pep_res.m_AlphaPeptide.m_tLength - 1))
				return false;
			else
				break;
		case 3:
			if((pep_res.m_XLink.m_tAlphaSite == pep_res.m_AlphaPeptide.m_tLength - 1)||(pep_res.m_XLink.m_tBetaSite == pep_res.m_BetaPeptide.m_tLength - 1))
				return false;
			else
				break;
		}
	}
	
	if(pep_res.m_lfEvalue > m_stFilterStandand.lfEvalue)
		return false;
	
	if(pep_res.m_lfScore < m_stFilterStandand.lfScore)
		return false;
	
	if(pep_res.m_AlphaPeptide.m_tLength < m_stFilterStandand.nPepLength)
		return false;
	if(pep_res.m_bPair)
	{
		
		if(pep_res.m_BetaPeptide.m_tLength < m_stFilterStandand.nPepLength)
			return false;
	}
	
	if(pep_res.m_stMatchInfo.lfMatchedSpecInt/pep_res.m_stMatchInfo.lfUnMatchedSpecInt < m_stFilterStandand.lfMatchOverUnMatch)
			return false;
	
	int nMaxLen = 0 , nCurLen = 0;
	for(size_t i=0;i<pep_res.m_AlphaPeptide.m_tLength;++i)
	{
		
		if(pep_res.m_stMatchInfo.aPepConf[i] >=  m_stFilterStandand.nAAConfLevel)
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
	
	//cout << "conti pep1 : " << nMaxLen << endl;
	if(nMaxLen < m_stFilterStandand.nContiAAnum)
		return false;

	if(pep_res.m_bPair)
	{
		nMaxLen = 0;
		nCurLen = 0;
		for(size_t i=0;i<pep_res.m_BetaPeptide.m_tLength;++i)
		{
			
			if(pep_res.m_stMatchInfo.aPepConf[pep_res.m_AlphaPeptide.m_tLength + i] >=  m_stFilterStandand.nAAConfLevel)
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
		
		//cout << "conti pep2 : " << nMaxLen << endl;
		if(nMaxLen < m_stFilterStandand.nContiAAnum)
			return false;
		
	}
	
	return true;
}
