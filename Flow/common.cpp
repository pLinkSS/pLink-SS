#include "../include/sdk.h"
#include "../include/interface.h"
#include "../include/option.h"
#include "common.h"
using namespace proteomics_sdk;
bool Score_Less(const CPeptideResult & pep1, const CPeptideResult & pep2)
{
	return pep1.m_lfScore < pep2.m_lfScore;
}

bool XLINK_Score_Less(const CXLinkPepResult & pep1, const CXLinkPepResult & pep2)
{
	if(pep1.m_lfScore < pep2.m_lfScore)
	{
		return true;
	}
	else if(pep1.m_lfScore > pep2.m_lfScore)
	{
		return false;
	}
	else
	{
		char szbuf[1024];
		sprintf(szbuf,"%s(%d)-%s(%d)",
			pep1.m_AlphaPeptide.m_szSequence,
			pep1.m_XLink.m_tAlphaSite,
			pep1.m_BetaPeptide.m_szSequence,
			pep1.m_XLink.m_tBetaSite);
		string strPep1 = szbuf;
		sprintf(szbuf,"%s(%d)-%s(%d)",
			pep2.m_AlphaPeptide.m_szSequence,
			pep2.m_XLink.m_tAlphaSite,
			pep2.m_BetaPeptide.m_szSequence,
			pep2.m_XLink.m_tBetaSite);
		string strPep2 = szbuf;
		return strPep1 < strPep2 ;
	}
	
}
bool XLINK_Score_More(const CXLinkPepResult & pep1, const CXLinkPepResult & pep2)
{
	return pep1.m_lfScore > pep2.m_lfScore;
}
bool MIN_MASS_LESS(const SPEC_SORTER_INDEX_INFO & a, const SPEC_SORTER_INDEX_INFO & b)
{
	return a.lfMin < b.lfMin;
}
bool XLINK_MIN_MASS_LESS(const SPEC_WND_INFO & a, const SPEC_WND_INFO & b)
{
	return a.lfMassMax < b.lfMassMax;
}
bool XLINK_OPENFLOW_Score_Less(const CXLinkOpenPepResult & pep1, const CXLinkOpenPepResult & pep2)
{
	return pep1.m_lfScore < pep2.m_lfScore;
}
bool XLINK_OPENFLOW_Score_More(const CXLinkOpenPepResult & pep1, const CXLinkOpenPepResult & pep2)
{
	return pep1.m_lfScore > pep2.m_lfScore;
}




