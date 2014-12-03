#include <string>
#include <vector>
#include "../include/predefine.h"
#include "ProteomicsSDK.h"
#include "Condition.h"
#include "SimpleMatchResult.h"

using namespace std;
using namespace proteomics_sdk;

CSimpleMatchResult::CSimpleMatchResult():
		m_tScore(0),m_tRealCandidate(0)
{};

void CSimpleMatchResult::remove_invalid(void)
{
		vector<CPeptideResult> vTemp;
		size_t i = 0;
		for(i = 0;i < m_vPeptideResults.size();++i)
		{
			if(m_vPeptideResults[i].m_bEV)
			{
				vTemp.push_back(m_vPeptideResults[i]);
			}
		}
		vTemp.swap(m_vPeptideResults);
};
	
double CSimpleMatchResult::Calc_Theoretical_MH(const CPeptideResult & pep_res, bool bPepMono)
{
	if(bPepMono)
		return pep_res.m_peptide.m_lfMass + IonMass_Mono_O + IonMass_Mono_H * 2.0 + IonMass_Proton;
	else
		return pep_res.m_peptide.m_lfMass + IonMass_Aver_O + IonMass_Aver_H * 3.0;

};
	
void CSimpleMatchResult::track(const CSpectrum & spectrum, const CPeptideResult & pep_res, size_t tCurrentPro, CCondition & cond)
{
	
	double lfCalc_MH = 0.0;
	lfCalc_MH = CSimpleMatchResult::Calc_Theoretical_MH(pep_res,
			cond.m_bPepMono); 
	double lfMZ = cond.GetMZ(spectrum.m_lfMH, spectrum.m_nCharge);
	double lfPepTol = cond.GetPepTol(lfMZ, spectrum.m_nCharge);
	// modify by emily
	double lfPepTolBase = cond.GetPepTolBase(lfMZ, spectrum.m_nCharge);
	
	if(fabs(spectrum.m_lfMH + lfPepTolBase - lfCalc_MH) > lfPepTol)
	{
		if(m_vlfScores.size() <= cond.m_tMinScoreNum)
		    	m_vlfScores.push_back(pep_res.m_lfScore);
		++m_tScore;
		return;
	}
	
	++m_tRealCandidate;
	++m_tScore;
	
	if(m_vPeptideResults.size() <= cond.m_nReportPep)
	{
		m_vPeptideResults.push_back(pep_res);
		m_vPeptideResults.back().m_vProteinID.push_back(tCurrentPro);
		sort(m_vPeptideResults.begin(), m_vPeptideResults.end(), CPeptideResult::Score_Greater);
	}
	else
	{
		bool bRedun = false;
		for(size_t j = m_vPeptideResults.size() - 1;;--j)
		{
			if(fabs(pep_res.m_lfScore - m_vPeptideResults[j].m_lfScore) <= EPS)
			{
				if(m_vPeptideResults[j].m_peptide.EqualPeptide(pep_res.m_peptide))
				{
					if(find(m_vPeptideResults[j].m_vProteinID.begin(), 
							m_vPeptideResults[j].m_vProteinID.end(),
							tCurrentPro) == m_vPeptideResults[j].m_vProteinID.end())
					m_vPeptideResults[j].m_vProteinID.push_back(tCurrentPro);
					bRedun = true;
					break;
				}
				
			}
			else
				if(pep_res.m_lfScore < m_vPeptideResults[j].m_lfScore)
				{
					break;
				}
			if(0 == j)
				break;
		}
		if(!bRedun)
		{
			m_vPeptideResults[cond.m_nReportPep] = pep_res;
			m_vPeptideResults.back().m_vProteinID.push_back(tCurrentPro);
			sort(m_vPeptideResults.begin(), m_vPeptideResults.end(), CPeptideResult::Score_Greater);
		}
	}
    if(m_vlfScores.size() <= cond.m_tMinScoreNum)
    	m_vlfScores.push_back(pep_res.m_lfScore);
}
