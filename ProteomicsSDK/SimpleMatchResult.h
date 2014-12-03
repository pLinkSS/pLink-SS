#ifndef SIMPLEMATCHRESULT_H_
#define SIMPLEMATCHRESULT_H_

namespace proteomics_sdk
{

class CSimpleMatchResult
{
public:
	CSimpleMatchResult();

	void remove_invalid(void);
	
	static double Calc_Theoretical_MH(const CPeptideResult & pep_res, bool bPepMono);
	
	void track(const CSpectrum & spectrum, const CPeptideResult & pep_res, size_t tCurrentPro, CCondition & cond);
	
	vector<double> m_vlfScores;
	
	vector<CPeptideResult> m_vPeptideResults;
	
	size_t m_tScore;

	size_t m_tRealCandidate;
};

}

#endif /*SIMPLEMATCHRESULT_H_*/
