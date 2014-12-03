#ifndef _XLINK_PEPTIDE_EVALUATER_H_INCLUDED_
#define _XLINK_PEPTIDE_EVALUATER_H_INCLUDED_

using namespace proteomics_sdk;

class proteomics_sdk::CSpectrum;
class proteomics_sdk::CCondition;

namespace proteomics_sdk
{

class CXLinkPeptideEvaluater
{
public:
	CXLinkPeptideEvaluater(void);
	virtual ~CXLinkPeptideEvaluater(void);

	virtual void Init(const CCondition &condition);

	virtual void Run(CSpectrum & spectrum, CXLinkMatchResult * pResult);
	
	virtual void Close(void);

protected:

	void _FitHighScores(void);

	void _ComputeExpectation(const CSpectrum & spectrum);

	void _Isqt( double x[], double y[], size_t n);

	CXLinkMatchResult * m_pResult;

	CCondition m_Condition;

	double m_Scoefficient[2];
	int m_nEvalueNo;
	double m_DT[6];
};

}
#endif
