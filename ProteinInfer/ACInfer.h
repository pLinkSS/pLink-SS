#ifndef ACINFER_H_
#define ACINFER_H_


//using namespace std;
using namespace proteomics_sdk;
namespace proteomics_sdk
{
class CSimpleMatchResult;
class CAssignedProtein;
class CACInfer:public CProteinInfer
{
	
public:
	CACInfer(void):m_pCond(NULL){;};
	virtual ~CACInfer(void){};
	virtual vector<CAssignedProtein> & Infer(vector<CSimpleMatchResult> & vResults, 
			const vector<CSpectrum> & vSpectra, FILTER_CRITERIA_INFO & cs);
	virtual ProteinInferType GetType(void){return PROTEIN_INFER_AC;};
	virtual void SetSign(string strFPSign){m_strFPSign = strFPSign;};
	virtual void Initialize(CCondition * pCond){m_pCond = pCond;}

protected:
	vector<CAssignedProtein>  m_vFoundProteins;
	string m_strFPSign;
	CCondition * m_pCond;
	double Filter(const vector<CSimpleMatchResult> & vResults, FILTER_CRITERIA_INFO & cs);
	static bool Less(const CAssignedProtein & a, const CAssignedProtein & b);

	static int _GetRoot(int  root[], int nOrd);

	static bool _IsSame(int x, int y, vector<string> * vAssignPep);

	static bool PeptideLess(const pair<const CPeptideResult *, size_t> & a, const pair<const CPeptideResult *, size_t> & b);

	
};
}

#endif /*ACINFER_H_*/
