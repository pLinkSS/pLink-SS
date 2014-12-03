#ifndef XLINK_ACINFER_H_
#define XLINK_ACINFER_H_


//using namespace std;
using namespace proteomics_sdk;
namespace proteomics_sdk
{
class CXLinkMatchResult;
class CAssignedProtein;
class CXLinkACInfer
{
	
public:
	CXLinkACInfer(void):m_pCond(NULL){
		m_pTrace = CTrace::GetInstance();
		m_pTrace->Debug("CXLinkACInfer::CXLinkACInfer()");
	}
	virtual ~CXLinkACInfer(void){};
	virtual vector<CAssignedProtein> & Infer(vector<CXLinkMatchResult> & vResults, 
			const vector<CSpectrum> & vSpectra, FILTER_CRITERIA_INFO & cs);
	virtual void SetSign(string strFPSign){m_strFPSign = strFPSign;};
	virtual void Initialize(CCondition * pCond) {
		m_pTrace->Debug("in Initialize");
		m_pCond = pCond;
	}
	
	
protected:
	vector<CAssignedProtein>  m_vFoundProteins;
	string m_strFPSign;
	CCondition * m_pCond;
	CTrace *m_pTrace;
	double Filter(const vector<CXLinkMatchResult> & vResults, FILTER_CRITERIA_INFO & cs);
	static bool Less(const CAssignedProtein & a, const CAssignedProtein & b);
	static bool PeptideLess(const pair<const CXLinkPepResult *, size_t> & a, const pair<const CXLinkPepResult *, size_t> & b);
	
};
}

#endif /*XLINK_ACINFER_H_*/
