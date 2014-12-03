#ifndef XLINKPEPRESULTFILTER_H_
#define XLINKPEPRESULTFILTER_H_

namespace proteomics_sdk
{

struct XLinkFilterStandand
{
	// maximum of e-value
	double lfEvalue;
	// minimum score
	double lfScore;
	// minimum pep length 
	int nPepLength;
	// minimum continuous identified amino length 
	int nContiAAnum;
	// minimum confidence level for identified amino 
	int nAAConfLevel;
	// minimum ratio of intensity sum of matched over unmatched
	double lfMatchOverUnMatch;
	// whether retain C-terminal linked peptide(s)
	bool bCterm;
	// xlink type
	int nXLinkType;
	
	void output()
	{
		cout << "XLINK FILTER STANDARD FOR PEPTIDE :" << endl		
		<< "evalue : " << lfEvalue << endl
		<< "score : " << lfScore << endl
		<< "nPepLength : " << nPepLength << endl
		<< "nContiAAnum : " << nContiAAnum << endl
		<< "nAAConfLevel : " << nAAConfLevel << endl
		<< "lfMatchOverUnMatch : " << lfMatchOverUnMatch << endl
		<< "bCterm : " << bCterm << endl;
		 
	}
		
	XLinkFilterStandand()
	{
		lfEvalue = 1.0;
		
		lfScore = 0.0;
		
		nPepLength = 0;
		 
		nContiAAnum = 1;
		 
		nAAConfLevel = 1;
		
		lfMatchOverUnMatch = 0.0;
		
		// retain it as default
		bCterm = true;
		
	}
};

class CXLinkPepResultFilter
{
public:
	CXLinkPepResultFilter();
	void SetFilterStandand(XLinkFilterStandand stFilterStandand);
	bool Filter(const CXLinkPepResult & pep_res);
	virtual ~CXLinkPepResultFilter();
private:
	XLinkFilterStandand m_stFilterStandand; 
	
};

}

#endif /*XLINKPEPRESULTFILTER_H_*/
