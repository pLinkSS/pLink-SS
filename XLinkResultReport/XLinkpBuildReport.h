#ifndef XLINKPBUILDREPORT_H_
#define XLINKPBUILDREPORT_H_

class CXLinkpBuildReport : public CXLinkResultReport
{
public:
	CXLinkpBuildReport();
	virtual ~CXLinkpBuildReport();
	virtual void Init(string strpFindFile,time_t tmStartTime = 0);
	virtual void GetLines(const vector<CXLinkMatchResult> & vResults , const vector<CSpectrum> & vSpectra ,string & strTXT);
	virtual bool LoadFile(string strReportFile, vector<CXLinkMatchResult> & vResults , vector<CSpectrum> & vSpectra );

	virtual void WriteFile(const vector<CXLinkMatchResult> & vResults , const vector<CSpectrum> & vSpectra , string strOutputPath , string strTitle = "");
	virtual void Close(void);
protected:
	double _GetDeltaMass(double lfCalcMH, double lfExpMH);
	double _GetXLinkerMass(const CXLinkPepResult & pep_res);
	bool _getline(FILE * fp, string strTitle,string & strValue);
	void _parseline(string strLine,char cSep,vector<string > & vStrs);
	int _getModifyId(string strModName);
	string _getPepString(const CXLinkPepResult & pep_res);
	void _SetPeptideInfo(CXLinkPepResult & pep_res,string & strPepSeq);
	void _SetProteinInfo(CXLinkPepResult & pep_res,string & strProSeq);
	void _SetProteinIDInfo(CXLinkPepResult & pep_res,string & strProSeq);
	
	string m_strpFindFile;
	string m_strWorkDir;
	string m_strSpectraPath;
	CCondition m_Condition;
	
	time_t m_tmStartTime;
};

#endif /*LINKPBUILDREPORT_H_*/
