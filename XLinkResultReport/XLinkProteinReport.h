#ifndef XLINKPROTEINREPORT_H_
#define XLINKPROTEINREPORT_H_

class CXLinkProteinReport : public CXLinkResultReport
{
public:
	CXLinkProteinReport();
	virtual ~CXLinkProteinReport();
	virtual void Init(string strpFindFile,time_t tmStartTime = 0);
	virtual void GetLines(const vector<CXLinkMatchResult> & vResults , const vector<CSpectrum> & vSpectra ,string & strTXT);
	virtual bool LoadFile(string strReportFile, vector<CXLinkMatchResult> & vResults , vector<CSpectrum> & vSpectra );
	virtual void WriteFile(const vector<CXLinkMatchResult> & vResults , const vector<CSpectrum> & vSpectra , string strOutputPath , string strTitle = "");
	virtual void Close(void);

protected:
	double _GetXLinkerMass(const CXLinkPepResult & pep_res);
	double _GetDeltaMass(double lfCalcMH, double lfExpMH);
	int _getModifyId(string strModName);
	void _parseline(string strLine,char cSep,vector<string > & vStrs);
	bool _getline(FILE * fp, string strTitle,string & strValue);

	CCondition m_Condition;
	CTrace *m_pTrace;
	time_t m_tmStartTime;
	string m_strpFindFile;
	string m_strWorkDir;
	string m_strSpectraPath;
};

#endif /*XLINKPROTEINREPORT_H_*/
