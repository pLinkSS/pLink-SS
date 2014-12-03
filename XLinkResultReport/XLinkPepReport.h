#ifndef XLINKPEPREPORT_H_
#define XLINKPEPREPORT_H_

class CXLinkPepReport : public CXLinkResultReport
{
public:
	CXLinkPepReport();
	virtual ~CXLinkPepReport();
	virtual void Init(string strpFindFile,time_t tmStartTime = 0);
	virtual void GetLines(const vector<CXLinkMatchResult> & vResults , const vector<CSpectrum> & vSpectra ,string & strTXT);
	virtual bool LoadFile(string strReportFile, vector<CXLinkMatchResult> & vResults , vector<CSpectrum> & vSpectra );
	virtual void WriteFile(const vector<CXLinkMatchResult> & vResults , const vector<CSpectrum> & vSpectra , string strOutputPath , string strTitle = "");
	virtual void Close(void);
	
protected:
	
	string _getPepString(const CXLinkPepResult & pep_res);

	string m_strpFindFile;
	string m_strWorkDir;
	CCondition m_Condition;
	time_t m_tmStartTime;
};

#endif /*XLINKPEPREPORT_H_*/
