#ifndef XLINKPLABELREPORT_H_
#define XLINKPLABELREPORT_H_

class CXLinkResultReport;

class CXLinkpLabelReport : public CXLinkResultReport
{
public:
	CXLinkpLabelReport();
	virtual ~CXLinkpLabelReport();
	virtual void Init(string strpFindFile,time_t tmStartTime = 0);
	virtual void GetLines(const vector<CXLinkMatchResult> & vResults , const vector<CSpectrum> & vSpectra ,string & strTXT);
	virtual bool LoadFile(string strReportFile, vector<CXLinkMatchResult> & vResults , vector<CSpectrum> & vSpectra );
	virtual void WriteFile(const vector<CXLinkMatchResult> & vResults , const vector<CSpectrum> & vSpectra , string strOutputPath , string strTitle = "");
	virtual void Close(void);
protected:
	string _2Uppercase(string str);
	
	CCondition m_Condition;
	string m_strWorkDir;
	string m_strpFindFile;
	string m_strSpectraPath;
	time_t m_tmStartTime;
};



#endif /*XLINKPLABELREPORT_H_*/
