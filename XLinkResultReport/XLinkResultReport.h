#ifndef XLINKRESULTREPORT_H_
#define XLINKRESULTREPORT_H_

enum XLINK_RESULT_REPORT_TYPE
{
	Xlink_STANDARD = 0,
	XLINK_PBUILD = 1,
	XLINK_PLABEL = 2,
	XLINK_PXBUILD = 3,
	XLINK_PQUANT = 4,
	XLINK_PEPTIDE =5,
};


class CXLinkResultReport
{
public:
	CXLinkResultReport(){} ;
	virtual ~CXLinkResultReport(){};
	virtual void Init(string strpFindFile,time_t tmStartTime = 0) = 0;
	virtual void GetLines(const vector<CXLinkMatchResult> & vResults , const vector<CSpectrum> & vSpectra ,string & strTXT) = 0;
	virtual bool LoadFile(string strReportFile, vector<CXLinkMatchResult> & vResults , vector<CSpectrum> & vSpectra ) = 0;
	virtual void WriteFile(const vector<CXLinkMatchResult> & vResults , const vector<CSpectrum> & vSpectra , string strOutputPath , string strTitle = "") = 0;
	virtual void Close(void) = 0;
};

#endif /*XLINKRESULTREPORT_H_*/
