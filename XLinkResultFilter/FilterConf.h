#ifndef CFILTERCONF_H_
#define CFILTERCONF_H_

/*
XLinkType=3
PepTolType=ppm
PepTol=50
PepTolBaseType=Da
PepTolBase=1
MaxEvalue=1.000000
MinScore=0.000000
MinPepLength=0
MinContiAAnum=0
AAConfLevel=1
MinMatchOverUnMatch=0.000000
SaveCTerm=0
FDR=0.01
ReverseTag=REVERSE_
*/

class CFilterConf
{
public:
	CFilterConf();
	virtual ~CFilterConf();
	void Load(string strConfigFile);
	void Output();
	
	string m_strWorkDir;
	string m_strOption;
	
	string m_strpFindFile;
	string m_strInputFile;
	string m_strOutputPath;
	string m_strTitle;
	
	double m_lfTol;
	string m_strTolType;
	vector<double> m_vTolBase;
	string m_strTolBaseType;
	
	int m_nXLinkType;
	int m_nLinkerId;
	
	double m_lfFDR;
	string m_strReverseTag;
	
	string m_strInclusionList;
	bool m_bInclusionListAvail;
	
	// condition for pfind
	CCondition m_Condition;
	
	// filter for pep context
	XLinkFilterStandand m_pepFilter;
	
protected:
	void _LoadIO();
	void _LoadFilter();
};

#endif /*CFILTERCONF_H_*/
