#include <iostream>
#include <string>
 
#include "../include/predefine.h"
#include "../include/sdk.h"
#include "../include/interface.h"
#include "../include/option.h"

#include "../XLinkResultReport/XLinkResultReport.h"
#include "../XLinkResultReport/XLinkResultReportFactory.h"

using namespace std;
using namespace proteomics_sdk;

bool LoadInput(string strOption);

string strpFindFile;
string strReportFile;
string strInputType;
string strOutputPath;
string strTitle;
string strFormat;

bool LoadInput(string strOption)
{
	COptionTool * pOption = new COptionTool("IO", strOption.c_str());
    if(!pOption)
    {
    	cout << "can't open pOption file : " << strOption << " on [IO]" << endl;
    	return false;
    }
    
    strpFindFile = pOption->GetString("pFindFile","");
    strReportFile = pOption->GetString("ReportFile","");
    strInputType = pOption->GetString("InputType","standard");
    strOutputPath = pOption->GetString("OutputPath","");
    strTitle = pOption->GetString("title","");
    strFormat = pOption->GetString("format","");
    delete pOption;
    return true;
	
}

int main(int argc , char * argv[])
{
	if(argc!=2)
	{
		cout << "usage : " << endl
		<< "XLinkResultConverter.exe converter.txt" << endl;
		return 0;
	}
	
	if(false == LoadInput(argv[1]))
	{
		cout << "Load Input File fail " << endl;
		return 0;
	}
	enum XLINK_RESULT_REPORT_TYPE reportType;

	if(strcmp(strFormat.c_str(),"standard") == 0)
	{
		reportType = Xlink_STANDARD;
	}
	else if (strcmp(strFormat.c_str(),"plabel") == 0)
	{
		reportType = XLINK_PLABEL;
	}
	else if (strcmp(strFormat.c_str(),"pbuild") == 0)
	{
		reportType = XLINK_PBUILD;
	}
	else if (strcmp(strFormat.c_str(),"pxbuild") == 0)
	{
		reportType = XLINK_PXBUILD;
	}
	else if (strcmp(strFormat.c_str(),"peptide") == 0)
	{
		reportType = XLINK_PEPTIDE;
	}
	else
		reportType = Xlink_STANDARD;
	
	CXLinkResultReportFactory reportFactory;
	
	//CXLinkResultReport * report = reportFactory.GetReport(XLINK_PBUILD);
	// todo : temp from pbuild -> standard output
	CXLinkResultReport * report = NULL;
	if(strInputType == "pbuild")
	{
		report = reportFactory.GetReport(XLINK_PBUILD);
	}
	else
	{
		report = reportFactory.GetReport(Xlink_STANDARD);
	}
	
	if(report == NULL)
	{
		cout << "can't get the report" << endl;
		return 0;
	}
	
	report->Init(strpFindFile);
	//cout << "report->init " << endl;
	vector<CXLinkMatchResult> vResults;
	vector<CSpectrum> vSpectra;

	//cout << "loading from input protein report.." << endl;
	bool bOK = report->LoadFile(strReportFile,vResults,vSpectra);
	
	if(bOK)
	{
		//cout << "loading complete .." << endl
			//<< "result num = " << vResults.size() << endl;
	}
	else
	{
		cout << "errors occur while loading .. " << endl;
		report->Close();
		return 0;
	}
	
	report->Close();
	delete report;
	
	//cout << "write to outputfile as the designated format .." << endl;
	report = reportFactory.GetReport(reportType);
	if(report == NULL)
	{
		cout << "can't get the report" << endl;
		return 0;
	}
	report->Init(strpFindFile);
	if(!strTitle.empty())
		strTitle = "_" + strTitle;
	report->WriteFile(vResults,vSpectra,strOutputPath,strTitle);
	report->Close();
	delete report;
	
}
