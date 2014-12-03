#include <iostream>
#include <string>
 
#include "../include/sdk.h"
#include "../include/interface.h"
#include "../include/option.h"
#include "Inferor.h"
#include "InferorFactory.h"

using namespace std;
using namespace proteomics_sdk;
int main(int argc, char * argv[])
{

	CTrace * pTrace = CTrace::GetInstance();
	
	string strArg("Preparing to infer proteins from peptides, the parameters are :\n");
	for(int i = 0;i < argc;++i)
		strArg += string("\t")+argv[i];
	strArg+=string("\n");
	pTrace->Info(strArg);

	if(argc < 6)
	{
		CErrInfo info("Inferor", "Main", "The parameters are invalid!");		
		ofstream ofs("Searcher.err.log", ios::app);
		ofs << info;	
		return 1;
	}
	
	string strpFindFile(argv[1]);
	string strOutputFile(argv[2]);
	int nFileTotal = atoi(argv[3]);
	int nSpectraTotal = atoi(argv[4]);
	string strIdentifier(argv[5]);
	
	// get time started
	time_t tmStartTime;
	if(argc == 7)
	{
		tmStartTime = atol(argv[6]); 
	}
	else
	{
		tmStartTime = time(NULL);
	}
	
	// judge whether to use xlink inferor ; added at 2013.6.18, judge whether use tri inferor
	INFEROR_TYPE eInferorType = COMMON_INFEROR;
	COptionTool * pOption = new COptionTool("flow", strpFindFile.c_str());
	COptionTool * pOption2 = new COptionTool("database", strpFindFile.c_str());
	if(pOption && pOption2)
	{
		string strFlow = pOption->GetString("flow", "FLOW_ALL_IDX");

		string strFlow2 = pOption2->GetString("index_content");//todo: here can't get, get NULL

		if("FLOW_XLINK" == strFlow)
		{
			if ("PEPTIDE_TRI_ALL" == strFlow2)
			{
				eInferorType = TRI_INFEROR;
			}
			else if ( "PEPTIDE_TRI_ION" == strFlow2 )
				eInferorType = TRI_INFEROR;
			else
				eInferorType = XLINK_INFEROR;
		}

		delete pOption;
		delete pOption2;
	}
	else
	{
		cout<<"error while read parameter."<<endl;
		return 0;
	}
	
	// todo: write pfd2plabel config file
	
	CInferorFactory inferorfactory;
	CInferor * pInferor = inferorfactory.GetInferor(eInferorType);
		
	try
	{
		if(!pInferor->Init(strpFindFile,tmStartTime))
		{
			pTrace->Alert("Initialization inferor module fail!");
			CErrInfo info("Inferor", "Main", "Init failed");
			throw invalid_argument(info.Get().c_str());
		}
	}
	catch(exception & e)
	{
		CErrInfo info("Inferor", "Main", "in the function Inferor::Init.");
		ofstream ofs("Searcher.err.log", ios::app);
		ofs << info;
		return 1;
	}
	catch(...)
	{
		CErrInfo info("Inferor", "Main", "caught an unknown exception in the function Inferor::Init.");
		ofstream ofs("Searcher.err.log", ios::app);
		ofs << info;
		return 1;
	}

	pTrace->Info("Complete initialization of inferor module.", MODULE_INFEROR);
	pTrace->Info("Start infering proteins from peptides ...", MODULE_INFEROR);
	
	try
	{
		pTrace->Debug("before pInferor->Run");
		pInferor->Run(strpFindFile, strOutputFile, nFileTotal, nSpectraTotal, strIdentifier);
		pTrace->Debug("after pInferor->Run");
	}
	catch(exception & e)
	{
		CErrInfo info("Inferor", "Main", "in the function Inferor::Run.");
		info.Append("strpFindFile=" + strpFindFile);
		info.Append("strOutputFile=" + strOutputFile);
		char temp[200] = {0};
		sprintf(temp, "nFileTotal=%d", nFileTotal);
		info.Append(temp);
		sprintf(temp, "nSpectraTotal=%d", nSpectraTotal);
		info.Append(temp);
		info.Append("strIdentifier=" + strIdentifier);
		ofstream ofs("Searcher.err.log", ios::app);
		ofs << info;
		return 1;
	}
	catch(...)
	{
		CErrInfo info("Inferor", "Main", "caught an unknown exception in the function Inferor::Run.");
		info.Append("strpFindFile=" + strpFindFile);
		info.Append("strOutputFile=" + strOutputFile);
		char temp[200] = {0};
		sprintf(temp, "nFileTotal=%d", nFileTotal);
		info.Append(temp);
		sprintf(temp, "nSpectraTotal=%d", nSpectraTotal);
		info.Append(temp);
		info.Append("strIdentifier=" + strIdentifier);
		ofstream ofs("Searcher.err.log", ios::app);
		ofs << info;
		return 1;
	}
	
	delete pInferor;
	
	pTrace->Info("Complete searching successfully.", MODULE_INFEROR);
	
	
	CTrace::DeleteInstance();

	return 0;
}
