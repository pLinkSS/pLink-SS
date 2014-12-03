#include <string>
#include <exception>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include "../include/sdk.h"
#include "Indexer.h"

using namespace std;

int main (int argc, char *argv[])
{
	string strParamFile= "./task.pindex";
	if(argc > 1)
	{
		strParamFile = argv[1];
		for(long t = 2; t < argc; ++t)
		{
			strParamFile += " " + string(argv[t]);
		}
	}

	try
	{
		CExcuteIndexer indexer;	
		indexer.Init(::Protein_Disk, ::ICT_SPSingleScan, strParamFile);
		indexer.Run();
	}
	catch(exception & e)
	{
		CErrInfo info("Index.exe", "Main", "failed.", e);
		ofstream ofs("Index.err.log", ios::app);
		ofs << info;
		cout << info << endl;
		return 1;
	}
	catch(...)
	{
		CErrInfo info("Index.exe", "Main", "caught an unknown exception.");
		ofstream ofs("Index.err.log", ios::app);
		ofs << info;
		cout << info << endl;
		return 1;
	}	
	
	return 0;
}
