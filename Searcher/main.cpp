#include <fstream>
#include <time.h>
#include <exception>
#include <string>
#include <iostream>
#include <stdexcept>
#include <limits>
#include "../include/sdk.h"
#include "../include/option.h"
#include "../include/interface.h"
#include "../SearchEngine/SearchEngine.h"
#include "Searcher.h"
using namespace proteomics_search;
using namespace std;
using namespace proteomics_sdk;

int main(int argc, char * argv[])
{
	CSearcher s;
	
	try
	{
		s.Init(argc, &argv[0]);
		s.Run();
	}
	catch(exception & e)
	{
		CErrInfo info("Searcher.exe", "Main", "failed.", e);
		ofstream ofs("Searcher.err.log", ios::app);
		ofs << info;
		cout << info << endl;
		return 1;
	}
	catch(...)
	{
		CErrInfo info("Searcher.exe", "Main", "caught an unknown exception.");
		ofstream ofs("Searcher.err.log", ios::app);
		ofs << info;
		cout << info << endl;
		return 1;
	}
	
	return 0;
}
