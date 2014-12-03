#include <fstream>
#include <time.h>
#include <string>
#include <exception>
#include <stdexcept>
#include <iostream>
#include <limits>

#ifdef WIN32
#include <direct.h>
#endif

#include "../include/sdk.h"
#include "../include/option.h"
#include "../include/interface.h"
#include "../SearchEngine/SearchEngine.h"
#include "../SpectraIO2.0/MS2Input/MS2InputFactory.h"

#include <pthread.h>
#include "Searcher.h"

using namespace std;
using namespace proteomics_sdk;
using namespace proteomics_search;
bool bAllBusy = false;

#ifdef WIN32
const char cSlash = '\\';
#else
const char cSlash = '/';
#endif
CSearcher::CSearcher()
{
}

CSearcher::~CSearcher()
{
    Close();
}

/* parse the command line 
 * input: command line parameters from the main function
 * output: true iff one paramter or the second parameter is like "xxxx.pfind"
 */
bool CSearcher::_ParseCmdLine(int nArgc, char ** pArgv)
{
	if(nArgc < 2)
	{
		m_strParamFile = "task.pfind";
		return true;
	}
	else
	{
		/*
		 * the valid command line is as "Searcher xxxx.pfind" in Windows
		 * or "./Searcher xxxx.pfind in Linux
		 */
		
		string str;
		m_strParamFile = pArgv[1];
		str = m_strParamFile.substr(m_strParamFile.length() - 6);
	    return (str == ".pfind" || str == ".PFIND");
	}
}

void CSearcher::Init(int nArgc, char ** pArgv)
{    	
	if(!_ParseCmdLine(nArgc, pArgv))
	{
		CErrInfo info("CSearcher", "Init", "parse the command line failed.");
		throw runtime_error(info.Get().c_str());
	}

    m_SearchEngine.Init(m_strParamFile);
}

void CSearcher::Run(void)
{   
	m_SearchEngine.Search();
	m_SearchEngine.Infer();
}

void CSearcher::Close(void)
{
}

