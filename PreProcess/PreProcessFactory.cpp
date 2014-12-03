#include "../include/sdk.h"
#include "../include/predefine.h"
#include "../include/interface.h"
#include "XLinkPreProc.h"
#include "PreProcessFactory.h"

using namespace proteomics_sdk;

CPreProcessFactory::CPreProcessFactory(void)
{
}

CPreProcessFactory::~CPreProcessFactory(void)
{
}

CPreProcess * CPreProcessFactory::GetPreProcessor(PreProcType ePP) const
{
	try
	{
	if(PRE_PROC_XLINK_HCD == ePP)
	{		
		return new CXLinkPreProc(); 
	}
	else
	{
		return new CXLinkPreProc();
	}	
	}
	catch(runtime_error & e)
	{
        CErrInfo info("CPreProcessFactory", "GetPreProcessor");

        throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
        CErrInfo info("CPreProcessFactory", "GetPreProcessor", "caught an unkown exception.");

        throw runtime_error(info.Get().c_str());
	}
}
