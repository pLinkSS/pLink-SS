#include <string>
#include "../include/predefine.h"
#include "../include/sdk.h"
#include "../include/interface.h"
#include "SimpleAppender.h"
#include "AppenderFactory.h"
using namespace std;
namespace proteomics_sdk
{

CAppenderFactory::CAppenderFactory()
{
}

CAppenderFactory::~CAppenderFactory()
{
}
CAppender *CAppenderFactory::GetAppender(LogAppenderType eType)const
{
	CAppender * p(NULL);
	try
	{
		if(LOG_APPENDER_CMD == eType)
		{
			p = new CSimpleAppender();
		}
	}
	catch(runtime_error & e)
	{
    CErrInfo info("CAppenderFactory", "GetAppender");

    throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
	    CErrInfo info("CAppenderFactory", "GetAppender", "caught an unkown exception.");
	
	    throw runtime_error(info.Get().c_str());
	}
	
	return p;
};
}

