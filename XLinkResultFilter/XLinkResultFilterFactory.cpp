#include "XLinkResultFilterInterface.h"
#include "XLinkResultFilter.h"
#include "XLinkResultFDRFilter.h"
#include "XLinkResultFilterFactory.h"

CXLinkResultFilterFactory::CXLinkResultFilterFactory()
{
}

CXLinkResultFilterFactory::~CXLinkResultFilterFactory()
{
}

CXLinkResultFilterInterface * CXLinkResultFilterFactory::GetFilter(int nType)
{
	if(nType == 0)
	{
		return new CXLinkResultFilter();
	}
	else if(nType == 1)
	{
		return new CXLinkResultFDRFilter();
	}
	else
		return new CXLinkResultFilter();
	
}
